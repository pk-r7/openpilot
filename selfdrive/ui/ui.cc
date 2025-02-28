#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <assert.h>
#include <sys/mman.h>
#include <sys/resource.h>

#include <cutils/properties.h>

#include <GLES3/gl3.h>
#include <EGL/egl.h>

#include <json.h>
#include <czmq.h>

#include "nanovg.h"
#define NANOVG_GLES3_IMPLEMENTATION
#include "nanovg_gl.h"
#include "nanovg_gl_utils.h"

#include "common/messaging.h"
#include "common/timing.h"
#include "common/util.h"
#include "common/swaglog.h"
#include "common/mat.h"

extern "C"{
#include "common/glutil.h"
}

#include "common/touch.h"
#include "common/framebuffer.h"
#include "common/visionipc.h"
#include "common/visionimg.h"
#include "common/modeldata.h"
#include "common/params.h"

#include "cereal/gen/c/log.capnp.h"

extern "C"{
#include "slplay.h"
}

#include "messaging.hpp"

//BB include BBUIState def
#include "bbuistate.h"
//BB end


#define STATUS_STOPPED 0
#define STATUS_DISENGAGED 1
#define STATUS_ENGAGED 2
#define STATUS_WARNING 3
#define STATUS_ALERT 4

#define ALERTSIZE_NONE 0
#define ALERTSIZE_SMALL 1
#define ALERTSIZE_MID 2
#define ALERTSIZE_FULL 3

#define UI_BUF_COUNT 4
//#define SHOW_SPEEDLIMIT 1
//#define DEBUG_TURN

//#define DEBUG_FPS

const int vwp_w = 1920;
const int vwp_h = 1080;
const int nav_w = 640;
const int nav_ww= 760;
const int sbr_w = 300;
const int bdr_s = 30;
const int box_x = sbr_w+bdr_s;
const int box_y = bdr_s;
const int box_w = vwp_w-sbr_w-(bdr_s*2);
const int box_h = vwp_h-(bdr_s*2);
const int viz_w = vwp_w-(bdr_s*2);
const int header_h = 420;
const int footer_h = 280;
const int footer_y = vwp_h-bdr_s-footer_h;

const int UI_FREQ = 30;   // Hz

const int MODEL_PATH_MAX_VERTICES_CNT = 98;
const int MODEL_LANE_PATH_CNT = 3;
const int TRACK_POINTS_MAX_CNT = 50 * 2;

const uint8_t bg_colors[][4] = {
  [STATUS_STOPPED] = {0x07, 0x23, 0x39, 0xff},
  [STATUS_DISENGAGED] = {0x17, 0x33, 0x49, 0xff},
  [STATUS_ENGAGED] = {0x17, 0x86, 0x44, 0xff},
  [STATUS_WARNING] = {0xDA, 0x6F, 0x25, 0xff},
  [STATUS_ALERT] = {0xC9, 0x22, 0x31, 0xff},
};

const uint8_t alert_colors[][4] = {
  [STATUS_STOPPED] = {0x07, 0x23, 0x39, 0xf1},
  [STATUS_DISENGAGED] = {0x17, 0x33, 0x49, 0xc8},
  [STATUS_ENGAGED] = {0x17, 0x86, 0x44, 0xf1},
  [STATUS_WARNING] = {0xDA, 0x6F, 0x25, 0xf1},
  [STATUS_ALERT] = {0xC9, 0x22, 0x31, 0xf1},
};

const int alert_sizes[] = {
  [ALERTSIZE_NONE] = 0,
  [ALERTSIZE_SMALL] = 241,
  [ALERTSIZE_MID] = 390,
  [ALERTSIZE_FULL] = vwp_h,
};

const int SET_SPEED_NA = 255;

// TODO: this is also hardcoded in common/transformations/camera.py
const mat3 intrinsic_matrix = (mat3){{
  910., 0., 582.,
  0., 910., 437.,
  0.,   0.,   1.
}};

typedef enum cereal_CarControl_HUDControl_AudibleAlert AudibleAlert;

typedef struct UIScene {
  int frontview;
  int fullview;

  int transformed_width, transformed_height;

  uint64_t model_ts;
  ModelData model;

  float mpc_x[50];
  float mpc_y[50];

  bool world_objects_visible;
  mat4 extrinsic_matrix;      // Last row is 0 so we can use mat4.

  float v_cruise;
  uint64_t v_cruise_update_ts;
  float v_ego;
  bool decel_for_model;

  float speedlimit;
  bool speedlimit_valid;
  bool map_valid;

  float curvature;
  int engaged;
  bool engageable;
  bool monitoring_active;

  bool uilayout_sidebarcollapsed;
  bool uilayout_mapenabled;
  // responsive layout
  int ui_viz_rx;
  int ui_viz_rw;
  int ui_viz_ro;

  int lead_status;
  float lead_d_rel, lead_y_rel, lead_v_rel;

  int front_box_x, front_box_y, front_box_width, front_box_height;

  uint64_t alert_ts;
  char alert_text1[1024];
  char alert_text2[1024];
  uint8_t alert_size;
  float alert_blinkingrate;

  float awareness_status;

  uint64_t started_ts;


  // Used to show gps planner status
  bool gps_planner_active;
} UIScene;

typedef struct {
  float x, y;
}vertex_data;

typedef struct {
  vertex_data v[MODEL_PATH_MAX_VERTICES_CNT];
  int cnt;
} model_path_vertices_data;

typedef struct {
  vertex_data v[TRACK_POINTS_MAX_CNT];
  int cnt;
} track_vertices_data;


typedef struct UIState {
  //BB define BBUIState
  BBUIState b;
  //BB end
  pthread_mutex_t lock;
  pthread_cond_t bg_cond;

  FramebufferState *fb;
  int fb_w, fb_h;
  EGLDisplay display;
  EGLSurface surface;

  NVGcontext *vg;

  int font_courbd;
  int font_sans_regular;
  int font_sans_semibold;
  int font_sans_bold;
  int img_wheel;
  int img_turn;
  int img_face;
  int img_map;

  // Sockets
  Context *ctx;
  SubSocket *thermal_sock;
  SubSocket *model_sock;
  SubSocket *controlsstate_sock;
  SubSocket *livecalibration_sock;
  SubSocket *radarstate_sock;
  SubSocket *plus_sock;
  SubSocket *map_data_sock;
  SubSocket *uilayout_sock;
  Poller * poller;

  int plus_state;

  // vision state
  bool vision_connected;
  bool vision_connect_firstrun;
  int ipc_fd;

  VIPCBuf bufs[UI_BUF_COUNT];
  VIPCBuf front_bufs[UI_BUF_COUNT];
  int cur_vision_idx;
  int cur_vision_front_idx;

  GLuint frame_program;
  GLuint frame_texs[UI_BUF_COUNT];
  EGLImageKHR khr[UI_BUF_COUNT];
  void *priv_hnds[UI_BUF_COUNT];
  GLuint frame_front_texs[UI_BUF_COUNT];
  EGLImageKHR khr_front[UI_BUF_COUNT];
  void *priv_hnds_front[UI_BUF_COUNT];

  GLint frame_pos_loc, frame_texcoord_loc;
  GLint frame_texture_loc, frame_transform_loc;

  GLuint line_program;
  GLint line_pos_loc, line_color_loc;
  GLint line_transform_loc;

  int rgb_width, rgb_height, rgb_stride;
  size_t rgb_buf_len;
  mat4 rgb_transform;

  int rgb_front_width, rgb_front_height, rgb_front_stride;
  size_t rgb_front_buf_len;

  UIScene scene;

  bool awake;
  int awake_timeout;

  int volume_timeout;
  int controls_timeout;
  int alert_sound_timeout;
  int speed_lim_off_timeout;
  int is_metric_timeout;
  int longitudinal_control_timeout;
  int limit_set_speed_timeout;

  bool controls_seen;

  int status;
  bool is_metric;
  bool longitudinal_control;
  bool limit_set_speed;
  float speed_lim_off;
  bool is_ego_over_limit;
  char alert_type[64];
  AudibleAlert alert_sound;
  int alert_size;
  float alert_blinking_alpha;
  bool alert_blinked;

  float light_sensor;

  int touch_fd;

  // Hints for re-calculations and redrawing
  bool model_changed;
  bool livempc_or_radarstate_changed;

  GLuint frame_vao[2], frame_vbo[2], frame_ibo[2];
  mat4 rear_frame_mat, front_frame_mat;

  model_path_vertices_data model_path_vertices[MODEL_LANE_PATH_CNT * 2];

  track_vertices_data track_vertices[2];
} UIState;



static int last_brightness = -1;
static void set_brightness(UIState *s, int brightness) {
  if (last_brightness != brightness && (s->awake || brightness == 0)) {
    FILE *f = fopen("/sys/class/leds/lcd-backlight/brightness", "wb");
    if (f != NULL) {
      fprintf(f, "%d", brightness);
      fclose(f);
      last_brightness = brightness;
    }
  }
}

static void set_awake(UIState *s, bool awake) {
  if (awake) {
    // 30 second timeout at 30 fps
    if (((s->b.tri_state_switch == 3) || (s->b.keepEonOff)) && !s->b.recording) {
      s->awake_timeout = 3*30;
    } else {
      s->awake_timeout = 30*30;
    }
  }
  if (s->awake != awake) {
    s->awake = awake;

    // TODO: replace command_awake and command_sleep with direct calls to android
    if (awake) {
      LOG("awake normal");
      set_brightness(s, 150);
      system("service call window 18 i32 1");  // enable event processing
      framebuffer_set_power(s->fb, HWC_POWER_MODE_NORMAL);
    } else {
      LOGW("awake off");
      set_brightness(s, 0);
      system("service call window 18 i32 0");  // disable event processing
      framebuffer_set_power(s->fb, HWC_POWER_MODE_OFF);
    }
  }
}

#include "dashcam.h"
#include "bbui.h"
#include "tbui.h"

static void set_volume(UIState *s, int volume) {
  char volume_change_cmd[64];
  sprintf(volume_change_cmd, "service call audio 3 i32 3 i32 %d i32 1", volume);

  // 5 second timeout at 60fps
  s->volume_timeout = 5 * UI_FREQ;
  int volume_changed = system(volume_change_cmd);
}

volatile sig_atomic_t do_exit = 0;
static void set_do_exit(int sig) {
  do_exit = 1;
}

static void read_param_bool(bool* param, const char* param_name) {
  char *s;
  const int result = read_db_value(NULL, param_name, &s, NULL);
  if (result == 0) {
    *param = s[0] == '1';
    free(s);
  }
}

static void read_param_float(float* param, const char* param_name) {
  char *s;
  const int result = read_db_value(NULL, param_name, &s, NULL);
  if (result == 0) {
    *param = strtod(s, NULL);
    free(s);
  }
}

static void read_param_bool_timeout(bool* param, const char* param_name, int* timeout) {
  if (*timeout > 0){
    (*timeout)--;
  } else {
    read_param_bool(param, param_name);
    *timeout = 2 * UI_FREQ; // 0.5Hz
  }
}

static void read_param_float_timeout(float* param, const char* param_name, int* timeout) {
  if (*timeout > 0){
    (*timeout)--;
  } else {
    read_param_float(param, param_name);
    *timeout = 2 * UI_FREQ; // 0.5Hz
  }
}

static const char frame_vertex_shader[] =
  "attribute vec4 aPosition;\n"
  "attribute vec4 aTexCoord;\n"
  "uniform mat4 uTransform;\n"
  "varying vec4 vTexCoord;\n"
  "void main() {\n"
  "  gl_Position = uTransform * aPosition;\n"
  "  vTexCoord = aTexCoord;\n"
  "}\n";

static const char frame_fragment_shader[] =
  "precision mediump float;\n"
  "uniform sampler2D uTexture;\n"
  "varying vec4 vTexCoord;\n"
  "void main() {\n"
  "  gl_FragColor = texture2D(uTexture, vTexCoord.xy);\n"
  "}\n";

static const char line_vertex_shader[] =
  "attribute vec4 aPosition;\n"
  "attribute vec4 aColor;\n"
  "uniform mat4 uTransform;\n"
  "varying vec4 vColor;\n"
  "void main() {\n"
  "  gl_Position = uTransform * aPosition;\n"
  "  vColor = aColor;\n"
  "}\n";

static const char line_fragment_shader[] =
  "precision mediump float;\n"
  "uniform sampler2D uTexture;\n"
  "varying vec4 vColor;\n"
  "void main() {\n"
  "  gl_FragColor = vColor;\n"
  "}\n";


static const mat4 device_transform = {{
  1.0,  0.0, 0.0, 0.0,
  0.0,  1.0, 0.0, 0.0,
  0.0,  0.0, 1.0, 0.0,
  0.0,  0.0, 0.0, 1.0,
}};

// frame from 4/3 to box size with a 2x zoom
static const mat4 frame_transform = {{
  2*(4./3.)/((float)viz_w/box_h), 0.0, 0.0, 0.0,
                                           0.0, 2.0, 0.0, 0.0,
                                           0.0, 0.0, 1.0, 0.0,
                                           0.0, 0.0, 0.0, 1.0,
}};

// frame from 4/3 to 16/9 display
static const mat4 full_to_wide_frame_transform = {{
  .75,  0.0, 0.0, 0.0,
  0.0,  1.0, 0.0, 0.0,
  0.0,  0.0, 1.0, 0.0,
  0.0,  0.0, 0.0, 1.0,
}};

typedef struct {
  AudibleAlert alert;
  const char* uri;
  bool loop;
} sound_file;

sound_file sound_table[] = {
  { cereal_CarControl_HUDControl_AudibleAlert_chimeDisengage, "../assets/sounds/disengaged.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_chimeEngage, "../assets/sounds/engaged.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_chimeWarning1, "../assets/sounds/warning_1.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_chimeWarning2, "../assets/sounds/warning_2.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_chimeWarningRepeat, "../assets/sounds/warning_repeat.wav", true },
  { cereal_CarControl_HUDControl_AudibleAlert_chimeError, "../assets/sounds/error.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_chimePrompt, "../assets/sounds/error.wav", false },
  { cereal_CarControl_HUDControl_AudibleAlert_none, NULL, false },
};

sound_file* get_sound_file(AudibleAlert alert) {
  for (sound_file *s = sound_table; s->alert != cereal_CarControl_HUDControl_AudibleAlert_none; s++) {
    if (s->alert == alert) {
      return s;
    }
  }

  return NULL;
}

void play_alert_sound(AudibleAlert alert) {
  sound_file* sound = get_sound_file(alert);
  char* error = NULL;

  slplay_play(sound->uri, sound->loop, &error);
  if(error) {
    LOGW("error playing sound: %s", error);
  }
}

void stop_alert_sound(AudibleAlert alert) {
  sound_file* sound = get_sound_file(alert);
  char* error = NULL;

  slplay_stop_uri(sound->uri, &error);
  if(error) {
    LOGW("error stopping sound: %s", error);
  }
}

void ui_sound_init(char **error) {
  slplay_setup(error);
  if (*error) return;

  for (sound_file *s = sound_table; s->alert != cereal_CarControl_HUDControl_AudibleAlert_none; s++) {
    slplay_create_player_for_uri(s->uri, error);
    if (*error) return;
  }
}

static void ui_init(UIState *s) {
  memset(s, 0, sizeof(UIState));

  pthread_mutex_init(&s->lock, NULL);
  pthread_cond_init(&s->bg_cond, NULL);

  s->ctx = Context::create();
  s->thermal_sock = SubSocket::create(s->ctx, "thermal");
  s->model_sock = SubSocket::create(s->ctx, "model");
  s->controlsstate_sock = SubSocket::create(s->ctx, "controlsState");
  s->uilayout_sock = SubSocket::create(s->ctx, "uiLayoutState");
  s->livecalibration_sock = SubSocket::create(s->ctx, "liveCalibration");
  s->radarstate_sock = SubSocket::create(s->ctx, "radarState");
  s->plus_sock = SubSocket::create(s->ctx, "plusFrame");
  s->map_data_sock = SubSocket::create(s->ctx, "liveMapData");
  s->poller = Poller::create({
                              s->thermal_sock,
                              s->model_sock,
                              s->controlsstate_sock,
                              s->uilayout_sock,
                              s->livecalibration_sock,
                              s->radarstate_sock,
                              s->plus_sock
                             });

#ifdef SHOW_SPEEDLIMIT
  s->map_data_sock = SubSock::create(s->ctx, "liveMapData");
  s->poller.registerSocket(s->map_data_sock);
#endif

  s->ipc_fd = -1;

  // init display
  s->fb = framebuffer_init("ui", 0x00010000, true,
                           &s->display, &s->surface, &s->fb_w, &s->fb_h);
  assert(s->fb);

  set_awake(s, true);

  // init drawing
  s->vg = nvgCreateGLES3(NVG_ANTIALIAS | NVG_STENCIL_STROKES | NVG_DEBUG);
  assert(s->vg);

  s->font_courbd = nvgCreateFont(s->vg, "courbd", "../assets/fonts/courbd.ttf");
  assert(s->font_courbd >= 0);
  s->font_sans_regular = nvgCreateFont(s->vg, "sans-regular", "../assets/fonts/opensans_regular.ttf");
  assert(s->font_sans_regular >= 0);
  s->font_sans_semibold = nvgCreateFont(s->vg, "sans-semibold", "../assets/fonts/opensans_semibold.ttf");
  assert(s->font_sans_semibold >= 0);
  s->font_sans_bold = nvgCreateFont(s->vg, "sans-bold", "../assets/fonts/opensans_bold.ttf");
  assert(s->font_sans_bold >= 0);

  assert(s->img_wheel >= 0);
  s->img_wheel = nvgCreateImage(s->vg, "../assets/img_chffr_wheel.png", 1);

  assert(s->img_turn >= 0);
  s->img_turn = nvgCreateImage(s->vg, "../assets/img_trafficSign_turn.png", 1);

  assert(s->img_face >= 0);
  s->img_face = nvgCreateImage(s->vg, "../assets/img_driver_face.png", 1);

  assert(s->img_map >= 0);
  s->img_map = nvgCreateImage(s->vg, "../assets/img_map.png", 1);

  // init gl
  s->frame_program = load_program(frame_vertex_shader, frame_fragment_shader);
  assert(s->frame_program);

  s->frame_pos_loc = glGetAttribLocation(s->frame_program, "aPosition");
  s->frame_texcoord_loc = glGetAttribLocation(s->frame_program, "aTexCoord");

  s->frame_texture_loc = glGetUniformLocation(s->frame_program, "uTexture");
  s->frame_transform_loc = glGetUniformLocation(s->frame_program, "uTransform");

  s->line_program = load_program(line_vertex_shader, line_fragment_shader);
  assert(s->line_program);

  s->line_pos_loc = glGetAttribLocation(s->line_program, "aPosition");
  s->line_color_loc = glGetAttribLocation(s->line_program, "aColor");
  s->line_transform_loc = glGetUniformLocation(s->line_program, "uTransform");

  glViewport(0, 0, s->fb_w, s->fb_h);

  glDisable(GL_DEPTH_TEST);

  assert(glGetError() == GL_NO_ERROR);

  for(int i = 0; i < 2; i++) {
    float x1, x2, y1, y2;
    if (i == 1) {
      // flip horizontally so it looks like a mirror
      x1 = 0.0;
      x2 = 1.0;
      y1 = 1.0;
      y2 = 0.0;
    } else {
      x1 = 1.0;
      x2 = 0.0;
      y1 = 1.0;
      y2 = 0.0;
    }
    const uint8_t frame_indicies[] = {0, 1, 2, 0, 2, 3};
    const float frame_coords[4][4] = {
      {-1.0, -1.0, x2, y1}, //bl
      {-1.0,  1.0, x2, y2}, //tl
      { 1.0,  1.0, x1, y2}, //tr
      { 1.0, -1.0, x1, y1}, //br
    };

    glGenVertexArrays(1,&s->frame_vao[i]);
    glBindVertexArray(s->frame_vao[i]);
    glGenBuffers(1, &s->frame_vbo[i]);
    glBindBuffer(GL_ARRAY_BUFFER, s->frame_vbo[i]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(frame_coords), frame_coords, GL_STATIC_DRAW);
    glEnableVertexAttribArray(s->frame_pos_loc);
    glVertexAttribPointer(s->frame_pos_loc, 2, GL_FLOAT, GL_FALSE,
                          sizeof(frame_coords[0]), (const void *)0);
    glEnableVertexAttribArray(s->frame_texcoord_loc);
    glVertexAttribPointer(s->frame_texcoord_loc, 2, GL_FLOAT, GL_FALSE,
                          sizeof(frame_coords[0]), (const void *)(sizeof(float) * 2));
    glGenBuffers(1, &s->frame_ibo[i]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s->frame_ibo[i]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(frame_indicies), frame_indicies, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindVertexArray(0);
  }

  s->model_changed = false;
  s->livempc_or_radarstate_changed = false;

  s->front_frame_mat = matmul(device_transform, full_to_wide_frame_transform);
  s->rear_frame_mat = matmul(device_transform, frame_transform);

  for(int i = 0;i < UI_BUF_COUNT; i++) {
    s->khr[i] = NULL;
    s->priv_hnds[i] = NULL;
    s->khr_front[i] = NULL;
    s->priv_hnds_front[i] = NULL;
  }
}

static void ui_init_vision(UIState *s, const VisionStreamBufs back_bufs,
                           int num_back_fds, const int *back_fds,
                           const VisionStreamBufs front_bufs, int num_front_fds,
                           const int *front_fds) {
  const VisionUIInfo ui_info = back_bufs.buf_info.ui_info;

  assert(num_back_fds == UI_BUF_COUNT);
  assert(num_front_fds == UI_BUF_COUNT);

  vipc_bufs_load(s->bufs, &back_bufs, num_back_fds, back_fds);
  vipc_bufs_load(s->front_bufs, &front_bufs, num_front_fds, front_fds);

  s->cur_vision_idx = -1;
  s->cur_vision_front_idx = -1;

  s->scene = (UIScene){
      .frontview = getenv("FRONTVIEW") != NULL,
      .fullview = getenv("FULLVIEW") != NULL,
      .transformed_width = ui_info.transformed_width,
      .transformed_height = ui_info.transformed_height,
      .front_box_x = ui_info.front_box_x,
      .front_box_y = ui_info.front_box_y,
      .front_box_width = ui_info.front_box_width,
      .front_box_height = ui_info.front_box_height,
      .world_objects_visible = false,  // Invisible until we receive a calibration message.
      .gps_planner_active = false,
  };

  s->rgb_width = back_bufs.width;
  s->rgb_height = back_bufs.height;
  s->rgb_stride = back_bufs.stride;
  s->rgb_buf_len = back_bufs.buf_len;

  s->rgb_front_width = front_bufs.width;
  s->rgb_front_height = front_bufs.height;
  s->rgb_front_stride = front_bufs.stride;
  s->rgb_front_buf_len = front_bufs.buf_len;

  s->rgb_transform = (mat4){{
    2.0f/s->rgb_width, 0.0f, 0.0f, -1.0f,
    0.0f, 2.0f/s->rgb_height, 0.0f, -1.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f,
  }};

  read_param_float(&s->speed_lim_off, "SpeedLimitOffset");
  read_param_bool(&s->is_metric, "IsMetric");
  read_param_bool(&s->longitudinal_control, "LongitudinalControl");
  read_param_bool(&s->limit_set_speed, "LimitSetSpeed");

  // Set offsets so params don't get read at the same time
  s->longitudinal_control_timeout = UI_FREQ / 3;
  s->is_metric_timeout = UI_FREQ / 2;
  s->limit_set_speed_timeout = UI_FREQ;
}

// Projects a point in car to space to the corresponding point in full frame
// image space.
vec3 car_space_to_full_frame(const UIState *s, vec4 car_space_projective) {
  const UIScene *scene = &s->scene;

  // We'll call the car space point p.
  // First project into normalized image coordinates with the extrinsics matrix.
  const vec4 Ep4 = matvecmul(scene->extrinsic_matrix, car_space_projective);

  // The last entry is zero because of how we store E (to use matvecmul).
  const vec3 Ep = {{Ep4.v[0], Ep4.v[1], Ep4.v[2]}};
  const vec3 KEp = matvecmul3(intrinsic_matrix, Ep);

  // Project.
  const vec3 p_image = {{KEp.v[0] / KEp.v[2], KEp.v[1] / KEp.v[2], 1.}};
  return p_image;
}

// Calculate an interpolation between two numbers at a specific increment
static float lerp(float v0, float v1, float t) {
  return (1 - t) * v0 + t * v1;
}

static void draw_chevron(UIState *s, float x_in, float y_in, float sz,
                          NVGcolor fillColor, NVGcolor glowColor) {
  const UIScene *scene = &s->scene;

  nvgSave(s->vg);

  nvgTranslate(s->vg, 240.0f, 0.0);
  nvgTranslate(s->vg, -1440.0f / 2, -1080.0f / 2);
  nvgScale(s->vg, 2.0, 2.0);
  nvgScale(s->vg, 1440.0f / s->rgb_width, 1080.0f / s->rgb_height);

  const vec4 p_car_space = (vec4){{x_in, y_in, 0., 1.}};
  const vec3 p_full_frame = car_space_to_full_frame(s, p_car_space);

  sz *= 30;
  sz /= (x_in / 3 + 30);
  if (sz > 30) sz = 30;
  if (sz < 15) sz = 15;

  float x = p_full_frame.v[0];
  float y = p_full_frame.v[1];

  // glow
  nvgBeginPath(s->vg);
  float g_xo = sz/5;
  float g_yo = sz/10;
  //BB added for printing the car
  //if position is 3 do nothing
  if (s->b.tri_state_switch == 3) {
    return;
  }
  if ((s->b.tri_state_switch == 2) && (s->b.icShowCar)) {
    nvgRestore(s->vg);
    bb_ui_draw_car(s);
  } else {
    if (x >= 0 && y >= 0.) {
        nvgMoveTo(s->vg, x+(sz*1.35)+g_xo, y+sz+g_yo);
        nvgLineTo(s->vg, x, y-g_xo);
        nvgLineTo(s->vg, x-(sz*1.35)-g_xo, y+sz+g_yo);
        nvgLineTo(s->vg, x+(sz*1.35)+g_xo, y+sz+g_yo);
        nvgClosePath(s->vg);
    }
    nvgFillColor(s->vg, glowColor);
    nvgFill(s->vg);

    // chevron
    nvgBeginPath(s->vg);
    if (x >= 0 && y >= 0.) {
        nvgMoveTo(s->vg, x+(sz*1.25), y+sz);
        nvgLineTo(s->vg, x, y);
        nvgLineTo(s->vg, x-(sz*1.25), y+sz);
        nvgLineTo(s->vg, x+(sz*1.25), y+sz);
        nvgClosePath(s->vg);
    }
    nvgFillColor(s->vg, fillColor);
    nvgFill(s->vg);
  }
  nvgRestore(s->vg);
}

static void ui_draw_lane_line(UIState *s, const model_path_vertices_data *pvd, NVGcolor color) {
  const UIScene *scene = &s->scene;

  nvgSave(s->vg);
  nvgTranslate(s->vg, 240.0f, 0.0); // rgb-box space
  nvgTranslate(s->vg, -1440.0f / 2, -1080.0f / 2); // zoom 2x
  nvgScale(s->vg, 2.0, 2.0);
  nvgScale(s->vg, 1440.0f / s->rgb_width, 1080.0f / s->rgb_height);
  nvgBeginPath(s->vg);

  bool started = false;
  for (int i=0; i<pvd->cnt; i++) {
    if (pvd->v[i].x < 0 || pvd->v[i].y < 0.) {
      continue;
    }
    if (!started) {
      nvgMoveTo(s->vg, pvd->v[i].x, pvd->v[i].y);
      started = true;
    } else {
      nvgLineTo(s->vg, pvd->v[i].x, pvd->v[i].y);
    }
  }

  nvgClosePath(s->vg);
  nvgFillColor(s->vg, color);
  nvgFill(s->vg);
  nvgRestore(s->vg);
}

static void update_track_data(UIState *s, bool is_mpc, track_vertices_data *pvd) {
  const UIScene *scene = &s->scene;
  const PathData path = scene->model.path;
  const float *mpc_x_coords = &scene->mpc_x[0];
  const float *mpc_y_coords = &scene->mpc_y[0];

  bool started = false;
  float off = is_mpc?0.3:0.5;
  float lead_d = scene->lead_d_rel*2.;
  float path_height = is_mpc?(lead_d>5.)?fmin(lead_d, 25.)-fmin(lead_d*0.35, 10.):20.
                            :(lead_d>0.)?fmin(lead_d, 50.)-fmin(lead_d*0.35, 10.):49.;
  pvd->cnt = 0;
  // left side up
  for (int i=0; i<=path_height; i++) {
    float px, py, mpx;
    if (is_mpc) {
      mpx = i==0?0.0:mpc_x_coords[i];
      px = lerp(mpx+1.0, mpx, i/100.0);
      py = mpc_y_coords[i] - off;
    } else {
      px = lerp(i+1.0, i, i/100.0);
      py = path.points[i] - off;
    }

    vec4 p_car_space = (vec4){{px, py, 0., 1.}};
    vec3 p_full_frame = car_space_to_full_frame(s, p_car_space);
    if (p_full_frame.v[0] < 0. || p_full_frame.v[1] < 0.) {
      continue;
    }
    pvd->v[pvd->cnt].x = p_full_frame.v[0];
    pvd->v[pvd->cnt].y = p_full_frame.v[1];
    pvd->cnt += 1;
  }

  // right side down
  for (int i=path_height; i>=0; i--) {
    float px, py, mpx;
    if (is_mpc) {
      mpx = i==0?0.0:mpc_x_coords[i];
      px = lerp(mpx+1.0, mpx, i/100.0);
      py = mpc_y_coords[i] + off;
    } else {
      px = lerp(i+1.0, i, i/100.0);
      py = path.points[i] + off;
    }

    vec4 p_car_space = (vec4){{px, py, 0., 1.}};
    vec3 p_full_frame = car_space_to_full_frame(s, p_car_space);
    pvd->v[pvd->cnt].x = p_full_frame.v[0];
    pvd->v[pvd->cnt].y = p_full_frame.v[1];
    pvd->cnt += 1;
  }
}

static void update_all_track_data(UIState *s) {
  const UIScene *scene = &s->scene;
  // Draw vision path
  update_track_data(s, false, &s->track_vertices[0]);

  if (scene->engaged) {
    // Draw MPC path when engaged
    update_track_data(s, true, &s->track_vertices[1]);
  }
}


static void ui_draw_track(UIState *s, bool is_mpc, track_vertices_data *pvd) {
const UIScene *scene = &s->scene;
  const PathData path = scene->model.path;
  const float *mpc_x_coords = &scene->mpc_x[0];
  const float *mpc_y_coords = &scene->mpc_y[0];

  nvgSave(s->vg);
  nvgTranslate(s->vg, 240.0f, 0.0); // rgb-box space
  nvgTranslate(s->vg, -1440.0f / 2, -1080.0f / 2); // zoom 2x
  nvgScale(s->vg, 2.0, 2.0);
  nvgScale(s->vg, 1440.0f / s->rgb_width, 1080.0f / s->rgb_height);
  nvgBeginPath(s->vg);

  bool started = false;
  float off = is_mpc?0.3:0.5;
  float lead_d = scene->lead_d_rel*2.;
  float path_height = is_mpc?(lead_d>5.)?fmin(lead_d, 25.)-fmin(lead_d*0.35, 10.):20.
                            :(lead_d>0.)?fmin(lead_d, 50.)-fmin(lead_d*0.35, 10.):49.;
  int vi = 0;
  for(int i = 0;i < pvd->cnt;i++) {
    if (pvd->v[i].x < 0 || pvd->v[i].y < 0) {
      continue;
    }

    if (!started) {
      nvgMoveTo(s->vg, pvd->v[i].x, pvd->v[i].y);
      started = true;
    } else {
      nvgLineTo(s->vg, pvd->v[i].x, pvd->v[i].y);
    }
  }

  nvgClosePath(s->vg);

  NVGpaint track_bg;
  if (is_mpc) {
    // Draw colored MPC track
    const uint8_t *clr = bg_colors[s->status];
    track_bg = nvgLinearGradient(s->vg, vwp_w, vwp_h, vwp_w, vwp_h*.4,
      nvgRGBA(clr[0], clr[1], clr[2], 255), nvgRGBA(clr[0], clr[1], clr[2], 255/2));
  } else {
    // Draw white vision track
    track_bg = nvgLinearGradient(s->vg, vwp_w, vwp_h, vwp_w, vwp_h*.4,
      nvgRGBA(255, 255, 255, 255), nvgRGBA(255, 255, 255, 0));
  }

  nvgFillPaint(s->vg, track_bg);
  nvgFill(s->vg);
  nvgRestore(s->vg);
}

static void draw_steering(UIState *s, float curvature) {
  float points[50];
  for (int i = 0; i < 50; i++) {
    float y_actual = i * tan(asin(clamp(i * curvature, -0.999, 0.999)) / 2.);
    points[i] = y_actual;
  }

  // ui_draw_lane_edge(s, points, 0.0, nvgRGBA(0, 0, 255, 128), 5);
}

static void draw_frame(UIState *s) {
  const UIScene *scene = &s->scene;

  float x1, x2, y1, y2;
  if (s->scene.frontview) {
    glBindVertexArray(s->frame_vao[1]);
  } else {
    glBindVertexArray(s->frame_vao[0]);
  }

  mat4 *out_mat;
  if (s->scene.frontview || s->scene.fullview) {
    out_mat = &s->front_frame_mat;
  } else {
    out_mat = &s->rear_frame_mat;
  }
  glActiveTexture(GL_TEXTURE0);
  //BB added to suppress video printing
  if ((s->b.tri_state_switch == 1) || (s->b.tri_state_switch == 4)) {
    if (s->scene.frontview && s->cur_vision_front_idx >= 0) {
      glBindTexture(GL_TEXTURE_2D, s->frame_front_texs[s->cur_vision_front_idx]);
    } else if (!scene->frontview && s->cur_vision_idx >= 0) {
      glBindTexture(GL_TEXTURE_2D, s->frame_texs[s->cur_vision_idx]);
    }
  }
  //BB end  

  glUseProgram(s->frame_program);
  glUniform1i(s->frame_texture_loc, 0);
  glUniformMatrix4fv(s->frame_transform_loc, 1, GL_TRUE, out_mat->v);

  assert(glGetError() == GL_NO_ERROR);
  glEnableVertexAttribArray(0);
  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, (const void*)0);
  glDisableVertexAttribArray(0);
  glBindVertexArray(0);
}

static inline bool valid_frame_pt(UIState *s, float x, float y) {
  return x >= 0 && x <= s->rgb_width && y >= 0 && y <= s->rgb_height;

}
static void update_lane_line_data(UIState *s, const float *points, float off, bool is_ghost, model_path_vertices_data *pvd) {
  pvd->cnt = 0;
  for (int i = 0; i < MODEL_PATH_MAX_VERTICES_CNT / 2; i++) {
    float px = (float)i;
    float py = points[i] - off;
    const vec4 p_car_space = (vec4){{px, py, 0., 1.}};
    const vec3 p_full_frame = car_space_to_full_frame(s, p_car_space);
    if(!valid_frame_pt(s, p_full_frame.v[0], p_full_frame.v[1]))
      continue;
    pvd->v[pvd->cnt].x = p_full_frame.v[0];
    pvd->v[pvd->cnt].y = p_full_frame.v[1];
    pvd->cnt += 1;
  }
  for (int i = MODEL_PATH_MAX_VERTICES_CNT / 2; i > 0; i--) {
    float px = (float)i;
    float py = is_ghost?(points[i]-off):(points[i]+off);
    const vec4 p_car_space = (vec4){{px, py, 0., 1.}};
    const vec3 p_full_frame = car_space_to_full_frame(s, p_car_space);
    if(!valid_frame_pt(s, p_full_frame.v[0], p_full_frame.v[1]))
      continue;
    pvd->v[pvd->cnt].x = p_full_frame.v[0];
    pvd->v[pvd->cnt].y = p_full_frame.v[1];
    pvd->cnt += 1;
  }
}

static void update_all_lane_lines_data(UIState *s, const PathData path, model_path_vertices_data *pstart) {
  update_lane_line_data(s, path.points, 0.025*path.prob, false, pstart);
  float var = fmin(path.std, 0.7);
  update_lane_line_data(s, path.points, -var, true, pstart + 1);
  update_lane_line_data(s, path.points, var, true, pstart + 2);
}

static void ui_draw_lane(UIState *s, const PathData *path, model_path_vertices_data *pstart, NVGcolor color) {
  ui_draw_lane_line(s, pstart, color);
  float var = fmin(path->std, 0.7);
  color.a /= 4;
  ui_draw_lane_line(s, pstart + 1, color);
  ui_draw_lane_line(s, pstart + 2, color);
}

static void ui_draw_vision_lanes(UIState *s) {
  const UIScene *scene = &s->scene;
  //draw nothing if position is 3
  if (s->b.tri_state_switch == 3) {
    return;
  }
  //BB add to draw our lanes
  if (s->b.tri_state_switch == 2) {
    bb_draw_lane_fill(s);
  }
  model_path_vertices_data *pvd = &s->model_path_vertices[0];
  if(s->model_changed) {
    update_all_lane_lines_data(s, scene->model.left_lane, pvd);
    update_all_lane_lines_data(s, scene->model.right_lane, pvd + MODEL_LANE_PATH_CNT);
    s->model_changed = false;
  }
  // Draw left lane edge
  ui_draw_lane(
      s, &scene->model.left_lane,
      pvd,
      nvgRGBAf(1.0, 1.0, 1.0, scene->model.left_lane.prob));

  // Draw right lane edge
  ui_draw_lane(
      s, &scene->model.right_lane,
      pvd + MODEL_LANE_PATH_CNT,
      nvgRGBAf(1.0, 1.0, 1.0, scene->model.right_lane.prob));

  if(s->livempc_or_radarstate_changed) {
    update_all_track_data(s);
    s->livempc_or_radarstate_changed = false;
  }
  // Draw vision path
  ui_draw_track(s, false, &s->track_vertices[0]);
  if (scene->engaged) {
    // Draw MPC path when engaged
    ui_draw_track(s, true, &s->track_vertices[1]);
  }
}

// Draw all world space objects.
static void ui_draw_world(UIState *s) {
  const UIScene *scene = &s->scene;
  if (!scene->world_objects_visible) {
    return;
  }

  if ((nanos_since_boot() - scene->model_ts) < 1000000000ULL) {
    // Draw lane edges and vision/mpc tracks
    ui_draw_vision_lanes(s);
  }

  if (scene->lead_status) {
    // Draw lead car indicator
    float fillAlpha = 0;
    float speedBuff = 10.;
    float leadBuff = 40.;
    if (scene->lead_d_rel < leadBuff) {
      fillAlpha = 255*(1.0-(scene->lead_d_rel/leadBuff));
      if (scene->lead_v_rel < 0) {
        fillAlpha += 255*(-1*(scene->lead_v_rel/speedBuff));
      }
      fillAlpha = (int)(fmin(fillAlpha, 255));
    }
    draw_chevron(s, scene->lead_d_rel+2.7, scene->lead_y_rel, 25,
                  nvgRGBA(201, 34, 49, fillAlpha), nvgRGBA(218, 202, 37, 255));
  }
}

static void ui_draw_vision_maxspeed(UIState *s) {
  /*if (!s->longitudinal_control){
    return;
  }*/

  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;

  char maxspeed_str[32];
  float maxspeed = s->scene.v_cruise;
  int maxspeed_calc = maxspeed * 0.6225 + 0.5;
  float speedlimit = s->scene.speedlimit;
  int speedlim_calc = speedlimit * 2.2369363 + 0.5;
  int speed_lim_off = s->speed_lim_off * 2.2369363 + 0.5;
  if (s->is_metric) {
    maxspeed_calc = maxspeed + 0.5;
    speedlim_calc = speedlimit * 3.6 + 0.5;
    speed_lim_off = s->speed_lim_off * 3.6 + 0.5;
  }

  bool is_cruise_set = (maxspeed != 0 && maxspeed != SET_SPEED_NA);
  bool is_speedlim_valid = s->scene.speedlimit_valid;
  bool is_set_over_limit = is_speedlim_valid && s->scene.engaged &&
                       is_cruise_set && maxspeed_calc > (speedlim_calc + speed_lim_off);

  int viz_maxspeed_w = 184;
  int viz_maxspeed_h = 202;
  int viz_maxspeed_x = (ui_viz_rx + (bdr_s*2));
  int viz_maxspeed_y = (box_y + (bdr_s*1.5));
  int viz_maxspeed_xo = 180;
  if (s->b.tri_state_switch != 2) {
    viz_maxspeed_xo = 0;
  }
  
  if (s->b.tri_state_switch == 2) {
    viz_maxspeed_w += viz_maxspeed_xo;
    viz_maxspeed_x += viz_maxspeed_w - (viz_maxspeed_xo * 2);
  }

  // Draw Background
  nvgBeginPath(s->vg);
  nvgRoundedRect(s->vg, viz_maxspeed_x, viz_maxspeed_y, viz_maxspeed_w, viz_maxspeed_h, 30);
  if (is_set_over_limit) {
    nvgFillColor(s->vg, nvgRGBA(218, 111, 37, 180));
  } else {
    nvgFillColor(s->vg, nvgRGBA(0, 0, 0, 100));
  }
  if (s->b.tri_state_switch != 2) {
    nvgFill(s->vg);
  }

  // Draw Border
  nvgBeginPath(s->vg);
  nvgRoundedRect(s->vg, viz_maxspeed_x, viz_maxspeed_y, viz_maxspeed_w, viz_maxspeed_h, 20);
  if (s->b.tri_state_switch != 2) {
  if (is_set_over_limit) {
    nvgStrokeColor(s->vg, nvgRGBA(218, 111, 37, 255));
  } else if (is_speedlim_valid && !s->is_ego_over_limit) {
    nvgStrokeColor(s->vg, nvgRGBA(255, 255, 255, 255));
  } else if (is_speedlim_valid && s->is_ego_over_limit) {
    nvgStrokeColor(s->vg, nvgRGBA(255, 255, 255, 20));
  } else {
    nvgStrokeColor(s->vg, nvgRGBA(255, 255, 255, 100));
  }
  nvgStrokeWidth(s->vg, 10);
  nvgStroke(s->vg);
  }
  // Draw "MAX" Text
  nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_BASELINE);
  nvgFontFace(s->vg, "sans-regular");
  nvgFontSize(s->vg, 26*2.5);
  if (is_cruise_set) {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 200));
  } else {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 100));
  }
  nvgText(s->vg, viz_maxspeed_x+(viz_maxspeed_xo/2)+(viz_maxspeed_w/2), 148, "MAX", NULL);

  // Draw Speed Text
  nvgFontFace(s->vg, "sans-bold");
  nvgFontSize(s->vg, 48*2.5);
  if (is_cruise_set) {
    snprintf(maxspeed_str, sizeof(maxspeed_str), "%d", maxspeed_calc);
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
    nvgText(s->vg, viz_maxspeed_x+(viz_maxspeed_xo/2)+(viz_maxspeed_w/2), 242, maxspeed_str, NULL);
  } else {
    nvgFontFace(s->vg, "sans-semibold");
    nvgFontSize(s->vg, 42*2.5);
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 100));
    nvgText(s->vg, viz_maxspeed_x+(viz_maxspeed_xo/2)+(viz_maxspeed_w/2), 242, "N/A", NULL);
  }
  //BB START: add new measures panel  const int bb_dml_w = 180;
	bb_ui_draw_UI(s) ;
  //BB END: add new measures panel

}

static void ui_draw_vision_speedlimit(UIState *s) {
  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;

  char speedlim_str[32];
  float speedlimit = s->scene.speedlimit;
  int speedlim_calc = speedlimit * 2.2369363 + 0.5;
  if (s->is_metric) {
    speedlim_calc = speedlimit * 3.6 + 0.5;
  }

  bool is_speedlim_valid = s->scene.speedlimit_valid;
  float hysteresis_offset = 0.5;
  if (s->is_ego_over_limit) {
    hysteresis_offset = 0.0;
  }
  s->is_ego_over_limit = is_speedlim_valid && s->scene.v_ego > (speedlimit + s->speed_lim_off + hysteresis_offset);

  int viz_speedlim_w = 180;
  int viz_speedlim_h = 202;
  int viz_speedlim_x = (ui_viz_rx + (bdr_s*2));
  int viz_speedlim_y = (box_y + (bdr_s*1.5));
  if (!is_speedlim_valid) {
    viz_speedlim_w -= 5;
    viz_speedlim_h -= 10;
    viz_speedlim_x += 9;
    viz_speedlim_y += 5;
  }
  int viz_speedlim_bdr = is_speedlim_valid ? 30 : 15;

  // Draw Background
  nvgBeginPath(s->vg);
  nvgRoundedRect(s->vg, viz_speedlim_x, viz_speedlim_y, viz_speedlim_w, viz_speedlim_h, viz_speedlim_bdr);
  if (is_speedlim_valid && s->is_ego_over_limit) {
    nvgFillColor(s->vg, nvgRGBA(218, 111, 37, 180));
  } else if (is_speedlim_valid) {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
  } else {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 100));
  }
  nvgFill(s->vg);

  // Draw Border
  if (is_speedlim_valid) {
    nvgStrokeWidth(s->vg, 10);
    nvgStroke(s->vg);
    nvgBeginPath(s->vg);
    nvgRoundedRect(s->vg, viz_speedlim_x, viz_speedlim_y, viz_speedlim_w, viz_speedlim_h, 20);
    if (s->is_ego_over_limit) {
      nvgStrokeColor(s->vg, nvgRGBA(218, 111, 37, 255));
    } else if (is_speedlim_valid) {
      nvgStrokeColor(s->vg, nvgRGBA(255, 255, 255, 255));
    }
  }

  // Draw "Speed Limit" Text
  nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_BASELINE);
  nvgFontFace(s->vg, "sans-semibold");
  nvgFontSize(s->vg, 50);
  nvgFillColor(s->vg, nvgRGBA(0, 0, 0, 255));
  if (is_speedlim_valid && s->is_ego_over_limit) {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
  }
  nvgText(s->vg, viz_speedlim_x+viz_speedlim_w/2 + (is_speedlim_valid ? 6 : 0), viz_speedlim_y + (is_speedlim_valid ? 50 : 45), "SMART", NULL);
  nvgText(s->vg, viz_speedlim_x+viz_speedlim_w/2 + (is_speedlim_valid ? 6 : 0), viz_speedlim_y + (is_speedlim_valid ? 90 : 85), "SPEED", NULL);

  // Draw Speed Text
  nvgFontFace(s->vg, "sans-bold");
  nvgFontSize(s->vg, 48*2.5);
  if (s->is_ego_over_limit) {
    nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
  } else {
    nvgFillColor(s->vg, nvgRGBA(0, 0, 0, 255));
  }
  if (is_speedlim_valid) {
    snprintf(speedlim_str, sizeof(speedlim_str), "%d", speedlim_calc);
    nvgText(s->vg, viz_speedlim_x+viz_speedlim_w/2, viz_speedlim_y + (is_speedlim_valid ? 170 : 165), speedlim_str, NULL);
  } else {
    nvgFontFace(s->vg, "sans-semibold");
    nvgFontSize(s->vg, 42*2.5);
    nvgText(s->vg, viz_speedlim_x+viz_speedlim_w/2, viz_speedlim_y + (is_speedlim_valid ? 170 : 165), "N/A", NULL);
  }
}

static void ui_draw_vision_speed(UIState *s) {
  if (s->b.tri_state_switch == 3) {
    return;
  }
  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;
  float speed = s->scene.v_ego;

  const int viz_speed_w = 280;
  const int viz_speed_x = ui_viz_rx+((ui_viz_rw/2)-(viz_speed_w/2));
  char speed_str[32];

  nvgBeginPath(s->vg);
  nvgRect(s->vg, viz_speed_x, box_y, viz_speed_w, header_h);
  nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_BASELINE);

  if (s->is_metric) {
    snprintf(speed_str, sizeof(speed_str), "%d", (int)(speed * 3.6 + 0.5));
  } else {
    snprintf(speed_str, sizeof(speed_str), "%d", (int)(speed * 2.2369363 + 0.5));
  }
  nvgFontFace(s->vg, "sans-bold");
  nvgFontSize(s->vg, 96*2.5);
  nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
  nvgText(s->vg, viz_speed_x+viz_speed_w/2, 240, speed_str, NULL);

  nvgFontFace(s->vg, "sans-regular");
  nvgFontSize(s->vg, 36*2.5);
  nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 200));

  if (s->is_metric) {
    nvgText(s->vg, viz_speed_x+viz_speed_w/2, 320, "kph", NULL);
  } else {
    nvgText(s->vg, viz_speed_x+viz_speed_w/2, 320, "mph", NULL);
  }
}

static void ui_draw_vision_event(UIState *s) {
  const UIScene *scene = &s->scene;
  const int ui_viz_rx = scene->ui_viz_rx;
  const int ui_viz_rw = scene->ui_viz_rw;
  const int viz_event_w = 220;
  const int viz_event_x = ((ui_viz_rx + ui_viz_rw) - (viz_event_w + (bdr_s*2)));
  const int viz_event_y = (box_y + (bdr_s*1.5));
  const int viz_event_h = (header_h - (bdr_s*1.5));
  if (s->scene.decel_for_model && s->scene.engaged) {
    // draw winding road sign
    const int img_turn_size = 160*1.5;
    const int img_turn_x = viz_event_x-(img_turn_size/4);
    const int img_turn_y = viz_event_y+bdr_s-25;
    float img_turn_alpha = 1.0f;
    nvgBeginPath(s->vg);
    NVGpaint imgPaint = nvgImagePattern(s->vg, img_turn_x, img_turn_y,
      img_turn_size, img_turn_size, 0, s->img_turn, img_turn_alpha);
    nvgRect(s->vg, img_turn_x, img_turn_y, img_turn_size, img_turn_size);
    nvgFillPaint(s->vg, imgPaint);
    nvgFill(s->vg);
  } else {
    // draw steering wheel
    const int bg_wheel_size = 96;
    const int bg_wheel_x = viz_event_x + (viz_event_w-bg_wheel_size);
    const int bg_wheel_y = viz_event_y + (bg_wheel_size/2);
    const int img_wheel_size = bg_wheel_size*1.5;
    const int img_wheel_x = bg_wheel_x-(img_wheel_size/2);
    const int img_wheel_y = bg_wheel_y-25;
    float img_wheel_alpha = 0.1f;
    bool is_engaged = (s->status == STATUS_ENGAGED);
    bool is_warning = (s->status == STATUS_WARNING);
    bool is_engageable = scene->engageable;
    if (is_engaged || is_warning || is_engageable) {
      nvgBeginPath(s->vg);
      nvgCircle(s->vg, bg_wheel_x, (bg_wheel_y + (bdr_s*1.5)), bg_wheel_size);
      if (is_engaged) {
        nvgFillColor(s->vg, nvgRGBA(23, 134, 68, 255));
      } else if (is_warning) {
        nvgFillColor(s->vg, nvgRGBA(218, 111, 37, 255));
      } else if (is_engageable) {
        nvgFillColor(s->vg, nvgRGBA(23, 51, 73, 255));
      }
      nvgFill(s->vg);
      img_wheel_alpha = 1.0f;
    }
    nvgBeginPath(s->vg);
    NVGpaint imgPaint = nvgImagePattern(s->vg, img_wheel_x, img_wheel_y,
      img_wheel_size, img_wheel_size, 0, s->img_wheel, img_wheel_alpha);
    nvgRect(s->vg, img_wheel_x, img_wheel_y, img_wheel_size, img_wheel_size);
    nvgFillPaint(s->vg, imgPaint);
    nvgFill(s->vg);
  }
}

static void ui_draw_vision_map(UIState *s) {
  const UIScene *scene = &s->scene;
  const int map_size = 96;
  const int map_x = (scene->ui_viz_rx + (map_size * 3) + (bdr_s * 3));
  const int map_y = (footer_y + ((footer_h - map_size) / 2));
  const int map_img_size = (map_size * 1.5);
  const int map_img_x = (map_x - (map_img_size / 2));
  const int map_img_y = (map_y - (map_size / 4));

  bool map_valid = s->scene.map_valid;
  float map_img_alpha = map_valid ? 1.0f : 0.15f;
  float map_bg_alpha = map_valid ? 0.3f : 0.1f;
  NVGcolor map_bg = nvgRGBA(0, 0, 0, (255 * map_bg_alpha));
  NVGpaint map_img = nvgImagePattern(s->vg, map_img_x, map_img_y,
    map_img_size, map_img_size, 0, s->img_map, map_img_alpha);

  nvgBeginPath(s->vg);
  nvgCircle(s->vg, map_x, (map_y + (bdr_s * 1.5)), map_size);
  nvgFillColor(s->vg, map_bg);
  nvgFill(s->vg);

  nvgBeginPath(s->vg);
  nvgRect(s->vg, map_img_x, map_img_y, map_img_size, map_img_size);
  nvgFillPaint(s->vg, map_img);
  nvgFill(s->vg);
}

static void ui_draw_vision_face(UIState *s) {
  const UIScene *scene = &s->scene;
  const int face_size = 96;
  const int face_x = (scene->ui_viz_rx + face_size + (bdr_s * 2));
  const int face_y = (footer_y + ((footer_h - face_size) / 2));
  const int face_img_size = (face_size * 1.5);
  const int face_img_x = (face_x - (face_img_size / 2));
  const int face_img_y = (face_y - (face_size / 4));
  float face_img_alpha = scene->monitoring_active ? 1.0f : 0.15f;
  float face_bg_alpha = scene->monitoring_active ? 0.3f : 0.1f;
  NVGcolor face_bg = nvgRGBA(0, 0, 0, (255 * face_bg_alpha));
  NVGpaint face_img = nvgImagePattern(s->vg, face_img_x, face_img_y,
    face_img_size, face_img_size, 0, s->img_face, face_img_alpha);

  nvgBeginPath(s->vg);
  nvgCircle(s->vg, face_x, (face_y + (bdr_s * 1.5)), face_size);
  nvgFillColor(s->vg, face_bg);
  nvgFill(s->vg);

  nvgBeginPath(s->vg);
  nvgRect(s->vg, face_img_x, face_img_y, face_img_size, face_img_size);
  nvgFillPaint(s->vg, face_img);
  nvgFill(s->vg);
}

static void ui_draw_vision_header(UIState *s) {
  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;

  nvgBeginPath(s->vg);
  NVGpaint gradient = nvgLinearGradient(s->vg, ui_viz_rx,
                        (box_y+(header_h-(header_h/2.5))),
                        ui_viz_rx, box_y+header_h,
                        nvgRGBAf(0,0,0,0.45), nvgRGBAf(0,0,0,0));
  nvgFillPaint(s->vg, gradient);
  nvgRect(s->vg, ui_viz_rx, box_y, ui_viz_rw, header_h);
  nvgFill(s->vg);

  ui_draw_vision_maxspeed(s);
  if (s->b.tri_state_switch == 2) {
    ui_draw_vision_speedlimit(s);
  }
  ui_draw_vision_speed(s);
  ui_draw_vision_event(s);
}

static void ui_draw_vision_footer(UIState *s) {
  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;

  nvgBeginPath(s->vg);
  nvgRect(s->vg, ui_viz_rx, footer_y, ui_viz_rw, footer_h);

  ui_draw_vision_face(s);

#ifdef SHOW_SPEEDLIMIT
  // ui_draw_vision_map(s);
#endif
	
  ui_draw_infobar(s);
}



static void ui_draw_vision_alert(UIState *s, int va_size, int va_color,
                                  const char* va_text1, const char* va_text2) {
  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;
  bool hasSidebar = !s->scene.uilayout_sidebarcollapsed;
  bool mapEnabled = s->scene.uilayout_mapenabled;
  bool longAlert1 = strlen(va_text1) > 15;

  const uint8_t *color = alert_colors[va_color];
  const int alr_s = alert_sizes[va_size];
  const int alr_x = ui_viz_rx-(mapEnabled?(hasSidebar?nav_w:(nav_ww)):0)-bdr_s;
  const int alr_w = ui_viz_rw+(mapEnabled?(hasSidebar?nav_w:(nav_ww)):0)+(bdr_s*2);
  const int alr_h = alr_s+(va_size==ALERTSIZE_NONE?0:bdr_s);
  const int alr_y = vwp_h-alr_h;

  nvgBeginPath(s->vg);
  nvgRect(s->vg, alr_x, alr_y, alr_w, alr_h);
  nvgFillColor(s->vg, nvgRGBA(color[0],color[1],color[2],(color[3]*s->alert_blinking_alpha)));
  nvgFill(s->vg);

  nvgBeginPath(s->vg);
  NVGpaint gradient = nvgLinearGradient(s->vg, alr_x, alr_y, alr_x, alr_y+alr_h,
                        nvgRGBAf(0.0,0.0,0.0,0.05), nvgRGBAf(0.0,0.0,0.0,0.35));
  nvgFillPaint(s->vg, gradient);
  nvgRect(s->vg, alr_x, alr_y, alr_w, alr_h);
  nvgFill(s->vg);

  nvgFillColor(s->vg, nvgRGBA(255, 255, 255, 255));
  nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_BASELINE);

  if (va_size == ALERTSIZE_SMALL) {
    nvgFontFace(s->vg, "sans-semibold");
    nvgFontSize(s->vg, 40*2.5);
    nvgText(s->vg, alr_x+alr_w/2, alr_y+alr_h/2+15, va_text1, NULL);
  } else if (va_size== ALERTSIZE_MID) {
    nvgFontFace(s->vg, "sans-bold");
    nvgFontSize(s->vg, 48*2.5);
    nvgText(s->vg, alr_x+alr_w/2, alr_y+alr_h/2-45, va_text1, NULL);
    nvgFontFace(s->vg, "sans-regular");
    nvgFontSize(s->vg, 36*2.5);
    nvgText(s->vg, alr_x+alr_w/2, alr_y+alr_h/2+75, va_text2, NULL);
  } else if (va_size== ALERTSIZE_FULL) {
    nvgFontSize(s->vg, (longAlert1?72:96)*2.5);
    nvgFontFace(s->vg, "sans-bold");
    nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
    nvgTextBox(s->vg, alr_x, alr_y+(longAlert1?360:420), alr_w-60, va_text1, NULL);
    nvgFontSize(s->vg, 48*2.5);
    nvgFontFace(s->vg, "sans-regular");
    nvgTextAlign(s->vg, NVG_ALIGN_CENTER | NVG_ALIGN_BOTTOM);
    nvgTextBox(s->vg, alr_x, alr_h-(longAlert1?300:360), alr_w-60, va_text2, NULL);
  }
  
}


static void ui_draw_vision(UIState *s) {
  //BB code added to only draw every other frame
  if (!s->b.shouldDrawFrame) {
    s->b.shouldDrawFrame = true;
    //return;
  }
  s->b.shouldDrawFrame = false;
  //BBEND

  const UIScene *scene = &s->scene;
  int ui_viz_rx = scene->ui_viz_rx;
  int ui_viz_rw = scene->ui_viz_rw;
  int ui_viz_ro = scene->ui_viz_ro;

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_STENCIL_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  // Draw video frames
  glEnable(GL_SCISSOR_TEST);
  glViewport(ui_viz_rx+ui_viz_ro, s->fb_h-(box_y+box_h), viz_w, box_h);
  glScissor(ui_viz_rx, s->fb_h-(box_y+box_h), ui_viz_rw, box_h);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  draw_frame(s);
  glViewport(0, 0, s->fb_w, s->fb_h);
  glDisable(GL_SCISSOR_TEST);

  glClear(GL_STENCIL_BUFFER_BIT);

  nvgBeginFrame(s->vg, s->fb_w, s->fb_h, 1.0f);
  nvgSave(s->vg);

  // Draw augmented elements
  const int inner_height = viz_w*9/16;
  nvgScissor(s->vg, ui_viz_rx, box_y, ui_viz_rw, box_h);
  nvgTranslate(s->vg, ui_viz_rx+ui_viz_ro, box_y + (box_h-inner_height)/2.0);
  nvgScale(s->vg, (float)viz_w / s->fb_w, (float)inner_height / s->fb_h);
  if (!scene->frontview && !scene->fullview) {
    ui_draw_world(s);
  }

  nvgRestore(s->vg);

  // Set Speed, Current Speed, Status/Events
  ui_draw_vision_header(s);

  if (s->scene.alert_size != ALERTSIZE_NONE) {
    // Controls Alerts
    ui_draw_vision_alert(s, s->scene.alert_size, s->status,
                            s->scene.alert_text1, s->scene.alert_text2);
  } else {
    ui_draw_vision_footer(s);
  }


  nvgEndFrame(s->vg);
  glDisable(GL_BLEND);
}

static void ui_draw_blank(UIState *s) {
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClear(GL_STENCIL_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
}

static void ui_draw(UIState *s) {
  if (s->vision_connected && s->plus_state == 0) {
    ui_draw_vision(s);
  } else {
    ui_draw_blank(s);
  }

  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClear(GL_STENCIL_BUFFER_BIT);

    nvgBeginFrame(s->vg, s->fb_w, s->fb_h, 1.0f);

    nvgEndFrame(s->vg);
    glDisable(GL_BLEND);
  }
}

static PathData read_path(cereal_ModelData_PathData_ptr pathp) {
  PathData ret = {0};

  struct cereal_ModelData_PathData pathd;
  cereal_read_ModelData_PathData(&pathd, pathp);

  ret.prob = pathd.prob;
  ret.std = pathd.std;

  capn_list32 polyp = pathd.poly;
  capn_resolve(&polyp.p);
  for (int i = 0; i < POLYFIT_DEGREE; i++) {
    ret.poly[i] = capn_to_f32(capn_get32(polyp, i));
  }

  // Compute points locations
  for (int i = 0; i < MODEL_PATH_DISTANCE; i++) {
    ret.points[i] = ret.poly[0] * (i*i*i) + ret.poly[1] * (i*i)+ ret.poly[2] * i + ret.poly[3];
  }

  return ret;
}

static ModelData read_model(cereal_ModelData_ptr modelp) {
  struct cereal_ModelData modeld;
  cereal_read_ModelData(&modeld, modelp);

  ModelData d = {0};

  d.path = read_path(modeld.path);
  d.left_lane = read_path(modeld.leftLane);
  d.right_lane = read_path(modeld.rightLane);

  struct cereal_ModelData_LeadData leadd;
  cereal_read_ModelData_LeadData(&leadd, modeld.lead);
  d.lead = (LeadData){
      .dist = leadd.dist, .prob = leadd.prob, .std = leadd.std,
  };

  return d;
}

static void update_status(UIState *s, int status) {
  //BB Variable for the old status
  int old_status = s->status;
  if (s->status != status) {
    s->status = status;
    set_awake(s, true);
    // wake up bg thread to change
    pthread_cond_signal(&s->bg_cond);
  }
}


void handle_message(UIState *s, Message * msg) {
  struct capn ctx;
  capn_init_mem(&ctx, (uint8_t*)msg->getData(), msg->getSize(), 0);

  cereal_Event_ptr eventp;
  eventp.p = capn_getp(capn_root(&ctx), 0, 1);
  struct cereal_Event eventd;
  cereal_read_Event(&eventd, eventp);
  double t = millis_since_boot();

  int bts = bb_get_button_status(s,(char *)"sound");
  if (eventd.which == cereal_Event_controlsState) {
    struct cereal_ControlsState datad;
    cereal_read_ControlsState(&datad, eventd.controlsState);

    s->controls_timeout = 1 * UI_FREQ;
    s->controls_seen = true;

    if (datad.vCruise != s->scene.v_cruise) {
      s->scene.v_cruise_update_ts = eventd.logMonoTime;
    }
    s->scene.v_cruise = datad.vCruise;
    s->scene.v_ego = datad.vEgo;
    s->scene.curvature = datad.curvature;
    //BB get angles
    s->b.angleSteers = datad.angleSteers;
    s->b.angleSteersDes = datad.angleSteersDes;
    //BB END
    s->scene.engaged = datad.enabled;
    s->scene.engageable = datad.engageable;
    s->scene.gps_planner_active = datad.gpsPlannerActive;
    s->scene.monitoring_active = datad.driverMonitoringOn;

    s->scene.frontview = datad.rearViewCam;

    s->scene.decel_for_model = datad.decelForModel;

    if (datad.alertSound != cereal_CarControl_HUDControl_AudibleAlert_none && datad.alertSound != s->alert_sound) {
      if (s->alert_sound != cereal_CarControl_HUDControl_AudibleAlert_none) {
        stop_alert_sound(s->alert_sound);
      }
      if (bts !=0) {
        play_alert_sound(datad.alertSound);
      }

      s->alert_sound = datad.alertSound;
      snprintf(s->alert_type, sizeof(s->alert_type), "%s", datad.alertType.str);
    } else if ((!datad.alertSound || datad.alertSound == cereal_CarControl_HUDControl_AudibleAlert_none)
                  && s->alert_sound != cereal_CarControl_HUDControl_AudibleAlert_none) {
      stop_alert_sound(s->alert_sound);
      s->alert_type[0] = '\0';
      s->alert_sound = cereal_CarControl_HUDControl_AudibleAlert_none;
    }

    if (datad.alertText1.str) {
      snprintf(s->scene.alert_text1, sizeof(s->scene.alert_text1), "%s", datad.alertText1.str);
    } else {
      s->scene.alert_text1[0] = '\0';
    }
    if (datad.alertText2.str) {
      snprintf(s->scene.alert_text2, sizeof(s->scene.alert_text2), "%s", datad.alertText2.str);
    } else {
      s->scene.alert_text2[0] = '\0';
    }
    s->scene.awareness_status = datad.awarenessStatus;

    s->scene.alert_ts = eventd.logMonoTime;

    s->scene.alert_size = datad.alertSize;
    if (datad.alertSize == cereal_ControlsState_AlertSize_none) {
      s->alert_size = ALERTSIZE_NONE;
    } else if (datad.alertSize == cereal_ControlsState_AlertSize_small) {
      s->alert_size = ALERTSIZE_SMALL;
    } else if (datad.alertSize == cereal_ControlsState_AlertSize_mid) {
      s->alert_size = ALERTSIZE_MID;
    } else if (datad.alertSize == cereal_ControlsState_AlertSize_full) {
      s->alert_size = ALERTSIZE_FULL;
    }

    if (datad.alertStatus == cereal_ControlsState_AlertStatus_userPrompt) {
      update_status(s, STATUS_WARNING);
    } else if (datad.alertStatus == cereal_ControlsState_AlertStatus_critical) {
      update_status(s, STATUS_ALERT);
    } else if (datad.enabled) {
      update_status(s, STATUS_ENGAGED);
    } else {
      update_status(s, STATUS_DISENGAGED);
    }

    s->scene.alert_blinkingrate = datad.alertBlinkingRate;
    if (datad.alertBlinkingRate > 0.) {
      if (s->alert_blinked) {
        if (s->alert_blinking_alpha > 0.0 && s->alert_blinking_alpha < 1.0) {
          s->alert_blinking_alpha += (0.05*datad.alertBlinkingRate);
        } else {
          s->alert_blinked = false;
        }
      } else {
        if (s->alert_blinking_alpha > 0.25) {
          s->alert_blinking_alpha -= (0.05*datad.alertBlinkingRate);
        } else {
          s->alert_blinking_alpha += 0.25;
          s->alert_blinked = true;
        }
      }
    }
  } else if (eventd.which == cereal_Event_radarState) {
    struct cereal_RadarState datad;
    cereal_read_RadarState(&datad, eventd.radarState);
    struct cereal_RadarState_LeadData leaddatad;
    cereal_read_RadarState_LeadData(&leaddatad, datad.leadOne);
    s->scene.lead_status = leaddatad.status;
    s->scene.lead_d_rel = leaddatad.dRel;
    s->scene.lead_y_rel = leaddatad.yRel;
    s->scene.lead_v_rel = leaddatad.vRel;
    s->livempc_or_radarstate_changed = true;
  } else if (eventd.which == cereal_Event_liveCalibration) {
    s->scene.world_objects_visible = true;
    struct cereal_LiveCalibrationData datad;
    cereal_read_LiveCalibrationData(&datad, eventd.liveCalibration);

    capn_list32 extrinsicl = datad.extrinsicMatrix;
    capn_resolve(&extrinsicl.p);  // is this a bug?
    for (int i = 0; i < 3 * 4; i++) {
      s->scene.extrinsic_matrix.v[i] =
          capn_to_f32(capn_get32(extrinsicl, i));
    }
  } else if (eventd.which == cereal_Event_model) {
    s->scene.model_ts = eventd.logMonoTime;
    s->scene.model = read_model(eventd.model);
    s->model_changed = true;
  } else if (eventd.which == cereal_Event_liveMpc) {
    struct cereal_LiveMpcData datad;
    cereal_read_LiveMpcData(&datad, eventd.liveMpc);

    capn_list32 x_list = datad.x;
    capn_resolve(&x_list.p);

    for (int i = 0; i < 50; i++){
      s->scene.mpc_x[i] = capn_to_f32(capn_get32(x_list, i));
    }

    capn_list32 y_list = datad.y;
    capn_resolve(&y_list.p);

    for (int i = 0; i < 50; i++){
      s->scene.mpc_y[i] = capn_to_f32(capn_get32(y_list, i));
    }
    s->livempc_or_radarstate_changed = true;
  } else if (eventd.which == cereal_Event_thermal) {
    struct cereal_ThermalData datad;
    cereal_read_ThermalData(&datad, eventd.thermal);

    if (!datad.started) {
      update_status(s, STATUS_STOPPED);
    } else if (s->status == STATUS_STOPPED) {
      // car is started but controls doesn't have fingerprint yet
      update_status(s, STATUS_DISENGAGED);
    }

    s->scene.started_ts = datad.startedTs;

    //BB CPU TEMP
    s->b.maxCpuTemp=datad.cpu0;
    if (s->b.maxCpuTemp<datad.cpu1)
    {
        s->b.maxCpuTemp=datad.cpu1;
    }
    else if (s->b.maxCpuTemp<datad.cpu2)
    {
        s->b.maxCpuTemp=datad.cpu2;
    }
    else if (s->b.maxCpuTemp<datad.cpu3)
    {
        s->b.maxCpuTemp=datad.cpu3;
    }
    s->b.maxBatTemp=datad.bat;
    s->b.freeSpace=datad.freeSpace;
    s->b.batteryPercent=datad.batteryPercent;
    s->b.chargingEnabled=!datad.chargingDisabled;
    s->b.fanSpeed=datad.fanSpeed;
    //BB END CPU TEMP
  } else if (eventd.which == cereal_Event_uiLayoutState) {
    struct cereal_UiLayoutState datad;
    cereal_read_UiLayoutState(&datad, eventd.uiLayoutState);
    s->scene.uilayout_sidebarcollapsed = datad.sidebarCollapsed;
    s->scene.uilayout_mapenabled = datad.mapEnabled;

    bool hasSidebar = !s->scene.uilayout_sidebarcollapsed;
    bool mapEnabled = s->scene.uilayout_mapenabled;
    if (mapEnabled) {
      s->scene.ui_viz_rx = hasSidebar ? (box_x+nav_w) : (box_x+nav_w-(bdr_s*4));
      s->scene.ui_viz_rw = hasSidebar ? (box_w-nav_w) : (box_w-nav_w+(bdr_s*4));
      s->scene.ui_viz_ro = -(sbr_w + 4*bdr_s);
    } else {
      s->scene.ui_viz_rx = hasSidebar ? box_x : (box_x-sbr_w+bdr_s*2);
      s->scene.ui_viz_rw = hasSidebar ? box_w : (box_w+sbr_w-(bdr_s*2));
      s->scene.ui_viz_ro = hasSidebar ? -(sbr_w - 6*bdr_s) : 0;
    }
  } else if (eventd.which == cereal_Event_liveMapData) {
    struct cereal_LiveMapData datad;
    cereal_read_LiveMapData(&datad, eventd.liveMapData);
    s->scene.map_valid = datad.mapValid;
  }
  capn_free(&ctx);
}

static void ui_update(UIState *s) {
  int err;

  if (s->vision_connect_firstrun) {
    // cant run this in connector thread because opengl.
    // do this here for now in lieu of a run_on_main_thread event

    for (int i=0; i<UI_BUF_COUNT; i++) {
      if(s->khr[i] != NULL) {
        visionimg_destroy_gl(s->khr[i], s->priv_hnds[i]);
        glDeleteTextures(1, &s->frame_texs[i]);
      }

      VisionImg img = {
        .fd = s->bufs[i].fd,
        .format = VISIONIMG_FORMAT_RGB24,
        .width = s->rgb_width,
        .height = s->rgb_height,
        .stride = s->rgb_stride,
        .bpp = 3,
        .size = s->rgb_buf_len,
      };
      s->frame_texs[i] = visionimg_to_gl(&img, &s->khr[i], &s->priv_hnds[i]);

      glBindTexture(GL_TEXTURE_2D, s->frame_texs[i]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      // BGR
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_BLUE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_GREEN);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    }

    for (int i=0; i<UI_BUF_COUNT; i++) {
      if(s->khr_front[i] != NULL) {
        visionimg_destroy_gl(s->khr_front[i], s->priv_hnds_front[i]);
        glDeleteTextures(1, &s->frame_front_texs[i]);
      }

      VisionImg img = {
        .fd = s->front_bufs[i].fd,
        .format = VISIONIMG_FORMAT_RGB24,
        .width = s->rgb_front_width,
        .height = s->rgb_front_height,
        .stride = s->rgb_front_stride,
        .bpp = 3,
        .size = s->rgb_front_buf_len,
      };
      s->frame_front_texs[i] = visionimg_to_gl(&img, &s->khr_front[i], &s->priv_hnds_front[i]);

      glBindTexture(GL_TEXTURE_2D, s->frame_front_texs[i]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      // BGR
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_BLUE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_GREEN);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
    }

    assert(glGetError() == GL_NO_ERROR);

    // Default UI Measurements (Assumes sidebar collapsed)
    s->scene.ui_viz_rx = (box_x-sbr_w+bdr_s*2);
    s->scene.ui_viz_rw = (box_w+sbr_w-(bdr_s*2));
    s->scene.ui_viz_ro = 0;

    s->vision_connect_firstrun = false;

    s->alert_blinking_alpha = 1.0;
    s->alert_blinked = false;
  }

  zmq_pollitem_t polls[1] = {{0}};
  // Wait for next rgb image from visiond
  while(true) {
    assert(s->ipc_fd >= 0);
    polls[0].fd = s->ipc_fd;
    polls[0].events = ZMQ_POLLIN;
    int ret = zmq_poll(polls, 1, 1000);
    if (ret < 0) {
      LOGW("poll failed (%d)", ret);
      close(s->ipc_fd);
      s->ipc_fd = -1;
      s->vision_connected = false;
      return;
    } else if (ret == 0)
      continue;
    // vision ipc event
    VisionPacket rp;
    err = vipc_recv(s->ipc_fd, &rp);
    if (err <= 0) {
      LOGW("vision disconnected");
      close(s->ipc_fd);
      s->ipc_fd = -1;
      s->vision_connected = false;
      return;
    }
    if (rp.type == VIPC_STREAM_ACQUIRE) {
      bool front = rp.d.stream_acq.type == VISION_STREAM_RGB_FRONT;
      int idx = rp.d.stream_acq.idx;

      int release_idx;
      if (front) {
        release_idx = s->cur_vision_front_idx;
      } else {
        release_idx = s->cur_vision_idx;
      }
      if (release_idx >= 0) {
        VisionPacket rep = {
          .type = VIPC_STREAM_RELEASE,
          .d = { .stream_rel = {
            .type = rp.d.stream_acq.type,
            .idx = release_idx,
          }},
        };
        vipc_send(s->ipc_fd, &rep);
      }

      if (front) {
        assert(idx < UI_BUF_COUNT);
        s->cur_vision_front_idx = idx;
      } else {
        assert(idx < UI_BUF_COUNT);
        s->cur_vision_idx = idx;
        // printf("v %d\n", ((uint8_t*)s->bufs[idx].addr)[0]);
      }
    } else {
      assert(false);
    }
    break;
  }
  // peek and consume all events in the zmq queue, then return.
  while(true) {
    bool awake = false;
    auto polls = s->poller->poll(0);

    if (polls.size() == 0)
      return;
      
    for (auto sock : polls){
      Message * msg = sock->receive();

      if (sock != s->thermal_sock){
        awake = true;
      }

      if (sock == s->plus_sock){
        s->plus_state = msg->getData()[0];
      } else {
        handle_message(s, msg);
      }

      delete msg;
    }

    if (((awake) && (s->b.tri_state_switch != 3) && (!s->b.keepEonOff)) || (s->b.recording)){
      set_awake(s, true);
    }
  }
}

static int vision_subscribe(int fd, VisionPacket *rp, VisionStreamType type) {
  int err;
  LOGW("vision_subscribe type:%d", type);

  VisionPacket p1 = {
    .type = VIPC_STREAM_SUBSCRIBE,
    .d = { .stream_sub = { .type = type, .tbuffer = true, }, },
  };
  err = vipc_send(fd, &p1);
  if (err < 0) {
    close(fd);
    return 0;
  }

  do {
    err = vipc_recv(fd, rp);
    if (err <= 0) {
      close(fd);
      return 0;
    }

    // release what we aren't ready for yet
    if (rp->type == VIPC_STREAM_ACQUIRE) {
      VisionPacket rep = {
        .type = VIPC_STREAM_RELEASE,
        .d = { .stream_rel = {
          .type = rp->d.stream_acq.type,
          .idx = rp->d.stream_acq.idx,
        }},
      };
      vipc_send(fd, &rep);
    }
  } while (rp->type != VIPC_STREAM_BUFS || rp->d.stream_bufs.type != type);

  return 1;
}

static void* vision_connect_thread(void *args) {
  int err;
  set_thread_name("vision_connect");

  UIState *s = (UIState*)args;
  while (!do_exit) {
    usleep(100000);
    pthread_mutex_lock(&s->lock);
    bool connected = s->vision_connected;
    pthread_mutex_unlock(&s->lock);
    if (connected) continue;

    int fd = vipc_connect();
    if (fd < 0) continue;

    VisionPacket back_rp, front_rp;
    if (!vision_subscribe(fd, &back_rp, VISION_STREAM_RGB_BACK)) continue;
    if (!vision_subscribe(fd, &front_rp, VISION_STREAM_RGB_FRONT)) continue;

    pthread_mutex_lock(&s->lock);
    assert(!s->vision_connected);
    s->ipc_fd = fd;

    ui_init_vision(s,
                   back_rp.d.stream_bufs, back_rp.num_fds, back_rp.fds,
                   front_rp.d.stream_bufs, front_rp.num_fds, front_rp.fds);

    s->vision_connected = true;
    s->vision_connect_firstrun = true;
    pthread_mutex_unlock(&s->lock);
  }
  return NULL;
}


#include <hardware/sensors.h>
#include <utils/Timers.h>

static void* light_sensor_thread(void *args) {
  int err;
  set_thread_name("light_sensor");

  UIState *s = (UIState*)args;
  s->light_sensor = 0.0;

  struct sensors_poll_device_t* device;
  struct sensors_module_t* module;

  hw_get_module(SENSORS_HARDWARE_MODULE_ID, (hw_module_t const**)&module);
  sensors_open(&module->common, &device);

  // need to do this
  struct sensor_t const* list;
  int count = module->get_sensors_list(module, &list);

  int SENSOR_LIGHT = 7;

  err = device->activate(device, SENSOR_LIGHT, 0);
  if (err != 0) goto fail;
  err = device->activate(device, SENSOR_LIGHT, 1);
  if (err != 0) goto fail;

  device->setDelay(device, SENSOR_LIGHT, ms2ns(100));

  while (!do_exit) {
    static const size_t numEvents = 1;
    sensors_event_t buffer[numEvents];

    int n = device->poll(device, buffer, numEvents);
    if (n < 0) {
      LOG_100("light_sensor_poll failed: %d", n);
    }
    if (n > 0) {
      s->light_sensor = buffer[0].light;
    }
  }

  return NULL;

fail:
  LOGE("LIGHT SENSOR IS MISSING");
  s->light_sensor = 255;
  return NULL;
}


static void* bg_thread(void* args) {
  UIState *s = (UIState*)args;
  set_thread_name("bg");

  EGLDisplay bg_display;
  EGLSurface bg_surface;

  FramebufferState *bg_fb = framebuffer_init("bg", 0x00001000, false,
                              &bg_display, &bg_surface, NULL, NULL);
  assert(bg_fb);

  int bg_status = -1;
  while(!do_exit) {
    pthread_mutex_lock(&s->lock);
    //BB Change of background based on our color
    int actual_status = bb_get_status(s);
    if (bg_status == actual_status) {
      // will always be signaled if it changes?
      pthread_cond_wait(&s->bg_cond, &s->lock);
    }
    bg_status = actual_status;
    //BB End of background color change
    pthread_mutex_unlock(&s->lock);

    assert(bg_status < ARRAYSIZE(bg_colors));
    const uint8_t *color = bg_colors[bg_status];

    glClearColor(color[0]/256.0, color[1]/256.0, color[2]/256.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    eglSwapBuffers(bg_display, bg_surface);
    assert(glGetError() == GL_NO_ERROR);
  }

  return NULL;
}

int is_leon() {
  #define MAXCHAR 1000
  FILE *fp;
  char str[MAXCHAR];
  const char* filename = "/proc/cmdline";

  fp = fopen(filename, "r");
  if (fp == NULL){
    printf("Could not open file %s",filename);
    return 0;
  }
  fgets(str, MAXCHAR, fp);
  fclose(fp);
  return strstr(str, "letv") != NULL;
}



int main(int argc, char* argv[]) {
  int err;
  setpriority(PRIO_PROCESS, 0, -14);

  zsys_handler_set(NULL);
  signal(SIGINT, (sighandler_t)set_do_exit);

  UIState uistate;
  UIState *s = &uistate;
  ui_init(s);
  //BB init our UI
  bb_ui_init(s);

  pthread_t connect_thread_handle;
  err = pthread_create(&connect_thread_handle, NULL,
                       vision_connect_thread, s);
  assert(err == 0);

  pthread_t light_sensor_thread_handle;
  err = pthread_create(&light_sensor_thread_handle, NULL,
                       light_sensor_thread, s);
  assert(err == 0);

  pthread_t bg_thread_handle;
  err = pthread_create(&bg_thread_handle, NULL,
                       bg_thread, s);
  assert(err == 0);
  s->b.touch_last_width = s->scene.ui_viz_rw;
  TouchState touch = {0};
  touch_init(&touch);
  s->touch_fd = touch.fd;

  char* error = NULL;
  ui_sound_init(&error);
  if (error) {
    LOGW(error);
    exit(1);
  }

  // light sensor scaling params
  const int LEON = is_leon();

  const float BRIGHTNESS_B = LEON ? 10.0 : 5.0;
  const float BRIGHTNESS_M = LEON ? 2.6 : 1.3;

  float smooth_brightness = BRIGHTNESS_B;

  const int MIN_VOLUME = LEON ? 12 : 9;
  const int MAX_VOLUME = LEON ? 15 : 12;

  set_volume(s, MIN_VOLUME);
#ifdef DEBUG_FPS
  vipc_t1 = millis_since_boot();
  double t1 = millis_since_boot();
  int draws = 0, old_draws = 0;
#endif //DEBUG_FPS
  while (!do_exit) {
    bool should_swap = false;
    if (!s->vision_connected) {
      // Delay a while to avoid 9% cpu usage while car is not started and user is keeping touching on the screen.
      // Don't hold the lock while sleeping, so that vision_connect_thread have chances to get the lock.
      usleep(30 * 1000);
    }
    pthread_mutex_lock(&s->lock);

    // light sensor is only exposed on EONs
    float clipped_brightness = (s->light_sensor*BRIGHTNESS_M) + BRIGHTNESS_B;
    if (clipped_brightness > 512) clipped_brightness = 512;
    smooth_brightness = clipped_brightness * 0.01 + smooth_brightness * 0.99;
    if (smooth_brightness > 255) smooth_brightness = 255;
    set_brightness(s, (int)smooth_brightness);

    int touched = 0;
    int touch_x = -1, touch_y = -1;
    int dc_touch_x = -1, dc_touch_y = -1;
    s->b.touch_timeout --;
    if (s->b.touch_timeout < 0) {
      s->b.touch_timeout = 0;
    }
    if (!s->vision_connected) {
      // Car is not started, keep in idle state and awake on touch events
      zmq_pollitem_t polls[1] = {{0}};
      polls[0].fd = s->touch_fd;
      polls[0].events = ZMQ_POLLIN;
      int ret = zmq_poll(polls, 1, 0);
      if (ret < 0)
        LOGW("poll failed (%d)", ret);
      else if (ret > 0) {
        // awake on any touch
        touched = touch_read(&touch, &touch_x, &touch_y);
      }
    } else {
      // Car started, fetch a new rgb image from ipc and peek for zmq events.
      touched = touch_poll(&touch, &touch_x, &touch_y, s->awake ? 2 : 500);
      //touched = touch_read(&touch, &touch_x, &touch_y);
      ui_update(s);
      if(!s->vision_connected) {
        // Visiond process is just stopped, force a redraw to make screen blank again.
        ui_draw(s);
        glFinish();
        should_swap = true;
      }
    }
    if (touched == 1) {
      set_awake(s, true);
      s->b.touch_last = true;
      s->b.touch_last_x = touch_x;
      s->b.touch_last_y = touch_y;
      s->b.touch_timeout = touch_timeout;
    }
    //BB check touch
    if ((s->b.touch_last) && (s->b.touch_last_width != s->scene.ui_viz_rw)) {
      s->b.touch_last_width=s->scene.ui_viz_rw;
      bb_handle_ui_touch(s,s->b.touch_last_x,s->b.touch_last_y);
      dc_touch_x = s->b.touch_last_x;
      dc_touch_y = s->b.touch_last_y;
      s->b.touch_last = false;
      s->b.touch_last_x = 0;
      s->b.touch_last_y = 0;
    }
    
    //s->b.touch_last_width = s->scene.ui_viz_rw;
    //BB Update our cereal polls
    bb_ui_poll_update(s);
    

    // manage wakefulness
    if (s->awake_timeout > 0) {
      s->awake_timeout--;
    } else {
      set_awake(s, false);
    }

    if (s->awake && s->vision_connected) {
      dashcam(s, dc_touch_x, dc_touch_y);
      ui_draw(s);
      glFinish();
      should_swap = true;
#ifdef DEBUG_FPS
      draws++;
      double t2 = millis_since_boot();
      const double interval = 30.;
      if(t2 - t1 >= interval * 1000.) {
        printf("ui draw fps: %.2f\n",((double)(draws - old_draws)) / interval) ;
        t1 = t2;
        old_draws = draws;
      }
#endif
    }

    if (s->volume_timeout > 0) {
      s->volume_timeout--;
    } else {
      int volume = fmin(MAX_VOLUME, MIN_VOLUME + s->scene.v_ego / 5);  // up one notch every 5 m/s
      set_volume(s, volume);
    }

    if (s->controls_timeout > 0) {
      s->controls_timeout--;
    } else {
      // stop playing alert sound
      if ((!s->vision_connected || (s->vision_connected && s->alert_sound_timeout == 0)) &&
            s->alert_sound != cereal_CarControl_HUDControl_AudibleAlert_none) {
        stop_alert_sound(s->alert_sound);
        s->alert_sound = cereal_CarControl_HUDControl_AudibleAlert_none;
      }

      // if visiond is still running and controlsState times out, display an alert
      if (s->controls_seen && s->vision_connected && strcmp(s->scene.alert_text2, "Controls Unresponsive") != 0) {
        s->scene.alert_size = ALERTSIZE_FULL;
        update_status(s, STATUS_ALERT);
        snprintf(s->scene.alert_text1, sizeof(s->scene.alert_text1), "%s", "TAKE CONTROL IMMEDIATELY");
        snprintf(s->scene.alert_text2, sizeof(s->scene.alert_text2), "%s", "Controls Unresponsive");
        ui_draw_vision_alert(s, s->scene.alert_size, s->status, s->scene.alert_text1, s->scene.alert_text2);

        s->alert_sound_timeout = 2 * UI_FREQ;

        s->alert_sound = cereal_CarControl_HUDControl_AudibleAlert_chimeWarningRepeat;
        play_alert_sound(s->alert_sound);
      }
      s->alert_sound_timeout--;
      s->controls_seen = false;
    }

    read_param_bool_timeout(&s->is_metric, "IsMetric", &s->is_metric_timeout);
    read_param_bool_timeout(&s->longitudinal_control, "LongitudinalControl", &s->longitudinal_control_timeout);
    read_param_bool_timeout(&s->limit_set_speed, "LimitSetSpeed", &s->limit_set_speed_timeout);
    read_param_float_timeout(&s->speed_lim_off, "SpeedLimitOffset", &s->limit_set_speed_timeout);

    pthread_mutex_unlock(&s->lock);

    // the bg thread needs to be scheduled, so the main thread needs time without the lock
    // safe to do this outside the lock?
    if (should_swap) {
      eglSwapBuffers(s->display, s->surface);
    }
  }

  set_awake(s, true);

  slplay_destroy();

  // wake up bg thread to exit
  pthread_mutex_lock(&s->lock);
  pthread_cond_signal(&s->bg_cond);
  pthread_mutex_unlock(&s->lock);
  err = pthread_join(bg_thread_handle, NULL);
  assert(err == 0);

  err = pthread_join(connect_thread_handle, NULL);
  assert(err == 0);

  return 0;
}
