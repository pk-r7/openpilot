#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5738963365238091084) {
   out_5738963365238091084[0] = delta_x[0] + nom_x[0];
   out_5738963365238091084[1] = delta_x[1] + nom_x[1];
   out_5738963365238091084[2] = delta_x[2] + nom_x[2];
   out_5738963365238091084[3] = delta_x[3] + nom_x[3];
   out_5738963365238091084[4] = delta_x[4] + nom_x[4];
   out_5738963365238091084[5] = delta_x[5] + nom_x[5];
   out_5738963365238091084[6] = delta_x[6] + nom_x[6];
   out_5738963365238091084[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6185348661339499235) {
   out_6185348661339499235[0] = -nom_x[0] + true_x[0];
   out_6185348661339499235[1] = -nom_x[1] + true_x[1];
   out_6185348661339499235[2] = -nom_x[2] + true_x[2];
   out_6185348661339499235[3] = -nom_x[3] + true_x[3];
   out_6185348661339499235[4] = -nom_x[4] + true_x[4];
   out_6185348661339499235[5] = -nom_x[5] + true_x[5];
   out_6185348661339499235[6] = -nom_x[6] + true_x[6];
   out_6185348661339499235[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2198506872339758253) {
   out_2198506872339758253[0] = 1.0;
   out_2198506872339758253[1] = 0.0;
   out_2198506872339758253[2] = 0.0;
   out_2198506872339758253[3] = 0.0;
   out_2198506872339758253[4] = 0.0;
   out_2198506872339758253[5] = 0.0;
   out_2198506872339758253[6] = 0.0;
   out_2198506872339758253[7] = 0.0;
   out_2198506872339758253[8] = 0.0;
   out_2198506872339758253[9] = 1.0;
   out_2198506872339758253[10] = 0.0;
   out_2198506872339758253[11] = 0.0;
   out_2198506872339758253[12] = 0.0;
   out_2198506872339758253[13] = 0.0;
   out_2198506872339758253[14] = 0.0;
   out_2198506872339758253[15] = 0.0;
   out_2198506872339758253[16] = 0.0;
   out_2198506872339758253[17] = 0.0;
   out_2198506872339758253[18] = 1.0;
   out_2198506872339758253[19] = 0.0;
   out_2198506872339758253[20] = 0.0;
   out_2198506872339758253[21] = 0.0;
   out_2198506872339758253[22] = 0.0;
   out_2198506872339758253[23] = 0.0;
   out_2198506872339758253[24] = 0.0;
   out_2198506872339758253[25] = 0.0;
   out_2198506872339758253[26] = 0.0;
   out_2198506872339758253[27] = 1.0;
   out_2198506872339758253[28] = 0.0;
   out_2198506872339758253[29] = 0.0;
   out_2198506872339758253[30] = 0.0;
   out_2198506872339758253[31] = 0.0;
   out_2198506872339758253[32] = 0.0;
   out_2198506872339758253[33] = 0.0;
   out_2198506872339758253[34] = 0.0;
   out_2198506872339758253[35] = 0.0;
   out_2198506872339758253[36] = 1.0;
   out_2198506872339758253[37] = 0.0;
   out_2198506872339758253[38] = 0.0;
   out_2198506872339758253[39] = 0.0;
   out_2198506872339758253[40] = 0.0;
   out_2198506872339758253[41] = 0.0;
   out_2198506872339758253[42] = 0.0;
   out_2198506872339758253[43] = 0.0;
   out_2198506872339758253[44] = 0.0;
   out_2198506872339758253[45] = 1.0;
   out_2198506872339758253[46] = 0.0;
   out_2198506872339758253[47] = 0.0;
   out_2198506872339758253[48] = 0.0;
   out_2198506872339758253[49] = 0.0;
   out_2198506872339758253[50] = 0.0;
   out_2198506872339758253[51] = 0.0;
   out_2198506872339758253[52] = 0.0;
   out_2198506872339758253[53] = 0.0;
   out_2198506872339758253[54] = 1.0;
   out_2198506872339758253[55] = 0.0;
   out_2198506872339758253[56] = 0.0;
   out_2198506872339758253[57] = 0.0;
   out_2198506872339758253[58] = 0.0;
   out_2198506872339758253[59] = 0.0;
   out_2198506872339758253[60] = 0.0;
   out_2198506872339758253[61] = 0.0;
   out_2198506872339758253[62] = 0.0;
   out_2198506872339758253[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6018312896326125070) {
   out_6018312896326125070[0] = state[0];
   out_6018312896326125070[1] = state[1];
   out_6018312896326125070[2] = state[2];
   out_6018312896326125070[3] = state[3];
   out_6018312896326125070[4] = state[4];
   out_6018312896326125070[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6018312896326125070[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6018312896326125070[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8176163857963979814) {
   out_8176163857963979814[0] = 1;
   out_8176163857963979814[1] = 0;
   out_8176163857963979814[2] = 0;
   out_8176163857963979814[3] = 0;
   out_8176163857963979814[4] = 0;
   out_8176163857963979814[5] = 0;
   out_8176163857963979814[6] = 0;
   out_8176163857963979814[7] = 0;
   out_8176163857963979814[8] = 0;
   out_8176163857963979814[9] = 1;
   out_8176163857963979814[10] = 0;
   out_8176163857963979814[11] = 0;
   out_8176163857963979814[12] = 0;
   out_8176163857963979814[13] = 0;
   out_8176163857963979814[14] = 0;
   out_8176163857963979814[15] = 0;
   out_8176163857963979814[16] = 0;
   out_8176163857963979814[17] = 0;
   out_8176163857963979814[18] = 1;
   out_8176163857963979814[19] = 0;
   out_8176163857963979814[20] = 0;
   out_8176163857963979814[21] = 0;
   out_8176163857963979814[22] = 0;
   out_8176163857963979814[23] = 0;
   out_8176163857963979814[24] = 0;
   out_8176163857963979814[25] = 0;
   out_8176163857963979814[26] = 0;
   out_8176163857963979814[27] = 1;
   out_8176163857963979814[28] = 0;
   out_8176163857963979814[29] = 0;
   out_8176163857963979814[30] = 0;
   out_8176163857963979814[31] = 0;
   out_8176163857963979814[32] = 0;
   out_8176163857963979814[33] = 0;
   out_8176163857963979814[34] = 0;
   out_8176163857963979814[35] = 0;
   out_8176163857963979814[36] = 1;
   out_8176163857963979814[37] = 0;
   out_8176163857963979814[38] = 0;
   out_8176163857963979814[39] = 0;
   out_8176163857963979814[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8176163857963979814[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8176163857963979814[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8176163857963979814[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8176163857963979814[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8176163857963979814[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8176163857963979814[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8176163857963979814[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8176163857963979814[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8176163857963979814[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8176163857963979814[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8176163857963979814[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8176163857963979814[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8176163857963979814[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8176163857963979814[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8176163857963979814[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8176163857963979814[56] = 0;
   out_8176163857963979814[57] = 0;
   out_8176163857963979814[58] = 0;
   out_8176163857963979814[59] = 0;
   out_8176163857963979814[60] = 0;
   out_8176163857963979814[61] = 0;
   out_8176163857963979814[62] = 0;
   out_8176163857963979814[63] = 1;
}
void h_25(double *state, double *unused, double *out_7694878293876356150) {
   out_7694878293876356150[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2104961253626386855) {
   out_2104961253626386855[0] = 0;
   out_2104961253626386855[1] = 0;
   out_2104961253626386855[2] = 0;
   out_2104961253626386855[3] = 0;
   out_2104961253626386855[4] = 0;
   out_2104961253626386855[5] = 0;
   out_2104961253626386855[6] = 1;
   out_2104961253626386855[7] = 0;
}
void h_24(double *state, double *unused, double *out_4391897424201540505) {
   out_4391897424201540505[0] = state[4];
   out_4391897424201540505[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4827567397836225366) {
   out_4827567397836225366[0] = 0;
   out_4827567397836225366[1] = 0;
   out_4827567397836225366[2] = 0;
   out_4827567397836225366[3] = 0;
   out_4827567397836225366[4] = 1;
   out_4827567397836225366[5] = 0;
   out_4827567397836225366[6] = 0;
   out_4827567397836225366[7] = 0;
   out_4827567397836225366[8] = 0;
   out_4827567397836225366[9] = 0;
   out_4827567397836225366[10] = 0;
   out_4827567397836225366[11] = 0;
   out_4827567397836225366[12] = 0;
   out_4827567397836225366[13] = 1;
   out_4827567397836225366[14] = 0;
   out_4827567397836225366[15] = 0;
}
void h_30(double *state, double *unused, double *out_3981030553482844530) {
   out_3981030553482844530[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4716502762485243472) {
   out_4716502762485243472[0] = 0;
   out_4716502762485243472[1] = 0;
   out_4716502762485243472[2] = 0;
   out_4716502762485243472[3] = 0;
   out_4716502762485243472[4] = 1;
   out_4716502762485243472[5] = 0;
   out_4716502762485243472[6] = 0;
   out_4716502762485243472[7] = 0;
}
void h_26(double *state, double *unused, double *out_1244876662352480417) {
   out_1244876662352480417[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7844754770306402535) {
   out_7844754770306402535[0] = 0;
   out_7844754770306402535[1] = 0;
   out_7844754770306402535[2] = 0;
   out_7844754770306402535[3] = 0;
   out_7844754770306402535[4] = 0;
   out_7844754770306402535[5] = 0;
   out_7844754770306402535[6] = 0;
   out_7844754770306402535[7] = 1;
}
void h_27(double *state, double *unused, double *out_3296926212594727138) {
   out_3296926212594727138[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1605727367337500656) {
   out_1605727367337500656[0] = 0;
   out_1605727367337500656[1] = 0;
   out_1605727367337500656[2] = 0;
   out_1605727367337500656[3] = 1;
   out_1605727367337500656[4] = 0;
   out_1605727367337500656[5] = 0;
   out_1605727367337500656[6] = 0;
   out_1605727367337500656[7] = 0;
}
void h_29(double *state, double *unused, double *out_4584346438159322047) {
   out_4584346438159322047[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2904085140005149528) {
   out_2904085140005149528[0] = 0;
   out_2904085140005149528[1] = 1;
   out_2904085140005149528[2] = 0;
   out_2904085140005149528[3] = 0;
   out_2904085140005149528[4] = 0;
   out_2904085140005149528[5] = 0;
   out_2904085140005149528[6] = 0;
   out_2904085140005149528[7] = 0;
}
void h_28(double *state, double *unused, double *out_6589850807251536648) {
   out_6589850807251536648[0] = state[5];
   out_6589850807251536648[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6601344197636831384) {
   out_6601344197636831384[0] = 0;
   out_6601344197636831384[1] = 0;
   out_6601344197636831384[2] = 0;
   out_6601344197636831384[3] = 0;
   out_6601344197636831384[4] = 0;
   out_6601344197636831384[5] = 1;
   out_6601344197636831384[6] = 0;
   out_6601344197636831384[7] = 0;
   out_6601344197636831384[8] = 0;
   out_6601344197636831384[9] = 0;
   out_6601344197636831384[10] = 0;
   out_6601344197636831384[11] = 0;
   out_6601344197636831384[12] = 0;
   out_6601344197636831384[13] = 0;
   out_6601344197636831384[14] = 1;
   out_6601344197636831384[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5738963365238091084) {
  err_fun(nom_x, delta_x, out_5738963365238091084);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6185348661339499235) {
  inv_err_fun(nom_x, true_x, out_6185348661339499235);
}
void car_H_mod_fun(double *state, double *out_2198506872339758253) {
  H_mod_fun(state, out_2198506872339758253);
}
void car_f_fun(double *state, double dt, double *out_6018312896326125070) {
  f_fun(state,  dt, out_6018312896326125070);
}
void car_F_fun(double *state, double dt, double *out_8176163857963979814) {
  F_fun(state,  dt, out_8176163857963979814);
}
void car_h_25(double *state, double *unused, double *out_7694878293876356150) {
  h_25(state, unused, out_7694878293876356150);
}
void car_H_25(double *state, double *unused, double *out_2104961253626386855) {
  H_25(state, unused, out_2104961253626386855);
}
void car_h_24(double *state, double *unused, double *out_4391897424201540505) {
  h_24(state, unused, out_4391897424201540505);
}
void car_H_24(double *state, double *unused, double *out_4827567397836225366) {
  H_24(state, unused, out_4827567397836225366);
}
void car_h_30(double *state, double *unused, double *out_3981030553482844530) {
  h_30(state, unused, out_3981030553482844530);
}
void car_H_30(double *state, double *unused, double *out_4716502762485243472) {
  H_30(state, unused, out_4716502762485243472);
}
void car_h_26(double *state, double *unused, double *out_1244876662352480417) {
  h_26(state, unused, out_1244876662352480417);
}
void car_H_26(double *state, double *unused, double *out_7844754770306402535) {
  H_26(state, unused, out_7844754770306402535);
}
void car_h_27(double *state, double *unused, double *out_3296926212594727138) {
  h_27(state, unused, out_3296926212594727138);
}
void car_H_27(double *state, double *unused, double *out_1605727367337500656) {
  H_27(state, unused, out_1605727367337500656);
}
void car_h_29(double *state, double *unused, double *out_4584346438159322047) {
  h_29(state, unused, out_4584346438159322047);
}
void car_H_29(double *state, double *unused, double *out_2904085140005149528) {
  H_29(state, unused, out_2904085140005149528);
}
void car_h_28(double *state, double *unused, double *out_6589850807251536648) {
  h_28(state, unused, out_6589850807251536648);
}
void car_H_28(double *state, double *unused, double *out_6601344197636831384) {
  H_28(state, unused, out_6601344197636831384);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
};

ekf_init(car);
