#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5738963365238091084);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6185348661339499235);
void car_H_mod_fun(double *state, double *out_2198506872339758253);
void car_f_fun(double *state, double dt, double *out_6018312896326125070);
void car_F_fun(double *state, double dt, double *out_8176163857963979814);
void car_h_25(double *state, double *unused, double *out_7694878293876356150);
void car_H_25(double *state, double *unused, double *out_2104961253626386855);
void car_h_24(double *state, double *unused, double *out_4391897424201540505);
void car_H_24(double *state, double *unused, double *out_4827567397836225366);
void car_h_30(double *state, double *unused, double *out_3981030553482844530);
void car_H_30(double *state, double *unused, double *out_4716502762485243472);
void car_h_26(double *state, double *unused, double *out_1244876662352480417);
void car_H_26(double *state, double *unused, double *out_7844754770306402535);
void car_h_27(double *state, double *unused, double *out_3296926212594727138);
void car_H_27(double *state, double *unused, double *out_1605727367337500656);
void car_h_29(double *state, double *unused, double *out_4584346438159322047);
void car_H_29(double *state, double *unused, double *out_2904085140005149528);
void car_h_28(double *state, double *unused, double *out_6589850807251536648);
void car_H_28(double *state, double *unused, double *out_6601344197636831384);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}