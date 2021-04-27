#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_err_fun(double *nom_x, double *delta_x, double *out_3677269693148534025);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3487540600943507420);
void live_H_mod_fun(double *state, double *out_2423990393477680069);
void live_f_fun(double *state, double dt, double *out_8275836173797875467);
void live_F_fun(double *state, double dt, double *out_6755434369386333780);
void live_h_3(double *state, double *unused, double *out_606102128611392106);
void live_H_3(double *state, double *unused, double *out_7126985997297192107);
void live_h_4(double *state, double *unused, double *out_396318269486243879);
void live_H_4(double *state, double *unused, double *out_2521285908091288524);
void live_h_9(double *state, double *unused, double *out_5904822057432972225);
void live_H_9(double *state, double *unused, double *out_4843791451962144494);
void live_h_10(double *state, double *unused, double *out_2507272995598371956);
void live_H_10(double *state, double *unused, double *out_6232543662739183234);
void live_h_12(double *state, double *unused, double *out_7549342054849768165);
void live_H_12(double *state, double *unused, double *out_4627841419830498093);
void live_h_31(double *state, double *unused, double *out_6596628051343678366);
void live_H_31(double *state, double *unused, double *out_6351109027671446400);
void live_h_32(double *state, double *unused, double *out_3379104654345507049);
void live_H_32(double *state, double *unused, double *out_6157287651965591246);
void live_h_13(double *state, double *unused, double *out_6785141495860880007);
void live_H_13(double *state, double *unused, double *out_2899070337993709067);
void live_h_14(double *state, double *unused, double *out_5904822057432972225);
void live_H_14(double *state, double *unused, double *out_4843791451962144494);
void live_h_19(double *state, double *unused, double *out_7217890003287938590);
void live_H_19(double *state, double *unused, double *out_4771915387660607813);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}