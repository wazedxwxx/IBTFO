// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "TimeAdvance.H"
void TimeAdvance(const double Psy_L,
                 const double Psy_H,
                 const int N_x,
                 const int N_y,
                 const int num_ghost_cell,
                 const double gamma,
                 double dt,
                 double *U_OLD,
                 double *U_TMP,
                 double *U_NEW,
                 double *XYCOORD,
                 double *SCHEME_IDX,
                 const int ndevices,
                 const int device)
{
    Advance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, U_TMP, U_NEW, XYCOORD, ndevices,device);
}