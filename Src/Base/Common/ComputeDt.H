// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _COMPUTEDT_H_
#define _COMPUTEDT_H_

#endif

// declare function
#include <algorithm>
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
#include <math.h>
using namespace std;

void ComputeDt(const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               const double gamma,
               const double CFL_number,
               double *U_OLD,
               double *XYCOORD,
               double *dt,
               const int ndevices,
               const int device);
