// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _WRITEDATA_H_
#include "EQDefine.H"
#include <string>
#include "CoordDefine.H"
#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "openacc.h"
using namespace std;
void WriteData(const double lo_x,
               const double lo_y,
               const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               int iter,
               double now_t,
               const double gamma,
               double *U_OLD,
               double *XYCOORD,
               const int ndevices,
               const int device);
#define _WRITEDATA_H_
#endif