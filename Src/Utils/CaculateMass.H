// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _CACULATEMASS_H_

#include "EQDefine.H"
#include <string>
#include "CoordDefine.H"
#include <iostream>
#include <math.h>
#include "openacc.h"

using namespace std;
void CaculateMass(const double Psy_L,
                  const double Psy_H,
                  const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *U_OLD,
                  double *XYCOORD,
                  double *MASS_DEVICE,
                  const int ndevices,
                  const int device);
#define _CACULATEMASS_H_
#endif