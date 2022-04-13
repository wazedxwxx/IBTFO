// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Slope_limiter.H"
#include "Riemann_solver.H"
#include "WriteData.H"
#include "Conserve2Flux.H"
#include <cstdlib>
#include <iostream>
#include "EQDefine.H"
#include "CoordDefine.H"
using namespace std;
void Advance(const double Psy_L,
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
             const int ndevices,
             const int device);
