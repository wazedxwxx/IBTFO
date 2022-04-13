// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Advance.H"
#include "Slope_limiter.H"
#include "Riemann_solver.H"
#include "WriteData.H"
#include "Conserve2Flux.H"
#include <cstdlib>
#include "EQDefine.H"
#include "CoordDefine.H"
using namespace std;
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
                 const int device);
