// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include <algorithm>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "openacc.h"
#include <algorithm>
using namespace std;
void Psy_coord(const double lo_x,
               const double lo_y,
               const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD,
               const int ndevices,
               const int device);
