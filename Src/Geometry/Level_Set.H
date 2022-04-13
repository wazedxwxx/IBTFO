// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _LEVELSET_H_
#include "Level_Set_function.H"
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
#include "openacc.h"
using namespace std;
void Level_Set(char *filename,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD,
               const int ndevices,
               const int device);
#define _LEVELSET_H_
#endif
               