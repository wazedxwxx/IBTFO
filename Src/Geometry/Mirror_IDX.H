// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _MIRROR_IDX_H_
#include <math.h>
#include <iostream>
#include "CoordDefine.H"
#include "EQDefine.H"
#include <algorithm>
#include "SchDefine.H"
using namespace std;
void Mirror_IDX(const int N_x,
                const int N_y,
                const int num_ghost_cell,
                double *GFM_Index,
                double *XYCOORD,
                int k_point,
                const int ndevices,
                const int device);
#define _MIRROR_IDX_H_
#endif