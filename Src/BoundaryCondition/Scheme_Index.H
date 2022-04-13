// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _SCHEME_INDEX_H_
#include "Mirror_IDX.H"
#include <iostream>
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include <algorithm>
using namespace std;
void Scheme_Index(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *XYCOORD,
                  double *GFM_Index,
                  const int ndevices,
                  const int device);
#define _SCHEME_INDEX_H_
#endif