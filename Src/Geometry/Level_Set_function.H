// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _LEVELSET_FUNC_H_

#include <math.h>
#include "ParamReader.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
using namespace std;

void Level_Set_function(char *filename,
                        const int N_x,
                        const int N_y,
                        const int num_ghost_cell,
                        double *XYCOORD,
                        const int ndevices,
                        const int device);

#define _LEVELSET_FUNC_H_
#endif
