// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "ParamReader.H"
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "openacc.h"
using namespace std;

void Initialize(char *filename,
                const double Psy_L,
                const double Psy_H,
                const int N_x,
                const int N_y,
                const int num_ghost_cell,
                const double gamma,
                double *U_OLD,
                double *U_NEW,
                double *XYCOORD,
               const int ndevices,
               const int device);