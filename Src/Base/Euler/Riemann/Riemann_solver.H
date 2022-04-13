// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Conserve2Flux.H"
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include <algorithm>
using namespace std;

void Riemann_solver(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device);

void Riemann_solverX(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device);

void Riemann_solverY(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device);