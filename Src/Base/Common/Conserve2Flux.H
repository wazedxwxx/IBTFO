// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#ifndef _CONSERVE2FLUX_H_
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
#include <algorithm>
using namespace std;

void Conserve2Flux(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *F,
                   double *G,
                   const int ndevices,
                   const int device);

void Conserve2FluxX(const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U,
                    double *G,
                    const int ndevices,
                    const int device);

void Conserve2FluxY(const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U,
                    double *F,
                    const int ndevices,
                    const int device);

void Conserve2FluxZ(const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U,
                    double *H,
                    const int ndevices,
                    const int device);
#define _CONSERVE2FLUX_H_
#endif





