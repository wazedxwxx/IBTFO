// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn
/* MUSCL Hancock Scheme Ref <<Riemann Solver>> Toro P505*/
#include "Advance.H"
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
             const int device)
{
    const double dx = Psy_L / N_x;
    const double dy = Psy_H / N_y;
    int lower = LOWER;
    int upper = UPPER;

#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord]) \
    present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size])          \
        present(U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])                  \
            present(U_NEW [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 1; i < N_x + 2 * num_ghost_cell - 1; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_TMP[Index_U_L(i, j, k)] = U_OLD[Index(i, j, k)];
                }
            }
        }

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell - 2; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_TMP[Index_U_R(i, j, k)] = U_OLD[Index(i + 1, j, k)];
                }
            }
        }

        Riemann_solverX(Psy_L,
                        Psy_H,
                        N_x,
                        N_y,
                        num_ghost_cell,
                        gamma,
                        U_OLD,
                        U_TMP,
                        ndevices,
                        device);

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                        U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)] -
                                                dt * (U_TMP[Index_F_OLD(i, j, k)] - U_TMP[Index_F_OLD(i - 1, j, k)]) / dx;
                }
            }
        }

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell - 1; j < upper - 1; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_TMP[Index_U_D(i, j, k)] = U_OLD[Index(i, j, k)];
                }
            }
        }

#pragma acc parallel loop async
        for (int j = lower; j < upper - 2; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_TMP[Index_U_U(i, j, k)] = U_OLD[Index(i, j + 1, k)];
                }
            }
        }

        Riemann_solverY(Psy_L,
                        Psy_H,
                        N_x,
                        N_y,
                        num_ghost_cell,
                        gamma,
                        U_OLD,
                        U_TMP,
                        ndevices,
                        device);

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                        U_NEW[Index(i, j, k)] = U_NEW[Index(i, j, k)] -
                                                dt * (U_TMP[Index_G_OLD(i, j, k)] - U_TMP[Index_G_OLD(i, j - 1, k)]) / dy;
                }
            }
        }
    }
}
