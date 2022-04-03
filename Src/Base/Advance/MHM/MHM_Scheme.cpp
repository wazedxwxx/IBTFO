/* MUSCL Hancock Scheme Ref <<Riemann Solver>> Toro P505*/
#include "Advance.H"
#include "Slope_limiter.H"
#include "Riemann_solver.H"
#include "WriteData.H"
#include "Conserve2Flux.H"
#include <cstdlib>
#include "EQDefine.H"
#include "CoordDefine.H"
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
             double *XYCOORD)
{
    const double dx = Psy_L / N_x;
    const double dy = Psy_H / N_y;

#pragma acc parallel loop
    for (int i = num_ghost_cell - 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_L(i, j, k)] = U_OLD[Index(i, j, k)] +
                                            0.5 * Slope_limiter(U_OLD[Index(i, j, k)] - U_OLD[Index(i - 1, j, k)],
                                                                U_OLD[Index(i + 1, j, k)] - U_OLD[Index(i, j, k)]);
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell - 2; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_R(i, j, k)] = U_OLD[Index(i + 1, j, k)] -
                                            0.5 * Slope_limiter(U_OLD[Index(i + 2, j, k)] - U_OLD[Index(i + 1, j, k)],
                                                                U_OLD[Index(i + 1, j, k)] - U_OLD[Index(i, j, k)]);
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell - 1; j < N_y + 2 * num_ghost_cell - 1; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_D(i, j, k)] = U_OLD[Index(i, j, k)] +
                                            0.5 * Slope_limiter(U_OLD[Index(i, j, k)] - U_OLD[Index(i, j - 1, k)],
                                                                U_OLD[Index(i, j + 1, k)] - U_OLD[Index(i, j, k)]);
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell - 2; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_U(i, j, k)] = U_OLD[Index(i, j + 1, k)] -
                                            0.5 * Slope_limiter(U_OLD[Index(i, j + 2, k)] - U_OLD[Index(i, j + 1, k)],
                                                                U_OLD[Index(i, j + 1, k)] - U_OLD[Index(i, j, k)]);
            }
        }
    }

    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_L(0, 0, 0)], &U_TMP[Index_F_L(0, 0, 0)], U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_R(0, 0, 0)], &U_TMP[Index_F_R(0, 0, 0)], U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_D(0, 0, 0)], U_TMP, &U_TMP[Index_G_D(0, 0, 0)]);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_U(0, 0, 0)], U_TMP, &U_TMP[Index_G_U(0, 0, 0)]);

#pragma acc parallel loop
    for (int i = num_ghost_cell - 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_L(i, j, k)] = U_TMP[Index_U_L(i, j, k)] + 0.5 *
                                                                            (U_TMP[Index_F_L(i, j, k)] -
                                                                             U_TMP[Index_F_R(i, j, k)]) *
                                                                            dt / dx;
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell - 2; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_R(i, j, k)] = U_TMP[Index_U_R(i, j, k)] + 0.5 *
                                                                            (U_TMP[Index_F_L(i, j, k)] -
                                                                             U_TMP[Index_F_R(i, j, k)]) *
                                                                            dt / dx;
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell - 1; j < N_y + 2 * num_ghost_cell - 1; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_D(i, j, k)] = U_TMP[Index_U_D(i, j, k)] + 0.5 *
                                                                            (U_TMP[Index_G_D(i, j, k)] -
                                                                             U_TMP[Index_G_U(i, j, k)]) *
                                                                            dt / dx;
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell - 2; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_TMP[Index_U_U(i, j, k)] = U_TMP[Index_U_U(i, j, k)] + 0.5 *
                                                                            (U_TMP[Index_G_D(i, j, k)] -
                                                                             U_TMP[Index_G_U(i, j, k)]) *
                                                                            dt / dx;
            }
        }
    }

    Riemann_solver(Psy_L,
                   Psy_H,
                   N_x,
                   N_y,
                   num_ghost_cell,
                   gamma,
                   U_OLD,
                   U_TMP);

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                    U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)] -
                                            dt * (U_TMP[Index_F_OLD(i, j, k)] - U_TMP[Index_F_OLD(i - 1, j, k)]) / dx -
                                            dt * (U_TMP[Index_G_OLD(i, j, k)] - U_TMP[Index_G_OLD(i, j - 1, k)]) / dy;
            }
        }
    }
}