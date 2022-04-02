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
             double *F_OLD,
             double *G_OLD,
             double *U_TMP,
             double *U_NEW,
             double *F_L,
             double *F_R,
             double *G_D,
             double *G_U,
             double *U_L,
             double *U_R,
             double *U_D,
             double *U_U,
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
                U_L[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] +
                                                                                          0.5 * Slope_limiter(U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i - 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)],
                                                                                                              U_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]);
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
                U_R[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                          0.5 * Slope_limiter(U_OLD[Index(i + 2, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)],
                                                                                                              U_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]);
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
                U_D[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] +
                                                                                          0.5 * Slope_limiter(U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j - 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)],
                                                                                                              U_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]);
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
                U_U[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                          0.5 * Slope_limiter(U_OLD[Index(i, j + 2, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)],
                                                                                                              U_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]);
            }
        }
    }

    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_L, F_L, U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_R, F_R, U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_D, U_TMP, G_D);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_U, U_TMP, G_U);

#pragma acc parallel loop
    for (int i = num_ghost_cell - 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_L[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_L[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + 0.5 *
                                                                                                                                                                        (F_L[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                                                                                         F_R[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) *
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
                U_R[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_R[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + 0.5 *
                                                                                                                                                                        (F_L[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                                                                                         F_R[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) *
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
                U_D[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_D[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + 0.5 *
                                                                                                                                                                        (G_D[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                                                                                         G_U[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) *
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
                U_U[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_U[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + 0.5 *
                                                                                                                                                                        (G_D[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                                                                                         G_U[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) *
                                                                                                                                                                        dt / dx;
            }
        }
    }

    Riemann_solver(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, U_OLD, F_OLD, G_OLD, F_L, F_R, G_D, G_U, U_L, U_R, U_D, U_U, U_TMP);

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] == 0)
                    U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell), N_y + 2 * num_ghost_cell] = U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                                          dt * (F_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - F_OLD[Index(i - 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) / dx -
                                                                                                                          dt * (G_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - G_OLD[Index(i, j - 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) / dy;
            }
        }
    }
}