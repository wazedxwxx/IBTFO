#include "Advance.H"
#include "Conserve2Flux.H"
#include "WriteData.H"
#include <iostream>
#include "EQDefine.H"
#include "CoordDefine.H"
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

    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_OLD, F_OLD, G_OLD);

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
                    U_TMP[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = 0.25 * (U_OLD[Index(i - 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] +
                                                                                                        U_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] +
                                                                                                        U_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] +
                                                                                                        U_OLD[Index(i, j - 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) -
                                                                                                0.5 * dt * (F_OLD[Index(i + 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - F_OLD[Index(i - 1, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) / dx;
                /*std::cout << " i " << i << " j " << j << " k " << k << " TYPE "<<XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)]
                                                                        << " U_OLD " << U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " U_OLD1 " << U_OLD[Index(i+1, j, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " U_OLD2 " << U_OLD[Index(i-1, j, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " U_OLD3 " << U_OLD[Index(i, j+1, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " U_OLD4 " << U_OLD[Index(i, j-1, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " F1 " << F_OLD[Index(i + 1, j, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " F2 " << F_OLD[Index(i - 1, j, 0, N_x + 2 * num_ghost_cell)]
                                                                        << " U_TMP " << U_TMP[Index(i - 1, j, 0, N_x + 2 * num_ghost_cell)]<<std::endl;*/
            }
        }
    }

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
                    U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_TMP[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] -
                                                                                                0.5 * dt * (G_OLD[Index(i, j + 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - G_OLD[Index(i, j - 1, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)]) / dy;
                else
                    U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)];
            }
        }
    }
}
