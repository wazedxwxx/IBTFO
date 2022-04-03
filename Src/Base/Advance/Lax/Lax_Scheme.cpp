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
             double *U_TMP,
             double *U_NEW,
             double *XYCOORD)
{

    const double dx = Psy_L / N_x;
    const double dy = Psy_H / N_y;

    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_OLD, &U_TMP[Index_F_OLD(0, 0, 0)], &U_TMP[Index_G_OLD(0, 0, 0)]);

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
                    U_TMP[Index_U_TMP(i, j, k)] = 0.25 * (U_OLD[Index(i - 1, j, k)] +
                                                          U_OLD[Index(i + 1, j, k)] +
                                                          U_OLD[Index(i, j + 1, k)] +
                                                          U_OLD[Index(i, j - 1, k)]) -
                                                  0.5 * dt * (U_TMP[Index_F_OLD(i + 1, j, k)] - U_TMP[Index_F_OLD(i - 1, j, k)]) / dx;
                /*std::cout << " i " << i << " j " << j << " k " << k << " TYPE "<<XYCOORD[Index_Coord(i, j, 5 )]
                                                                        << " U_OLD " << U_OLD[Index(i, j, 0 )]
                                                                        << " U_OLD1 " << U_OLD[Index(i+1, j, 0 )]
                                                                        << " U_OLD2 " << U_OLD[Index(i-1, j, 0 )]
                                                                        << " U_OLD3 " << U_OLD[Index(i, j+1, 0 )]
                                                                        << " U_OLD4 " << U_OLD[Index(i, j-1, 0 )]
                                                                        << " F1 " << F_OLD[Index(i + 1, j, 0 )]
                                                                        << " F2 " << F_OLD[Index(i - 1, j, 0 )]
                                                                        << " U_TMP " << U_TMP[Index(i - 1, j, 0 )]<<std::endl;*/
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
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                    U_NEW[Index(i, j, k)] = U_TMP[Index_U_TMP(i, j, k)] -
                                            0.5 * dt * (U_TMP[Index_G_OLD(i, j + 1, k)] - U_TMP[Index_G_OLD(i, j - 1, k)]) / dy;
                /*                 else
                                    U_NEW[Index(i, j, k )] = U_OLD[Index(i, j, k )]; */
            }
        }
    }
}
