#include "Scheme_Index.H"
#include "Mirror_IDX.H"
#include <iostream>
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
void Scheme_Index(const int N_x,
                const int N_y,
                const int num_ghost_cell,
                double *XYCOORD,
                double *GFM_Index)
{
// Calculate Mirror point location
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
            {
                GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] = XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)] +
                                                                          2 * XYCOORD[Index_Coord(i, j, 3, N_x + 2 * num_ghost_cell)] *
                                                                              std::abs(XYCOORD[Index_Coord(i, j, 2, N_x + 2 * num_ghost_cell)]); // x_Mirror = x_GFM + 2*n_x*Phi

                /*std::cout << " i " << i << " j " << j
                          << " x " << XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)]
                          << " y " << XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)]
                          << " Phi " << std::abs(XYCOORD[Index_Coord(i, j, 2, N_x + 2 * num_ghost_cell)])
                          << " nx " << XYCOORD[Index_Coord(i, j, 3, N_x + 2 * num_ghost_cell)]
                          << " MirrorX " << GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]
                          << " " << std::endl;*/

                GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] = XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)] +
                                                                          2 * XYCOORD[Index_Coord(i, j, 4, N_x + 2 * num_ghost_cell)] *
                                                                              std::abs(XYCOORD[Index_Coord(i, j, 2, N_x + 2 * num_ghost_cell)]); // y_Mirror = y_GFM + 2*n_y*Phi
            }
            else
            {
                GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] = 0;
                GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] = 0;
            }
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
            {
                int IDX;
                int IDY;
                double x_mirror = GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)];
                double y_mirror = GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)];
                Mirror_IDX(N_x, N_y, num_ghost_cell, &IDX, &IDY, XYCOORD, x_mirror, y_mirror);
                GFM_Index[Index_sch(i, j, 2, N_x + 2 * num_ghost_cell)] = IDX;
                GFM_Index[Index_sch(i, j, 3, N_x + 2 * num_ghost_cell)] = IDY;
            }
        }
    }

// Calculate Inverse distance - Distance
//           i,j+1     i+1,j+1             a2=1/d2      a3=1/d3
//                                =>
//           i,j        i+1,j              a1=1/d1      a4=1/d4
//
//
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
            {
                double d1 = 0;
                double d2 = 0;
                double d3 = 0;
                double d4 = 0;
                double a1 = 0;
                double a2 = 0;
                double a3 = 0;
                double a4 = 0;
                int IDX = GFM_Index[Index_sch(i, j, 2, N_x + 2 * num_ghost_cell)];
                int IDY = GFM_Index[Index_sch(i, j, 3, N_x + 2 * num_ghost_cell)];

                if (XYCOORD[Index_Coord(IDX, IDY, 5, N_x + 2 * num_ghost_cell)] < 1)
                {
                    d1 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]));
                    a1 = 1 / d1;
                }
                if (XYCOORD[Index_Coord(IDX, IDY + 1, 5, N_x + 2 * num_ghost_cell)] < 1)
                {
                    d2 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY + 1, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY + 1, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]));
                    a2 = 1 / d2;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5, N_x + 2 * num_ghost_cell)] < 1)
                {
                    d3 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]));
                    a3 = 1 / d3;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY, 5, N_x + 2 * num_ghost_cell)] < 1)
                {
                    d4 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 0, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 1, N_x + 2 * num_ghost_cell)] - GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]));
                    a4 = 1 / d4;
                }

                GFM_Index[Index_sch(i, j, 4, N_x + 2 * num_ghost_cell)] = a1;
                GFM_Index[Index_sch(i, j, 5, N_x + 2 * num_ghost_cell)] = a2;
                GFM_Index[Index_sch(i, j, 6, N_x + 2 * num_ghost_cell)] = a3;
                GFM_Index[Index_sch(i, j, 7, N_x + 2 * num_ghost_cell)] = a4;

                /*         std::cout << " i " << i << " j " << j
                          << " x " << XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)]
                          << " y " << XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)]
                          << " Phi " << std::abs(XYCOORD[Index_Coord(i, j, 2, N_x + 2 * num_ghost_cell)])
                          << " nx " << XYCOORD[Index_Coord(i, j, 3, N_x + 2 * num_ghost_cell)]
                          << " ny " << XYCOORD[Index_Coord(i, j, 4, N_x + 2 * num_ghost_cell)]
                          << " MirrorX " << GFM_Index[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]
                          << " MirrorY " << GFM_Index[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]
                          << " MirrorIDX " << GFM_Index[Index_sch(i, j, 2, N_x + 2 * num_ghost_cell)]
                          << " MirrorIDY " << GFM_Index[Index_sch(i, j, 3, N_x + 2 * num_ghost_cell)]
                          << " 1/D1 " << a1
                          << " 1/D2 " << a2
                          << " 1/D3 " << a3
                          << " 1/D4 " << a4
                          << " " << std::endl;*/
            }
        }
    }
}