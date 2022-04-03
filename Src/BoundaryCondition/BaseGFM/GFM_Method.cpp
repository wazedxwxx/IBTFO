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
            if (XYCOORD[Index_Coord(i, j, 5 )] > 1)
            {
                GFM_Index[Index_sch(i, j, 0 )] = XYCOORD[Index_Coord(i, j, 0 )] +
                                                                          2 * XYCOORD[Index_Coord(i, j, 3 )] *
                                                                              std::abs(XYCOORD[Index_Coord(i, j, 2 )]); // x_Mirror = x_GFM + 2*n_x*Phi

                /*std::cout << " i " << i << " j " << j
                          << " x " << XYCOORD[Index_Coord(i, j, 0 )]
                          << " y " << XYCOORD[Index_Coord(i, j, 1 )]
                          << " Phi " << std::abs(XYCOORD[Index_Coord(i, j, 2 )])
                          << " nx " << XYCOORD[Index_Coord(i, j, 3 )]
                          << " MirrorX " << GFM_Index[Index_sch(i, j, 0 )]
                          << " " << std::endl;*/

                GFM_Index[Index_sch(i, j, 1 )] = XYCOORD[Index_Coord(i, j, 1 )] +
                                                                          2 * XYCOORD[Index_Coord(i, j, 4 )] *
                                                                              std::abs(XYCOORD[Index_Coord(i, j, 2 )]); // y_Mirror = y_GFM + 2*n_y*Phi
            }
            else
            {
                GFM_Index[Index_sch(i, j, 0 )] = 0;
                GFM_Index[Index_sch(i, j, 1 )] = 0;
            }
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5 )] > 1)
            {
                int IDX;
                int IDY;
                double x_mirror = GFM_Index[Index_sch(i, j, 0 )];
                double y_mirror = GFM_Index[Index_sch(i, j, 1 )];
                Mirror_IDX(N_x, N_y, num_ghost_cell, &IDX, &IDY, XYCOORD, x_mirror, y_mirror);
                GFM_Index[Index_sch(i, j, 2 )] = IDX;
                GFM_Index[Index_sch(i, j, 3 )] = IDY;
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
            if (XYCOORD[Index_Coord(i, j, 5 )] > 1)
            {
                double d1 = 0;
                double d2 = 0;
                double d3 = 0;
                double d4 = 0;
                double a1 = 0;
                double a2 = 0;
                double a3 = 0;
                double a4 = 0;
                int IDX = GFM_Index[Index_sch(i, j, 2 )];
                int IDY = GFM_Index[Index_sch(i, j, 3 )];

                if (XYCOORD[Index_Coord(IDX, IDY, 5 )] < 1)
                {
                    d1 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) +
                                   (XYCOORD[Index_Coord(IDX, IDY, 1 )] - GFM_Index[Index_sch(i, j, 1 )]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 1 )] - GFM_Index[Index_sch(i, j, 1 )]));
                    a1 = 1 / d1;
                }
                if (XYCOORD[Index_Coord(IDX, IDY + 1, 5 )] < 1)
                {
                    d2 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY + 1, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) +
                                   (XYCOORD[Index_Coord(IDX, IDY + 1, 1 )] - GFM_Index[Index_sch(i, j, 1 )]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 1 )] - GFM_Index[Index_sch(i, j, 1 )]));
                    a2 = 1 / d2;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5 )] < 1)
                {
                    d3 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1 )] - GFM_Index[Index_sch(i, j, 1 )]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1 )] - GFM_Index[Index_sch(i, j, 1 )]));
                    a3 = 1 / d3;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY, 5 )] < 1)
                {
                    d4 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 0 )] - GFM_Index[Index_sch(i, j, 0 )]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY, 1 )] - GFM_Index[Index_sch(i, j, 1 )]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 1 )] - GFM_Index[Index_sch(i, j, 1 )]));
                    a4 = 1 / d4;
                }

                GFM_Index[Index_sch(i, j, 4 )] = a1;
                GFM_Index[Index_sch(i, j, 5 )] = a2;
                GFM_Index[Index_sch(i, j, 6 )] = a3;
                GFM_Index[Index_sch(i, j, 7 )] = a4;

                /*         std::cout << " i " << i << " j " << j
                          << " x " << XYCOORD[Index_Coord(i, j, 0 )]
                          << " y " << XYCOORD[Index_Coord(i, j, 1 )]
                          << " Phi " << std::abs(XYCOORD[Index_Coord(i, j, 2 )])
                          << " nx " << XYCOORD[Index_Coord(i, j, 3 )]
                          << " ny " << XYCOORD[Index_Coord(i, j, 4 )]
                          << " MirrorX " << GFM_Index[Index_sch(i, j, 0 )]
                          << " MirrorY " << GFM_Index[Index_sch(i, j, 1 )]
                          << " MirrorIDX " << GFM_Index[Index_sch(i, j, 2 )]
                          << " MirrorIDY " << GFM_Index[Index_sch(i, j, 3 )]
                          << " 1/D1 " << a1
                          << " 1/D2 " << a2
                          << " 1/D3 " << a3
                          << " 1/D4 " << a4
                          << " " << std::endl;*/
            }
        }
    }
}