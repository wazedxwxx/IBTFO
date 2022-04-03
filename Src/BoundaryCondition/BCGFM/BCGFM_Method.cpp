#include "Scheme_Index.H"
#include "Mirror_IDX.H"
#include <iostream>
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include "ParamReader.H"
void Scheme_Index(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *XYCOORD,
                  double *GFM_Index)
{

    /*
    Index Define
    0:  Image point x coord
    1:  Image point y coord
    2:  Image point neighbor point x index
    3:  Image point neighbor point y index
    4:  Image point inverse distance coffcient
    5:  Image point inverse distance coffcient
    6:  Image point inverse distance coffcient
    7:  Image point inverse distance coffcient
    8:  Extra point x coord
    9:  Extra point y coord
    10: Extra point neighbor point x index
    11: Extra point neighbor point y index
    12: Extra point inverse distance coffcient
    13: Extra point inverse distance coffcient
    14: Extra point inverse distance coffcient
    15: Extra point inverse distance coffcient

    */

    const double dx = XYCOORD[Index_Coord(1, 0, 0)] - XYCOORD[Index_Coord(0, 0, 0)];
    const double dy = XYCOORD[Index_Coord(0, 1, 0)] - XYCOORD[Index_Coord(0, 0, 0)];

    const double delta = sqrt(dx * dx + dy * dy) / 2;

// Calculate Mirror point location
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {

                // Image point location
                GFM_Index[Index_sch(i, j, 0)] = XYCOORD[Index_Coord(i, j, 0)] +
                                                XYCOORD[Index_Coord(i, j, 3)] *
                                                    (std::abs(XYCOORD[Index_Coord(i, j, 2)]) + delta); // x_Image = x_GFM + n_x*(Phi+0.5*norm(dx,dy))

                GFM_Index[Index_sch(i, j, 1)] = XYCOORD[Index_Coord(i, j, 1)] +
                                                XYCOORD[Index_Coord(i, j, 4)] *
                                                    (std::abs(XYCOORD[Index_Coord(i, j, 2)]) + delta); // y_Image = y_GFM + n_y*(Phi+0.5*norm(dx,dy)

                GFM_Index[Index_sch(i, j, 8)] = XYCOORD[Index_Coord(i, j, 0)] +
                                                XYCOORD[Index_Coord(i, j, 3)] *
                                                    (std::abs(XYCOORD[Index_Coord(i, j, 2)]) + 3.0 * delta); // x_Image = x_GFM + n_x*(Phi+Phi+0.5*norm(dx,dy))

                GFM_Index[Index_sch(i, j, 9)] = XYCOORD[Index_Coord(i, j, 1)] +
                                                XYCOORD[Index_Coord(i, j, 4)] *
                                                    (std::abs(XYCOORD[Index_Coord(i, j, 2)]) + 3.0 * delta); // y_Image = y_GFM + n_y*(Phi+Phi+0.5*norm(dx,dy)
            }
            else
            {
                GFM_Index[Index_sch(i, j, 0)] = 0;
                GFM_Index[Index_sch(i, j, 1)] = 0;
                GFM_Index[Index_sch(i, j, 8)] = 0;
                GFM_Index[Index_sch(i, j, 9)] = 0;
            }
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {
                int IDX;
                int IDY;
                double x_image = GFM_Index[Index_sch(i, j, 0)];
                double y_image = GFM_Index[Index_sch(i, j, 1)];
                double x_extra = GFM_Index[Index_sch(i, j, 8)];
                double y_extra = GFM_Index[Index_sch(i, j, 9)];
                Mirror_IDX(N_x, N_y, num_ghost_cell, &IDX, &IDY, XYCOORD, x_image, y_image);
                GFM_Index[Index_sch(i, j, 2)] = IDX;
                GFM_Index[Index_sch(i, j, 3)] = IDY;

                Mirror_IDX(N_x, N_y, num_ghost_cell, &IDX, &IDY, XYCOORD, x_extra, y_extra);
                GFM_Index[Index_sch(i, j, 10)] = IDX;
                GFM_Index[Index_sch(i, j, 11)] = IDY;
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
            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {
                double d1 = 0;
                double d2 = 0;
                double d3 = 0;
                double d4 = 0;
                double a1 = 0;
                double a2 = 0;
                double a3 = 0;
                double a4 = 0;
                int IDX = GFM_Index[Index_sch(i, j, 2)];
                int IDY = GFM_Index[Index_sch(i, j, 3)];

                if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1)
                {
                    d1 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY, 0)] - GFM_Index[Index_sch(i, j, 0)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 0)] - GFM_Index[Index_sch(i, j, 0)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY, 1)] - GFM_Index[Index_sch(i, j, 1)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 1)] - GFM_Index[Index_sch(i, j, 1)]));
                    a1 = 1 / d1;
                }
                if (XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1)
                {
                    d2 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 0)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 0)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 1)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 1)]));
                    a2 = 1 / d2;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                {
                    d3 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 0)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 0)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 1)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 1)]));
                    a3 = 1 / d3;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1)
                {
                    d4 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY, 0)] - GFM_Index[Index_sch(i, j, 0)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 0)] - GFM_Index[Index_sch(i, j, 0)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY, 1)] - GFM_Index[Index_sch(i, j, 1)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 1)] - GFM_Index[Index_sch(i, j, 1)]));
                    a4 = 1 / d4;
                }

                GFM_Index[Index_sch(i, j, 4)] = a1;
                GFM_Index[Index_sch(i, j, 5)] = a2;
                GFM_Index[Index_sch(i, j, 6)] = a3;
                GFM_Index[Index_sch(i, j, 7)] = a4;
            }
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {
                double d1 = 0;
                double d2 = 0;
                double d3 = 0;
                double d4 = 0;
                double a1 = 0;
                double a2 = 0;
                double a3 = 0;
                double a4 = 0;
                int IDX = GFM_Index[Index_sch(i, j, 10)];
                int IDY = GFM_Index[Index_sch(i, j, 11)];

                if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1)
                {
                    d1 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY, 0)] - GFM_Index[Index_sch(i, j, 8)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 0)] - GFM_Index[Index_sch(i, j, 8)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY, 1)] - GFM_Index[Index_sch(i, j, 9)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY, 1)] - GFM_Index[Index_sch(i, j, 9)]));
                    a1 = 1 / d1;
                }
                if (XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1)
                {
                    d2 = std::sqrt((XYCOORD[Index_Coord(IDX, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 8)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 8)]) +
                                   (XYCOORD[Index_Coord(IDX, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 9)]) *
                                       (XYCOORD[Index_Coord(IDX, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 9)]));
                    a2 = 1 / d2;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                {
                    d3 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 8)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 0)] - GFM_Index[Index_sch(i, j, 8)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 9)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY + 1, 1)] - GFM_Index[Index_sch(i, j, 9)]));
                    a3 = 1 / d3;
                }
                if (XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1)
                {
                    d4 = std::sqrt((XYCOORD[Index_Coord(IDX + 1, IDY, 0)] - GFM_Index[Index_sch(i, j, 8)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 0)] - GFM_Index[Index_sch(i, j, 8)]) +
                                   (XYCOORD[Index_Coord(IDX + 1, IDY, 1)] - GFM_Index[Index_sch(i, j, 9)]) *
                                       (XYCOORD[Index_Coord(IDX + 1, IDY, 1)] - GFM_Index[Index_sch(i, j, 9)]));
                    a4 = 1 / d4;
                }

                GFM_Index[Index_sch(i, j, 12)] = a1;
                GFM_Index[Index_sch(i, j, 13)] = a2;
                GFM_Index[Index_sch(i, j, 14)] = a3;
                GFM_Index[Index_sch(i, j, 15)] = a4;
            }
        }
    }
}