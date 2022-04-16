// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Scheme_Index.H"
void Scheme_Index(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *XYCOORD,
                  double *GFM_Index,
                  const int ndevices,
                  const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord]) \
    present(GFM_Index [(M)*lower * num_sch:(M) * (upper - lower) * num_sch])
    {

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] > 1)
                {
                    GFM_Index[Index_sch(i, j, 0)] = XYCOORD[Index_Coord(i, j, 0)] +
                                                    2 * XYCOORD[Index_Coord(i, j, 3)] *
                                                        std::abs(XYCOORD[Index_Coord(i, j, 2)]); // x_Mirror = x_GFM + 2*n_x*Phi

                    GFM_Index[Index_sch(i, j, 1)] = XYCOORD[Index_Coord(i, j, 1)] +
                                                    2 * XYCOORD[Index_Coord(i, j, 4)] *
                                                        std::abs(XYCOORD[Index_Coord(i, j, 2)]); // y_Mirror = y_GFM + 2*n_y*Phi
                }
                else
                {
                    GFM_Index[Index_sch(i, j, 0)] = 0;
                    GFM_Index[Index_sch(i, j, 1)] = 0;
                }
            }
        }

  Mirror_IDX(N_x, N_y, num_ghost_cell, GFM_Index, XYCOORD, 0, ndevices, device);


    }
}
