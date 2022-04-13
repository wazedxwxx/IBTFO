// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Mirror_IDX.H"

void Mirror_IDX(const int N_x,
                const int N_y,
                const int num_ghost_cell,
                double *GFM_Index,
                double *XYCOORD,
                int k_point,
                const int ndevices,
                const int device)
{
    int lower = LOWER;
    int upper = UPPER;
// x direction search
#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord], GFM_Index [(M)*lower * num_sch:(M) * (upper - lower) * num_sch])
    {
#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] > 1)
                {
                    double x_mirror = GFM_Index[Index_sch(i, j, k_point)];
#pragma acc loop
                    for (int ii = 0; ii < N_x + 2 * num_ghost_cell - 1; ii++)
                    {
                        if (x_mirror < XYCOORD[Index_Coord(ii + 1, j, 0)] && x_mirror >= XYCOORD[Index_Coord(ii, j, 0)])
                        {
                            GFM_Index[Index_sch(i, j, k_point + 2)] = ii;
                        }
                    }
                }
            }
        }
        


#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] > 1)
                {
                    double y_mirror = GFM_Index[Index_sch(i, j, k_point + 1)];
#pragma acc loop
                    for (int jj = lower; jj < upper - 1; jj++)
                    {
                        if (y_mirror < XYCOORD[Index_Coord(0, jj + 1, 1)] && y_mirror >= XYCOORD[Index_Coord(0, jj, 1)])
                        {
                            GFM_Index[Index_sch(i, j, k_point + 3)] = jj;
                        }
                    }
                }
            }
        }
    }
}
