// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Level_Set.H"
// 0 x_coord 1 y_coord 2 phi 3 n_x 4 n_y 5 cell type
void Level_Set(char *filename,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD,
               const int ndevices,
               const int device)
{
    int lower = LOWER;
    int upper = UPPER;

#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord])
    {
/* #pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                // XYCOORD[Index_Coord(i, j, 2)] =
                // Level_Set_function(gemo_factor, XYCOORD[Index_Coord(i, j, 0)], XYCOORD[Index_Coord(i, j, 1)]); // Phi
            }
        } */

        Level_Set_function(filename, N_x, N_y, num_ghost_cell, XYCOORD, ndevices, device);

#pragma acc parallel loop async
        for (int j = lower + 1; j < upper - 1; j++)
        {
#pragma acc loop
            for (int i = 1; i < N_x + 2 * num_ghost_cell - 1; i++)
            {
                XYCOORD[Index_Coord(i, j, 3)] = (XYCOORD[Index_Coord(i + 1, j, 2)] -
                                                 XYCOORD[Index_Coord(i - 1, j, 2)]) /
                                                ((XYCOORD[Index_Coord(i + 1, j, 0)] -
                                                  XYCOORD[Index_Coord(i - 1, j, 0)])); // N_x
                XYCOORD[Index_Coord(i, j, 4)] = (XYCOORD[Index_Coord(i, j + 1, 2)] -
                                                 XYCOORD[Index_Coord(i, j - 1, 2)]) /
                                                ((XYCOORD[Index_Coord(i, j + 1, 1)] -
                                                  XYCOORD[Index_Coord(i, j - 1, 1)])); // N_y
            }
        }

        // Identify each cell type, fluid cell(0), solid cell(1), and ghost cell(2)

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 2)] >= 0)
                    XYCOORD[Index_Coord(i, j, 5)] = 0; // fluid cell
                else
                    XYCOORD[Index_Coord(i, j, 5)] = 1; // solid cell
            }
        }

#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                {
#pragma acc loop
                    for (int k = 1; k < num_ghost_cell + 1; k++)
                    {
                        XYCOORD[Index_Coord(i + k, j, 5)] = 2 * XYCOORD[Index_Coord(i + k, j, 5)];
                        XYCOORD[Index_Coord(i - k, j, 5)] = 2 * XYCOORD[Index_Coord(i - k, j, 5)];
                        XYCOORD[Index_Coord(i, j + k, 5)] = 2 * XYCOORD[Index_Coord(i, j + k, 5)];
                        XYCOORD[Index_Coord(i, j - k, 5)] = 2 * XYCOORD[Index_Coord(i, j - k, 5)];
                    }
                }
            }
        }
        acc_async_wait_all();

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] > 1)
                    XYCOORD[Index_Coord(i, j, 5)] = 2; // ghost cell
            }
        }

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                if (i < num_ghost_cell || i >= N_x + num_ghost_cell)
                {
                    XYCOORD[Index_Coord(i, j, 5)] = 2; // ghost cell
                }
            }
        }

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                if (j < num_ghost_cell || j >= N_y + num_ghost_cell)
                {
                    XYCOORD[Index_Coord(i, j, 5)] = 2; // ghost cell
                }
            }
        }

    }


    

    std::cout << " ====  Geometry Initialize complete ====" << std::endl;
}
