#include "Level_Set.H"
#define Index(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
#include <iostream>
// 0 x_coord 1 y_coord 2 phi 3 n_x 4 n_y 5 cell type
void Level_Set(const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD)
{
#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = Level_Set_function(XYCOORD[Index(i, j, 0, N_x + 2 * num_ghost_cell)],
                                                                                   XYCOORD[Index(i, j, 1, N_x + 2 * num_ghost_cell)]); // Phi
        }
    }

#pragma acc parallel loop
    for (int i = 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 1; j < N_y + 2 * num_ghost_cell - 1; j++)
        {
            XYCOORD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = (XYCOORD[Index(i + 1, j, 2, N_x + 2 * num_ghost_cell)] -
                                                                 XYCOORD[Index(i - 1, j, 2, N_x + 2 * num_ghost_cell)]) /
                                                                ((XYCOORD[Index(i + 1, j, 0, N_x + 2 * num_ghost_cell)] -
                                                                  XYCOORD[Index(i - 1, j, 0, N_x + 2 * num_ghost_cell)])); // N_x
            XYCOORD[Index(i, j, 4, N_x + 2 * num_ghost_cell)] = (XYCOORD[Index(i, j + 1, 2, N_x + 2 * num_ghost_cell)] -
                                                                 XYCOORD[Index(i, j - 1, 2, N_x + 2 * num_ghost_cell)]) /
                                                                ((XYCOORD[Index(i, j + 1, 1, N_x + 2 * num_ghost_cell)] -
                                                                  XYCOORD[Index(i, j - 1, 1, N_x + 2 * num_ghost_cell)])); // N_y
        }
    }

    // Identify each cell type, fluid cell(0), solid cell(1), and ghost cell(2)

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] >= 0)
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 0; // fluid cell
            else
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 1; // solid cell
        }
    }

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] == 0)
            {
                for (int k = 1; k < num_ghost_cell + 1; k++)
                {
                    XYCOORD[Index(i + k, j, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i + k, j, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i - k, j, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i - k, j, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i, j + k, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i, j + k, 5, N_x + 2 * num_ghost_cell)];
                    XYCOORD[Index(i, j - k, 5, N_x + 2 * num_ghost_cell)] = 2 * XYCOORD[Index(i, j - k, 5, N_x + 2 * num_ghost_cell)];
                }
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            if (XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
                XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = N_x + num_ghost_cell; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = N_y + num_ghost_cell; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = N_y + num_ghost_cell; j < N_y + 2 * num_ghost_cell; j++)
        {
            XYCOORD[Index(i, j, 5, N_x + 2 * num_ghost_cell)] = 2; // ghost cell
        }
    }

}