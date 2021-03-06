// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Scheme_Index.H"
void Scheme_Index(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *XYCOORD,
                  double *SCHEME_IDX)
{

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {

            SCHEME_IDX[Index_sch(i, j, 0)] = 0; // x-direction index
            SCHEME_IDX[Index_sch(i, j, 1)] = 0; // y-direction index
            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {

                for (int k = 0; k < num_ghost_cell + 1; k++)
                {
                    if (XYCOORD[Index_Coord(i - k, j, 5)] < 1) // left search
                    {
                        SCHEME_IDX[Index_sch(i, j, 0)] = 1 - 2 * k;
                        break;
                    }
                }
                if (SCHEME_IDX[Index_sch(i, j, 0)] == 0)
                {
                    for (int k = 0; k < num_ghost_cell + 1; k++)
                    {
                        if (XYCOORD[Index_Coord(i + k, j, 5)] < 1) // right search
                        {
                            SCHEME_IDX[Index_sch(i, j, 0)] = 2 * k - 1;
                            break;
                        }
                    }
                }

                for (int k = 0; k < num_ghost_cell + 1; k++)
                {
                    if (XYCOORD[Index_Coord(i, j - k, 5)] < 1) // down search
                    {
                        SCHEME_IDX[Index_sch(i, j, 1)] = 1 - 2 * k;
                        break;
                    }
                }
                if (SCHEME_IDX[Index_sch(i, j, 1)] == 0)
                {
                    for (int k = 0; k < num_ghost_cell + 1; k++)
                    {
                        if (XYCOORD[Index_Coord(i, j + k, 5)] < 1) // right search
                        {
                            SCHEME_IDX[Index_sch(i, j, 1)] = 2 * k - 1;
                            break;
                        }
                    }
                }
            }
        }
    }
};


