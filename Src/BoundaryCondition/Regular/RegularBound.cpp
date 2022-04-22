// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Boundary.H"

void Boundary(const int N_x,
              const int N_y,
              const int num_ghost_cell,
              const double gamma,
              double *U_OLD,
              double *U_NEW,
              double *XYCOORD,
              double *SCHEME_IDX,
              const int ndevices,
              const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq]) \
    present(U_NEW [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_OLD[Index(i, j, k)] = U_NEW[Index(i, j, k)];
                }
            }
        }
    }
}
