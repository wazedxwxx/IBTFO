#include "Boundary.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include <iostream>
#include <math.h>

using namespace std;
void Global_Boundary(const int N_x,
                     const int N_y,
                     const int num_ghost_cell,
                     const double gamma,
                     double *U_OLD,
                     double *U_NEW,
                     const int ndevices,
                     const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U_OLD [(M) * lower * num_eq:(M) * (upper - lower) * num_eq]) \
    present(U_NEW [(M) * lower * num_eq:(M) * (upper - lower) * num_eq])
    {
        // Upon and Down boundary
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                    for (int k = 0; k < num_eq; k++)
                    {
                        if (j < num_ghost_cell)
                        {
                            U_OLD[Index(i, j, k)] = U_NEW[Index(i, 2 * num_ghost_cell - 1, k)];
                        }
                        if (j >= N_y + num_ghost_cell)
                        {
                            U_OLD[Index(i, j, k)] = U_NEW[Index(i, N_y + num_ghost_cell-1, k)];
                        }
                    }
            }
        }
        // Reflect
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (j < num_ghost_cell)
                {
                    U_OLD[Index(i, j, 2)] = -U_NEW[Index(i, 2 * num_ghost_cell - j - 1, 2)];
                }
                if (j >= N_y + num_ghost_cell)
                {
                    U_OLD[Index(i, j , 2)] = -U_NEW[Index(i, 2 * (N_y + num_ghost_cell) - 1 - j, 2)];
                }
            }
        }

// Right and Left Boundary
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    U_OLD[Index(i, j, k)] = U_NEW[Index(num_ghost_cell, j, k)];
                    U_OLD[Index(N_x + num_ghost_cell + i, j, k)] = U_NEW[Index(N_x + num_ghost_cell - 1, j, k)];
                }
            }
        }
    }
}