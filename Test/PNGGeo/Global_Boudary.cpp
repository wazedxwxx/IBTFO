#include "Boundary.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include <iostream>
#include <math.h>
void Global_Boundary(const int N_x,
                     const int N_y,
                     const int num_ghost_cell,
                     const double gamma,
                     double *U_OLD,
                     double *U_NEW)
{
    // Upon and Down boundary
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            for (int k = 0; k < num_eq; k++)
            {
                U_OLD[Index(i, j, k)] = U_NEW[Index(i, 2 * num_ghost_cell - 1, k)];
                U_OLD[Index(i, N_y + num_ghost_cell + j, k)] = U_NEW[Index(i, N_y + num_ghost_cell - 1, k)];
            }
        }
    }
    // Reflect
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            U_OLD[Index(i, j, 2)] = -U_NEW[Index(i, 2 * num_ghost_cell - j - 1, 2)];
            U_OLD[Index(i, N_y + num_ghost_cell + j, 2)] = -U_NEW[Index(i, N_y + num_ghost_cell - 1 - j, 2)];
        }
    }

// Right and Left Boundary
#pragma acc parallel loop
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
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

// Right and Left Boundary
#pragma acc parallel loop
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
#pragma acc loop
        for (int i = 0; i < num_ghost_cell; i++)
        {
#pragma acc loop

            U_OLD[Index(i, j, 0)] = 1.0;
            U_OLD[Index(i, j, 1)] = 3.0;
            U_OLD[Index(i, j, 2)] = 0.0;
            U_OLD[Index(i, j, 3)] = 0.71429 / (gamma - 1) + 0.5 * 1.0 * (3.0 * 3.0 + 0.0 * 0.0);

        }
    }
}