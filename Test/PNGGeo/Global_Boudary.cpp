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
            double rho1 = 1.0;
            double p1 = 0.71429;
            double Ms = 1.3;

            double rho2 = rho1 * (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms);
            double p2 = p1 * ((2 * gamma) * Ms * Ms / (gamma + 1) - (gamma - 1) / (gamma + 1));
            double v2 = 0;
            double u2 = Ms * (1 - rho1 / rho2);

            U_OLD[Index(i, j, 0)] = rho2;
            U_OLD[Index(i, j, 1)] = rho2*u2;
            U_OLD[Index(i, j, 2)] = 0.0;
            U_OLD[Index(i, j, 3)] = p2 / (gamma - 1) + 0.5 * rho2 * (u2 * u2);
        }
    }
}