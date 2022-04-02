#include "Conserve2Flux.H"
#include "EQDefine.H"
#include "CoordDefine.H"
void Conserve2Flux(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *F,
                   double *G)
{
#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            double rho = U[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)];
            double u = U[Index(i, j, 1, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] / U[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)];
            double v = U[Index(i, j, 2, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] / U[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)];
            double p = (gamma - 1) * (U[Index(i, j, 3, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] - 0.5 * rho * (u * u + v * v));


            F[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * u;
            F[Index(i, j, 1, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * u * u + p;
            F[Index(i, j, 2, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * u * v;
            F[Index(i, j, 3, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = (U[Index(i, j, 3, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + p) * u;

            G[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * v;
            G[Index(i, j, 1, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * u * v;
            G[Index(i, j, 2, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = rho * v * v + p;
            G[Index(i, j, 3, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] = (U[Index(i, j, 3, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] + p) * v;

        }
    }
}