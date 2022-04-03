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
            double rho = U[Index(i, j, 0)];
            double u = U[Index(i, j, 1)] / U[Index(i, j, 0)];
            double v = U[Index(i, j, 2)] / U[Index(i, j, 0)];
            double p = (gamma - 1) * (U[Index(i, j, 3)] - 0.5 * rho * (u * u + v * v));

            F[Index(i, j, 0)] = rho * u;
            F[Index(i, j, 1)] = rho * u * u + p;
            F[Index(i, j, 2)] = rho * u * v;
            F[Index(i, j, 3)] = (U[Index(i, j, 3)] + p) * u;

            G[Index(i, j, 0)] = rho * v;
            G[Index(i, j, 1)] = rho * u * v;
            G[Index(i, j, 2)] = rho * v * v + p;
            G[Index(i, j, 3)] = (U[Index(i, j, 3)] + p) * v;
        }
    }
}