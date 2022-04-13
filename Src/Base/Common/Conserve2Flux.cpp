#include "Conserve2Flux.H"

void Conserve2Flux(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *F,
                   double *G,
                   const int ndevices,
                   const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U[:(M) * (upper - lower - 1) * num_tmp_size]) \
    present(F[:(M) * (upper - lower - 1) * num_tmp_size])              \
        present(G[:(M) * (upper - lower - 1) * num_tmp_size])
    {
#pragma acc parallel loop async
        for (int j = 0; j < upper - lower; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho = U[Index_U_TMP(i, j, 0)];
                double u = U[Index_U_TMP(i, j, 1)] / rho;
                double v = U[Index_U_TMP(i, j, 2)] / rho;
                double p = (gamma - 1) * (U[Index_U_TMP(i, j, 3)] - 0.5 * rho * (u * u + v * v));


                F[Index_U_TMP(i, j, 0)] = rho * u;
                F[Index_U_TMP(i, j, 1)] = rho * u * u + p;
                F[Index_U_TMP(i, j, 2)] = rho * u * v;
                F[Index_U_TMP(i, j, 3)] = (U[Index_U_TMP(i, j, 3)] + p) * u;

                G[Index_U_TMP(i, j, 0)] = rho * v;
                G[Index_U_TMP(i, j, 1)] = rho * u * v;
                G[Index_U_TMP(i, j, 2)] = rho * v * v + p;
                G[Index_U_TMP(i, j, 3)] = (U[Index_U_TMP(i, j, 3)] + p) * v;
            }
        }
    }
}

void Conserve2FluxX(const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U,
                    double *F,
                    const int ndevices,
                    const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U[:(M) * (upper - lower - 1) * num_tmp_size]) \
    present(F[:(M) * (upper - lower - 1) * num_tmp_size])
    {
#pragma acc parallel loop async
        for (int j = 0; j < upper - lower; j++)
        {

#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho = U[Index_U_TMP(i, j, 0)];
                double u = U[Index_U_TMP(i, j, 1)] / rho;
                double v = U[Index_U_TMP(i, j, 2)] / rho;
                double p = (gamma - 1) * (U[Index_U_TMP(i, j, 3)] - 0.5 * rho * (u * u + v * v));

                F[Index_U_TMP(i, j, 0)] = rho * u;
                F[Index_U_TMP(i, j, 1)] = rho * u * u + p;
                F[Index_U_TMP(i, j, 2)] = rho * u * v;
                F[Index_U_TMP(i, j, 3)] = (U[Index_U_TMP(i, j, 3)] + p) * u;
            }
        }
    }
}

void Conserve2FluxY(const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U,
                    double *G,
                    const int ndevices,
                    const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U[:(M) * (upper - lower - 1) * num_tmp_size]) \
    present(G[:(M) * (upper - lower - 1) * num_tmp_size])
    {
#pragma acc parallel loop async
        for (int j = 0; j < upper - lower; j++)
        {

#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho = U[Index_U_TMP(i, j, 0)];
                double u = U[Index_U_TMP(i, j, 1)] / rho;
                double v = U[Index_U_TMP(i, j, 2)] / rho;
                double p = (gamma - 1) * (U[Index_U_TMP(i, j, 3)] - 0.5 * rho * (u * u + v * v));

                G[Index_U_TMP(i, j, 0)] = rho * v;
                G[Index_U_TMP(i, j, 1)] = rho * u * v;
                G[Index_U_TMP(i, j, 2)] = rho * v * v + p;
                G[Index_U_TMP(i, j, 3)] = (U[Index_U_TMP(i, j, 3)] + p) * v;
            }
        }
    }
}