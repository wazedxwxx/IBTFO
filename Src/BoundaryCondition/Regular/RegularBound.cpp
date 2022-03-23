#include "Boundary.H"
#include "WriteData.H"
#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
using namespace std;
void Boundary(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  const double gamma,
                  double *U_OLD,
                  double *U_NEW,
                  double *XYCOORD,
                  double *SCHEME_IDX)
{
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell)];
            }
        }
    }

    // Upon and Down boundary

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            for (int k = 0; k < num_eq; k++)
            {
                U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(i, 2 * num_ghost_cell - 1, k, N_x + 2 * num_ghost_cell)];
                U_OLD[Index(i, N_y + num_ghost_cell + j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(i, N_y + num_ghost_cell - 1, k, N_x + 2 * num_ghost_cell)];
            }
        }
    }
 // Reflect
/*#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < num_ghost_cell; j++)
        {
            U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = -U_NEW[Index(i, 2 * num_ghost_cell - j - 1, 2, N_x + 2 * num_ghost_cell)];
            U_OLD[Index(i, N_y + num_ghost_cell + j, 2, N_x + 2 * num_ghost_cell)] = -U_NEW[Index(i, N_y + num_ghost_cell - 1 - j, 2, N_x + 2 * num_ghost_cell)];
        }
    }*/

    // Right and Left Boundary

#pragma acc parallel loop
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
#pragma acc loop
        for (int i = 0; i < num_ghost_cell; i++)
            for (int k = 0; k < num_eq; k++)
            {
                U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(2 * num_ghost_cell - 1, j, k, N_x + 2 * num_ghost_cell)];
                U_OLD[Index(N_x + num_ghost_cell + i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(N_x + num_ghost_cell - 1, j, k, N_x + 2 * num_ghost_cell)];
            }
    }
}
