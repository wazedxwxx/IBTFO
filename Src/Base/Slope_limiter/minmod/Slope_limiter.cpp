#include "Slope_limiter.H"
#include <math.h>
#include "EQDefine.H"
using namespace std;
void Slope_limiter(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   double *U_TMP)
{
#pragma acc parallel loop
    for (int i = num_ghost_cell - 1; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {

                U_TMP[Index_U_SLOUT(i, j, k)] = 0;
                if (U_TMP[Index_U_SLIN1(i, j, k)] * U_TMP[Index_U_SLIN2(i, j, k)] > 0)
                {
                    double r = U_TMP[Index_U_SLIN1(i, j, k)] / U_TMP[Index_U_SLIN2(i, j, k)];

                    if (r > 1)
                    {
                        U_TMP[Index_U_SLOUT(i, j, k)] = U_TMP[Index_U_SLIN2(i, j, k)];
                    }
                    else
                    {
                        U_TMP[Index_U_SLOUT(i, j, k)] = r * U_TMP[Index_U_SLIN2(i, j, k)];
                    }
                }
            }
        }
    }
}