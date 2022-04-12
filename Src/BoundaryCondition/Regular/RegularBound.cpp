#include "Boundary.H"
#include "WriteData.H"
#include <algorithm>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
using namespace std;
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
#pragma acc data present(U_OLD[(M) * LOWER * num_eq:(M) * (UPPER-LOWER) * num_eq]) \
                 present(U_NEW[(M) * LOWER * num_eq:(M) * (UPPER-LOWER) * num_eq])
    {
#pragma acc parallel loop async
        for (int j = LOWER + num_ghost_cell; j < UPPER - num_ghost_cell; j++)
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
