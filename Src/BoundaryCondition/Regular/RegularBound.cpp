#include "Boundary.H"
#include "WriteData.H"
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
                U_OLD[Index(i, j, k )] = U_NEW[Index(i, j, k )];
            }
        }
    }
}



