#include "Advance.H"
#include "TimeAdvance.H"
#include "Slope_limiter.H"
#include "Riemann_solver.H"
#include "WriteData.H"
#include "Conserve2Flux.H"
#include <cstdlib>
#include "EQDefine.H"
#include "CoordDefine.H"
#include "Boundary.H"
using namespace std;
void TimeAdvance(const double Psy_L,
                 const double Psy_H,
                 const int N_x,
                 const int N_y,
                 const int num_ghost_cell,
                 const double gamma,
                 double dt,
                 double *U_OLD,
                 double *U_TMP,
                 double *U_NEW,
                 double *XYCOORD,
                 double *SCHEME_IDX)
{

    Advance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, U_TMP, U_NEW, XYCOORD);
    Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, SCHEME_IDX);
    Global_Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW);
    Advance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, U_TMP, &U_TMP[Index_U_TMPRK(0, 0, 0)], XYCOORD);

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                    U_NEW[Index(i, j, k)] = 0.5 *
                                            (U_OLD[Index(i, j, k)] +
                                             U_TMP[Index_U_TMPRK(i, j, k)]);
                else
                    U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)];
            }
        }
    }
}