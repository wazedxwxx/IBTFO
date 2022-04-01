#include "Advance.H"
#include "TimeAdvance.H"
#include "Slope_limiter.H"
#include "Riemann_solver.H"
#include "WriteData.H"
#include "Conserve2Flux.H"
#include <cstdlib>
#include "EQDefine.H"
#include "CoordDefine.H"
using namespace std;
void TimeAdvance(const double Psy_L,
                  const double Psy_H,
                  const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  const double gamma,
                  double dt,
                  double *U_OLD,
                  double *F_OLD,
                  double *G_OLD,
                  double *U_TMP,
                  double *U_NEW,
                  double *U_TMPRK,
                  double *F_L,
                  double *F_R,
                  double *G_D,
                  double *G_U,
                  double *U_L,
                  double *U_R,
                  double *U_D,
                  double *U_U,
                  double *XYCOORD){
Advance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, F_OLD, G_OLD, U_TMP, U_NEW, F_L, F_R, G_D, G_U, U_L, U_R, U_D, U_U,XYCOORD);
                  }