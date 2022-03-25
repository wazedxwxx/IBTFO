#include <math.h>
#include "ComputeDt.H"

#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
#define Index_Coord(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
#include <iostream>
using namespace std;

void ComputeDt(const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               const double gamma,
               const double CFL_number,
               double *U_OLD,
               double *XYCOORD,
               double *dt)
{
    double a_max = 0.0;
    double u_max = 0.0;
    double v_max = 0.0;
    const double dx = Psy_L / N_x;
    const double dy = Psy_H / N_y;
#pragma acc parallel loop reduction(max:a_max) reduction(max:u_max) reduction(max:v_max)
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 1)
            {
                double a = 0.0;
                double rho = U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double u = U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double v = U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double p = (gamma - 1) * (U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho * (u * u + v * v));
                a = std::pow((gamma * p / rho), 0.5);
                a_max = max(a_max, a);
                u_max = max(u_max, abs(u));
                v_max = max(v_max, abs(v));
                
            }
        }
    }
  

    *dt = CFL_number * dx / (a_max +  max(u_max, v_max));
}
