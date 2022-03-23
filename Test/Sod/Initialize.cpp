#include "Initialize.H"
#include <math.h>
#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
#define Index_Coord(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
using namespace std;
void Initialize(const double Psy_L,
                const double Psy_H,
                const int N_x,
                const int N_y,
                const int num_ghost_cell,
                const double rho_L,
                const double u_L,
                const double v_L,
                const double p_L,
                const double rho_R,
                const double u_R,
                const double v_R,
                const double p_R,
                const double gamma,
                double *U_OLD,
                double *U_NEW,
                double *XYCOORD)
{
#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            if (j < (N_y + 2 * num_ghost_cell) / 2)
            //   if (XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)] <
            //       -pow(3, 0.5) * XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)] + 0.5 * pow(3, 0.5) + 0.5)
            {
                U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho_L;
                U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho_L * u_L;
                U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho_L * v_L;
                U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p_L / (gamma - 1) + 0.5 * rho_L * (u_L * u_L + v_L * v_L);
            }
            else
            {
                U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho_R;
                U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho_R * u_R;
                U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho_R * v_R;
                U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p_R / (gamma - 1) + 0.5 * rho_R * (u_R * u_R + v_R * v_R);
            }
            for (int k = 0; k < num_eq; k++)
            {
                U_NEW[Index(i, j, i, N_x + 2 * num_ghost_cell)] = U_OLD[Index(i, j, i, N_x + 2 * num_ghost_cell)];
            }
        }
    }
}