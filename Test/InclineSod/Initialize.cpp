#include "Initialize.H"
#include <math.h>
#include "ParamReader.H"
#include "EQDefine.H"
#include "CoordDefine.H"
using namespace std;
void Initialize(char *filename,
                const double Psy_L,
                const double Psy_H,
                const int N_x,
                const int N_y,
                const int num_ghost_cell,
                const double gamma,
                double *U_OLD,
                double *U_NEW,
                double *XYCOORD)
{
    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    const double angle = 3.1415926535 * para.get("angle", 30) / 180; // tube angle
    const double rho_L = para.get("rho_L", 1);                       // left side density
    const double u_L = para.get("u_L", 0);                           // left side x-vel
    const double v_L = para.get("v_L", 0);                           // left side y-vel
    const double p_L = para.get("p_L", 1);                           // left side pressure
    const double rho_R = para.get("rho_R", 0.125);                   // right side density
    const double u_R = para.get("u_R", 0);                           // right side x-vel
    const double v_R = para.get("v_R", 0);                           // right side y-vel
    const double p_R = para.get("p_R", 0.1);                         // right side pressure

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5 || XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1.5)
            {
                if (XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)] <
                    -(1 / tan(angle)) * XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)] + 0.5 * (1 / tan(angle)) + 0.5)
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
            }
            else
            {
                U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = -1;
                U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = 0;
                U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = 0;
                U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = -1;
            }
            for (int k = 0; k < num_eq; k++)
            {
                U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)];
            }
        }
    }
}
