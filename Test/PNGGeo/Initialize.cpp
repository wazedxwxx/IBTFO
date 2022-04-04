#include "Initialize.H"
#include "ParamReader.H"
#include <math.h>
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
    const double rho_1 = para.get("rho_L", 1);     // left side density
    const double u_1 = para.get("u_L", 0);         // left side x-vel
    const double v_1 = para.get("v_L", 0);         // left side y-vel
    const double p_1 = para.get("p_L", 1);         // left side pressure
    const double rho_R = para.get("rho_R", 0.125); // right side density
    const double u_R = para.get("u_R", 0);         // right side x-vel
    const double v_R = para.get("v_R", 0);         // right side y-vel
    const double p_R = para.get("p_R", 0.1);       // right side pressure
#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            if (XYCOORD[Index_Coord(i, j, 0)] > 0.1)
            {
                U_OLD[Index(i, j, 0)] = rho_1;
                U_OLD[Index(i, j, 1)] = rho_1 * u_1;
                U_OLD[Index(i, j, 2)] = rho_1 * v_1;
                U_OLD[Index(i, j, 3)] = p_1 / (gamma - 1) + 0.5 * rho_1 * (u_1 * u_1 + v_1 * v_1);
            }
            else
            {
                double Ms = 1.3;
                double rho_2 = rho_1 * (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms);
                double p_2 = p_1 * ((2 * gamma) * Ms * Ms / (gamma + 1) - (gamma - 1) / (gamma + 1));
                double v_2 = 0;
                double u_2 = Ms * (1 - rho_1 / rho_2);
                U_OLD[Index(i, j, 0)] = rho_2;
                U_OLD[Index(i, j, 1)] = rho_2 * u_2;
                U_OLD[Index(i, j, 2)] = 0;
                U_OLD[Index(i, j, 3)] = p_2 / (gamma - 1) + 0.5 * rho_2 * (u_2 * u_2);
            }
            for (int k = 0; k < num_eq; k++)
            {
                U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)];
            }
        }
    }
}
