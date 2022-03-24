#include "Initialize.H"
#include <math.h>
#include "ParamReader.H"
#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
#define Index_Coord(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
using namespace std;
void Initialize(const double Psy_L,
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
    Params<double> para(DetectParams.open("vortex.inp").numbers());
    const double fa = para.get("fa", 0);
    const double fw = para.get("fw", 0);
    const double f0 = para.get("f0", 0);
    const double f1 = para.get("f1", 1);
    const double rho0 = 1.0;
    const double p0 = 2.0;
    const double xlc = 0.50;
    const double ylc = 0.50;
    const double rvortex = 0.40;
    const double pi = 3.1415926535;
    const double ua = fa * rvortex * pi;
    const double uw = fw * rvortex * pi;
    const double rl0 = f0 * rvortex;
    const double rl1 = f1 * rvortex;
    double ur = 0;
    double p = 0;
    double u = 0;
    double v = 0;

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = -1;
            U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = 0;
            U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = 0;
            U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = -1;
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5 || XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1.5)
            {
                double xp = XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)];
                double yp = XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)];
                double r = pow((xp - xlc) * (xp - xlc) + (yp - ylc) * (yp - ylc), 0.5);
                double phi = atan((yp - ylc) / (xp - xlc));
                if (xp - xlc < 0.0)
                    phi = phi + pi;

                if (r < rvortex)
                {
                    if (r < rvortex / 2.0)
                    {
                        ur = ua * 2.0 * r / rvortex;
                        p = (r / rvortex) * (r / rvortex) + 1.0 - 2.0 * log(2.0);
                    }
                    else
                    {
                        ur = ua * 2.0 * (1.0 - r / rvortex);
                        p = (r / rvortex) * (r / rvortex) + 3.0 - 4.0 * r / rvortex + 2.0 * log(r / rvortex);
                    }
                    double u = -sin(phi) * ur;
                    double v = cos(phi) * ur;
                }
                else
                {
                    p = 0.0;
                    u = 0.0;
                    v = 0.0;
                }

                U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho0;
                U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho0 * u;
                U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho0 * v;
                U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p / (gamma - 1) + 0.5 * rho0 * (u * u + v * v);
            }
        }
    }
