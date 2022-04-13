#include "Initialize.H"

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
    const double fa = para.get("fa", 0);
    const double fw = para.get("fw", 0);

    const double rho0 = 1.0;
    const double p0 = 2.0;
    const double xlc = 0.50;
    const double ylc = 0.50;
    const double rvortex = 0.40;
    const double pi = M_PI;
    const double ua = fa * rvortex * pi;
    const double uw = fw * rvortex * pi;

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            double ur = 0;
            double p = 0;
            double u = 0;
            double v = 0;
            double xp = 0;
            double yp = 0;
            double r = 0;
            double phi = 0;

            U_OLD[Index(i, j, 0)] = rho0;
            U_OLD[Index(i, j, 1)] = rho0 * u;
            U_OLD[Index(i, j, 2)] = rho0 * v;
            U_OLD[Index(i, j, 3)] = (p0 + 2.0 * rho0 * ua * ua * p) / (gamma - 1) + 0.5 * rho0 * (u * u + v * v);
            if (XYCOORD[Index_Coord(i, j, 5)] < 0.5 || XYCOORD[Index_Coord(i, j, 5)] > 1.5)
            {
                xp = XYCOORD[Index_Coord(i, j, 0)];
                yp = XYCOORD[Index_Coord(i, j, 1)];
                r = pow((xp - xlc) * (xp - xlc) + (yp - ylc) * (yp - ylc), 0.5);
                phi = atan((yp - ylc) / (xp - xlc));

                // cout << " x_coord " << xp << " y_coord " << yp << " phi " << phi << endl;

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
                    u = -sin(phi) * ur;
                    v = cos(phi) * ur;
                }
                else
                {
                    p = 0.0;
                    u = 0.0;
                    v = 0.0;
                }

                U_OLD[Index(i, j, 0)] = rho0;
                U_OLD[Index(i, j, 1)] = rho0 * u;
                U_OLD[Index(i, j, 2)] = rho0 * v;
                U_OLD[Index(i, j, 3)] = (p0 + 2.0 * rho0 * ua * ua * p) / (gamma - 1) + 0.5 * rho0 * (u * u + v * v);
            }
            for (int k = 0; k < num_eq; k++)
            {
                U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)];
            }
        }
    }
}
