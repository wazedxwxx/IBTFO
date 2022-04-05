#include "Boundary.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
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
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                    U_OLD[Index(i, j, k)] = U_NEW[Index(i, j, k)];
            }
        }
    }

// ghost-cell
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] > 1) // ghost-cell
                {

                    if (SCHEME_IDX[Index_sch(i, j, 0)] != 0 && SCHEME_IDX[Index_sch(i, j, 1)] != 0)
                    {
                        int IDX = SCHEME_IDX[Index_sch(i, j, 0)];
                        int IDY = SCHEME_IDX[Index_sch(i, j, 1)];
                        double rhox = U_OLD[Index(i + IDX, j, 0)];
                        double ux = U_OLD[Index(i + IDX, j, 1)] / rhox;
                        double vx = U_OLD[Index(i + IDX, j, 2)] / rhox;
                        double px = (gamma - 1) * (U_OLD[Index(i + IDX, j, 3)] - 0.5 * rhox * (ux * ux + vx * vx));

                        double rhoy = U_OLD[Index(i, j + IDY, 0)];
                        double uy = U_OLD[Index(i, j + IDY, 1)] / rhoy;
                        double vy = U_OLD[Index(i, j + IDY, 2)] / rhoy;
                        double py = (gamma - 1) * (U_OLD[Index(i, j + IDY, 3)] - 0.5 * rhoy * (uy * uy + vy * vy));

                        double rho = 0.5 * (rhox + rhoy);
                        double u = 0.5 * (-ux + uy);
                        double v = 0.5 * (vx - vy);
                        double p = 0.5 * (px + py);
                        U_OLD[Index(i, j, 0)] = rho;
                        U_OLD[Index(i, j, 1)] = rho * u;
                        U_OLD[Index(i, j, 2)] = rho * v;
                        U_OLD[Index(i, j, 3)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                    }
                    else
                    {
                        if (SCHEME_IDX[Index_sch(i, j, 0)] != 0)
                        {
                            int IDX = SCHEME_IDX[Index_sch(i, j, 0)];
                            int IDY = SCHEME_IDX[Index_sch(i, j, 1)];
                            double rhox = U_OLD[Index(i + IDX, j, 0)];
                            double ux = U_OLD[Index(i + IDX, j, 1)] / rhox;
                            double vx = U_OLD[Index(i + IDX, j, 2)] / rhox;
                            double px = (gamma - 1) * (U_OLD[Index(i + IDX, j, 3)] - 0.5 * rhox * (ux * ux + vx * vx));
                            double rho = rhox;
                            double u = -ux;
                            double v = vx;
                            double p = px;
                            U_OLD[Index(i, j, 0)] = rho;
                            U_OLD[Index(i, j, 1)] = rho * u;
                            U_OLD[Index(i, j, 2)] = rho * v;
                            U_OLD[Index(i, j, 3)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                        }

                        if (SCHEME_IDX[Index_sch(i, j, 1)] != 0)
                        {
                            int IDX = SCHEME_IDX[Index_sch(i, j, 0)];
                            int IDY = SCHEME_IDX[Index_sch(i, j, 1)];
                            double rhoy = U_OLD[Index(i, j + IDY, 0)];
                            double uy = U_OLD[Index(i, j + IDY, 1)] / rhoy;
                            double vy = U_OLD[Index(i, j + IDY, 2)] / rhoy;
                            double py = (gamma - 1) * (U_OLD[Index(i, j + IDY, 3)] - 0.5 * rhoy * (uy * uy + vy * vy));

                            double rho = rhoy;
                            double u = uy;
                            double v = -vy;
                            double p = py;
                            U_OLD[Index(i, j, 0)] = rho;
                            U_OLD[Index(i, j, 1)] = rho * u;
                            U_OLD[Index(i, j, 2)] = rho * v;
                            U_OLD[Index(i, j, 3)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                        }
                    }
                }
            }
        }
    }

}


