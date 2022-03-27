#include "Boundary.H"
#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
#define Index_Coord(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
#define Index_sch(a, b, c, N) ((N) * (b) + (a)) * 8 + (c)
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
                if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] == 0)
                    U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(i, j, k, N_x + 2 * num_ghost_cell)];
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
                if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1) // ghost-cell
                {

                    if (SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] != 0 && SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] != 0)
                    {
                        int IDX = SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)];
                        int IDY = SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)];
                        double rhox = U_OLD[Index(i +IDX, j, 0, N_x + 2 * num_ghost_cell)];
                        double ux = U_OLD[Index(i + IDX, j, 1, N_x + 2 * num_ghost_cell)] / rhox;
                        double vx = U_OLD[Index(i + IDX, j, 2, N_x + 2 * num_ghost_cell)] / rhox;
                        double px = (gamma - 1) * (U_OLD[Index(i + IDX, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rhox * (ux * ux + vx * vx));

                        double rhoy = U_OLD[Index(i, j + IDY, 0, N_x + 2 * num_ghost_cell)];
                        double uy = U_OLD[Index(i, j + IDY, 1, N_x + 2 * num_ghost_cell)] / rhoy;
                        double vy = U_OLD[Index(i, j + IDY, 2, N_x + 2 * num_ghost_cell)] / rhoy;
                        double py = (gamma - 1) * (U_OLD[Index(i, j + IDY, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rhoy * (uy * uy + vy * vy));

                        double rho = 0.5 * (rhox + rhoy);
                        double u = 0.5 * (-ux + uy);
                        double v = 0.5 * (vx - vy);
                        double p = 0.5 * (px + py);
                        U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho;
                        U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho * u;
                        U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho * v;
                        U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                    }
                    else
                    {
                        if (SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] != 0)
                        {
                            int IDX = SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)];
                            int IDY = SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)];
                            double rhox = U_OLD[Index(i + IDX, j, 0, N_x + 2 * num_ghost_cell)];
                            double ux = U_OLD[Index(i + IDX, j, 1, N_x + 2 * num_ghost_cell)] / rhox;
                            double vx = U_OLD[Index(i + IDX, j, 2, N_x + 2 * num_ghost_cell)] / rhox;
                            double px = (gamma - 1) * (U_OLD[Index(i + IDX, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rhox * (ux * ux + vx * vx));
                            double rho = rhox;
                            double u = -ux;
                            double v = vx;
                            double p = px;
                            U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho;
                            U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho * u;
                            U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho * v;
                            U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                        }

                        if (SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] != 0)
                        {
                            int IDX = SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)];
                            int IDY = SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)];
                            double rhoy = U_OLD[Index(i, j + IDY, 0, N_x + 2 * num_ghost_cell)];
                            double uy = U_OLD[Index(i, j + IDY, 1, N_x + 2 * num_ghost_cell)] / rhoy;
                            double vy = U_OLD[Index(i, j + IDY, 2, N_x + 2 * num_ghost_cell)] / rhoy;
                            double py = (gamma - 1) * (U_OLD[Index(i, j + IDY, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rhoy * (uy * uy + vy * vy));

                            double rho = rhoy;
                            double u = uy;
                            double v = -vy;
                            double p = py;
                            U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] = rho;
                            U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] = rho * u;
                            U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] = rho * v;
                            U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] = p / (gamma - 1) + 0.5 * rho * (u * u + v * v);
                        }
                    }
                }
            }
        }
    }

// Right and Left Boundary
#pragma acc parallel loop
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
#pragma acc loop
        for (int i = 0; i < num_ghost_cell; i++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(2 * num_ghost_cell - 1, j, k, N_x + 2 * num_ghost_cell)];
                U_OLD[Index(N_x + num_ghost_cell + i, j, k, N_x + 2 * num_ghost_cell)] = U_NEW[Index(N_x + num_ghost_cell - 1, j, k, N_x + 2 * num_ghost_cell)];
            }
        }
    }
}