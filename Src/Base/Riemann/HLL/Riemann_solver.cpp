#include "Conserve2Flux.H"
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
using namespace std;
void Riemann_solver(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *F_OLD,
                    double *G_OLD,
                    double *F_L,
                    double *F_R,
                    double *G_D,
                    double *G_U,
                    double *U_L,
                    double *U_R,
                    double *U_D,
                    double *U_U,
                    double *U_TMP)
{
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_L, F_L, U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_R, F_R, U_TMP);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_D, U_TMP, G_D);
    Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, U_U, U_TMP, G_U);

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
            double rho_L = U_L[Index(i, j, 0, N_x + 2 * num_ghost_cell)]; // density
            double rho_R = U_R[Index(i, j, 0, N_x + 2 * num_ghost_cell)]; // density
            double rho_U = U_U[Index(i, j, 0, N_x + 2 * num_ghost_cell)]; // density
            double rho_D = U_D[Index(i, j, 0, N_x + 2 * num_ghost_cell)]; // density

/*             rho_L = rho_L > 1e-1 ? rho_L : 1e-1;
            rho_R = rho_R > 1e-1 ? rho_R : 1e-1;
            rho_U = rho_U > 1e-1 ? rho_U : 1e-1;
            rho_D = rho_D > 1e-1 ? rho_D : 1e-1; */

            double u_L = U_L[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / rho_L;                                                 // x-velocity
            double v_L = U_L[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / rho_L;                                                 // y-velocity
            double p_L = (gamma - 1) * (U_L[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho_L * (u_L * u_L + v_L * v_L)); // pressure

            double u_R = U_R[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / rho_R;                                                 // x-velocity
            double v_R = U_R[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / rho_R;                                                 // y-velocity
            double p_R = (gamma - 1) * (U_R[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho_R * (u_R * u_R + v_R * v_R)); // pressure

            double u_U = U_U[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / rho_U;                                                 // x-velocity
            double v_U = U_U[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / rho_U;                                                 // y-velocity
            double p_U = (gamma - 1) * (U_U[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho_U * (u_U * u_U + v_U * v_U)); // pressure

            double u_D = U_D[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / rho_D;                                                 // x-velocity
            double v_D = U_D[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / rho_D;                                                 // y-velocity
            double p_D = (gamma - 1) * (U_D[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho_D * (u_D * u_D + v_D * v_D)); // pressure

/*             p_L = p_L > 1e-1 ? p_L : 1e-1;
            p_R = p_R > 1e-1 ? p_R : 1e-1;
            p_U = p_U > 1e-1 ? p_U : 1e-1;
            p_D = p_D > 1e-1 ? p_D : 1e-1; */

            double a_L = std::pow((gamma * p_L / rho_L), 0.5);
            double a_R = std::pow((gamma * p_R / rho_R), 0.5);
            double a_U = std::pow((gamma * p_U / rho_U), 0.5);
            double a_D = std::pow((gamma * p_D / rho_D), 0.5);
            double rho_bar_LR = 0.5 * (rho_L + rho_R);
            double rho_bar_DU = 0.5 * (rho_U + rho_D);
            double a_bar_LR = 0.5 * (a_L + a_R);
            double a_bar_DU = 0.5 * (a_U + a_D);
            double p_star_LR = 0.5 * (p_L + p_R) + 0.5 * (u_L - u_R) * (rho_bar_LR * a_bar_LR);
            double p_star_DU = 0.5 * (p_U + p_D) + 0.5 * (v_D - v_U) * (rho_bar_DU * a_bar_DU);
            double S_L = u_L - a_L * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_L + (gamma - 1) / 2 / gamma);
            double S_R = u_R + a_R * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_R + (gamma - 1) / 2 / gamma);
            double S_D = v_D - a_D * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_D + (gamma - 1) / 2 / gamma);
            double S_U = v_U + a_U * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_U + (gamma - 1) / 2 / gamma);

#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                if (S_L >= 0)
                {
                    F_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = F_L[Index(i, j, k, N_x + 2 * num_ghost_cell)];
                }
                else
                {
                    if (S_R <= 0)
                        F_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = F_R[Index(i, j, k, N_x + 2 * num_ghost_cell)];
                    else
                        F_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = (S_R * F_L[Index(i, j, k, N_x + 2 * num_ghost_cell)] -
                                                                           S_L * F_R[Index(i, j, k, N_x + 2 * num_ghost_cell)] +
                                                                           S_R * S_L * (U_R[Index(i, j, k, N_x + 2 * num_ghost_cell)] - U_L[Index(i, j, k, N_x + 2 * num_ghost_cell)])) /
                                                                          (S_R - S_L);
                }

                if (S_D >= 0)
                {
                    G_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = G_D[Index(i, j, k, N_x + 2 * num_ghost_cell)];
                }
                else
                {
                    if (S_U <= 0)
                        G_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = G_U[Index(i, j, k, N_x + 2 * num_ghost_cell)];
                    else
                        G_OLD[Index(i, j, k, N_x + 2 * num_ghost_cell)] = (S_U * G_D[Index(i, j, k, N_x + 2 * num_ghost_cell)] -
                                                                           S_D * G_U[Index(i, j, k, N_x + 2 * num_ghost_cell)] +
                                                                           S_U * S_D * (U_U[Index(i, j, k, N_x + 2 * num_ghost_cell)] - U_D[Index(i, j, k, N_x + 2 * num_ghost_cell)])) /
                                                                          (S_U - S_D);
                }
            }
        }
    }
}