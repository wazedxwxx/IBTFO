// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Riemann_solver.H"
void Riemann_solver(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device)
{
    long int lower = LOWER;
    long int upper = UPPER;
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size],   \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])   
    {
        Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_L(0, lower, 0)], &U_TMP[Index_F_L(0, lower, 0)], &U_TMP[Index_U_TMP(0, lower, 0)], ndevices,device);
        Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_R(0, lower, 0)], &U_TMP[Index_F_R(0, lower, 0)], &U_TMP[Index_U_TMP(0, lower, 0)], ndevices,device);
        Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_D(0, lower, 0)], &U_TMP[Index_U_TMP(0, lower, 0)], &U_TMP[Index_G_D(0, lower, 0)], ndevices,device);
        Conserve2Flux(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_U(0, lower, 0)], &U_TMP[Index_U_TMP(0, lower, 0)], &U_TMP[Index_G_U(0, lower, 0)], ndevices,device);

#pragma acc parallel loop async
for (int j = lower; j < upper; j++)
        {
#pragma acc loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho_L = U_TMP[Index_U_L(i, j, 0)]; // density
                double rho_R = U_TMP[Index_U_R(i, j, 0)]; // density
                double rho_U = U_TMP[Index_U_U(i, j, 0)]; // density
                double rho_D = U_TMP[Index_U_D(i, j, 0)]; // density

                double u_L = U_TMP[Index_U_L(i, j, 1)] / rho_L;                                                 // x-velocity
                double v_L = U_TMP[Index_U_L(i, j, 2)] / rho_L;                                                 // y-velocity
                double p_L = (gamma - 1) * (U_TMP[Index_U_L(i, j, 3)] - 0.5 * rho_L * (u_L * u_L + v_L * v_L)); // pressure

                double u_R = U_TMP[Index_U_R(i, j, 1)] / rho_R;                                                 // x-velocity
                double v_R = U_TMP[Index_U_R(i, j, 2)] / rho_R;                                                 // y-velocity
                double p_R = (gamma - 1) * (U_TMP[Index_U_R(i, j, 3)] - 0.5 * rho_R * (u_R * u_R + v_R * v_R)); // pressure

                double u_U = U_TMP[Index_U_U(i, j, 1)] / rho_U;                                                 // x-velocity
                double v_U = U_TMP[Index_U_U(i, j, 2)] / rho_U;                                                 // y-velocity
                double p_U = (gamma - 1) * (U_TMP[Index_U_U(i, j, 3)] - 0.5 * rho_U * (u_U * u_U + v_U * v_U)); // pressure

                double u_D = U_TMP[Index_U_D(i, j, 1)] / rho_D;                                                 // x-velocity
                double v_D = U_TMP[Index_U_D(i, j, 2)] / rho_D;                                                 // y-velocity
                double p_D = (gamma - 1) * (U_TMP[Index_U_D(i, j, 3)] - 0.5 * rho_D * (u_D * u_D + v_D * v_D)); // pressure

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
                        U_TMP[Index_F_OLD(i, j, k)] = U_TMP[Index_F_L(i, j, k)];
                    }
                    else
                    {
                        if (S_R <= 0)
                            U_TMP[Index_F_OLD(i, j, k)] = U_TMP[Index_F_R(i, j, k)];
                        else
                            U_TMP[Index_F_OLD(i, j, k)] = (S_R * U_TMP[Index_F_L(i, j, k)] -
                                                           S_L * U_TMP[Index_F_R(i, j, k)] +
                                                           S_R * S_L * (U_TMP[Index_U_R(i, j, k)] - U_TMP[Index_U_L(i, j, k)])) /
                                                          (S_R - S_L);
                    }

                    if (S_D >= 0)
                    {
                    }
                    else
                    {
                        if (S_U <= 0)
                            U_TMP[Index_G_OLD(i, j, k)] = U_TMP[Index_G_U(i, j, k)];
                        else
                            U_TMP[Index_G_OLD(i, j, k)] = (S_U * U_TMP[Index_G_D(i, j, k)] -
                                                           S_D * U_TMP[Index_G_U(i, j, k)] +
                                                           S_U * S_D * (U_TMP[Index_U_U(i, j, k)] - U_TMP[Index_U_D(i, j, k)])) /
                                                          (S_U - S_D);
                    }
                }
            }
        }
    }
}


void Riemann_solverX(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size],   \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])   
    {
        Conserve2FluxX(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_L(0, lower, 0)], &U_TMP[Index_F_L(0, lower, 0)], ndevices,device);
        Conserve2FluxX(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_R(0, lower, 0)], &U_TMP[Index_F_R(0, lower, 0)], ndevices,device);

#pragma acc parallel loop async
for (int j = lower; j < upper; j++)
        {
#pragma acc loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho_L = U_TMP[Index_U_L(i, j, 0)]; // density
                double rho_R = U_TMP[Index_U_R(i, j, 0)]; // density


                double u_L = U_TMP[Index_U_L(i, j, 1)] / rho_L;                                                 // x-velocity
                double v_L = U_TMP[Index_U_L(i, j, 2)] / rho_L;                                                 // y-velocity
                double p_L = (gamma - 1) * (U_TMP[Index_U_L(i, j, 3)] - 0.5 * rho_L * (u_L * u_L + v_L * v_L)); // pressure

                double u_R = U_TMP[Index_U_R(i, j, 1)] / rho_R;                                                 // x-velocity
                double v_R = U_TMP[Index_U_R(i, j, 2)] / rho_R;                                                 // y-velocity
                double p_R = (gamma - 1) * (U_TMP[Index_U_R(i, j, 3)] - 0.5 * rho_R * (u_R * u_R + v_R * v_R)); // pressure


                double a_L = std::pow((gamma * p_L / rho_L), 0.5);
                double a_R = std::pow((gamma * p_R / rho_R), 0.5);

                double rho_bar_LR = 0.5 * (rho_L + rho_R);
                double a_bar_LR = 0.5 * (a_L + a_R);
                double p_star_LR = 0.5 * (p_L + p_R) + 0.5 * (u_L - u_R) * (rho_bar_LR * a_bar_LR);
                double S_L = u_L - a_L * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_L + (gamma - 1) / 2 / gamma);
                double S_R = u_R + a_R * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_R + (gamma - 1) / 2 / gamma);

#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    if (S_L >= 0)
                    {
                        U_TMP[Index_F_OLD(i, j, k)] = U_TMP[Index_F_L(i, j, k)];
                    }
                    else
                    {
                        if (S_R <= 0)
                            U_TMP[Index_F_OLD(i, j, k)] = U_TMP[Index_F_R(i, j, k)];
                        else
                            U_TMP[Index_F_OLD(i, j, k)] = (S_R * U_TMP[Index_F_L(i, j, k)] -
                                                           S_L * U_TMP[Index_F_R(i, j, k)] +
                                                           S_R * S_L * (U_TMP[Index_U_R(i, j, k)] - U_TMP[Index_U_L(i, j, k)])) /
                                                          (S_R - S_L);
                    }
                }
            }
        }
    }
}


void Riemann_solverY(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    const double gamma,
                    double *U_OLD,
                    double *U_TMP,
                    const int ndevices,
                    const int device)
{
    long int lower = LOWER;
    long int upper = UPPER;
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size],   \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])   
    {
        Conserve2FluxY(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_D(0, lower, 0)], &U_TMP[Index_G_D(0, lower, 0)], ndevices,device);
        Conserve2FluxY(N_x, N_y, num_ghost_cell, gamma, &U_TMP[Index_U_U(0, lower, 0)], &U_TMP[Index_G_U(0, lower, 0)], ndevices,device);

#pragma acc parallel loop async
for (int j = lower; j < upper; j++)
        {
#pragma acc loop
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double rho_U = U_TMP[Index_U_U(i, j, 0)]; // density
                double rho_D = U_TMP[Index_U_D(i, j, 0)]; // density

                double u_U = U_TMP[Index_U_U(i, j, 1)] / rho_U;                                                 // x-velocity
                double v_U = U_TMP[Index_U_U(i, j, 2)] / rho_U;                                                 // y-velocity
                double p_U = (gamma - 1) * (U_TMP[Index_U_U(i, j, 3)] - 0.5 * rho_U * (u_U * u_U + v_U * v_U)); // pressure

                double u_D = U_TMP[Index_U_D(i, j, 1)] / rho_D;                                                 // x-velocity
                double v_D = U_TMP[Index_U_D(i, j, 2)] / rho_D;                                                 // y-velocity
                double p_D = (gamma - 1) * (U_TMP[Index_U_D(i, j, 3)] - 0.5 * rho_D * (u_D * u_D + v_D * v_D)); // pressure

                double a_U = std::pow((gamma * p_U / rho_U), 0.5);
                double a_D = std::pow((gamma * p_D / rho_D), 0.5);
                double rho_bar_DU = 0.5 * (rho_U + rho_D);
                double a_bar_DU = 0.5 * (a_U + a_D);
                double p_star_DU = 0.5 * (p_U + p_D) + 0.5 * (v_D - v_U) * (rho_bar_DU * a_bar_DU);
                double S_D = v_D - a_D * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_D + (gamma - 1) / 2 / gamma);
                double S_U = v_U + a_U * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_U + (gamma - 1) / 2 / gamma);

#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    if (S_D >= 0)
                    {
                    }
                    else
                    {
                        if (S_U <= 0)
                            U_TMP[Index_G_OLD(i, j, k)] = U_TMP[Index_G_U(i, j, k)];
                        else
                            U_TMP[Index_G_OLD(i, j, k)] = (S_U * U_TMP[Index_G_D(i, j, k)] -
                                                           S_D * U_TMP[Index_G_U(i, j, k)] +
                                                           S_U * S_D * (U_TMP[Index_U_U(i, j, k)] - U_TMP[Index_U_D(i, j, k)])) /
                                                          (S_U - S_D);
                    }
                }
            }
        }
    }
}