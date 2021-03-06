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
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size], \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
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
                double u_star = 0.5 * (u_L + u_R) + 0.5 * (p_L - p_R) / (rho_bar_LR * a_bar_LR);
                double v_star = 0.5 * (v_D + v_U) + 0.5 * (p_D - p_U) / (rho_bar_DU * a_bar_DU);
                double rhoL_star = rho_L + (p_star_LR - p_L) / a_L / a_L;
                double rhoR_star = rho_R + (p_star_LR - p_R) / a_R / a_R;
                double rhoD_star = rho_D + (p_star_DU - p_D) / a_D / a_D;
                double rhoU_star = rho_U + (p_star_DU - p_U) / a_U / a_U;
                double S_L = u_L - a_L * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_L + (gamma - 1) / 2 / gamma);
                double S_R = u_R + a_R * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_R + (gamma - 1) / 2 / gamma);

                double S_HL = u_L - a_L;
                double S_TL = u_star - pow(a_L * (p_star_LR / p_R), ((gamma - 1) / 2 / gamma));

                double S_HR = u_R + a_R;
                double S_TR = u_star + pow(a_R * (p_star_LR / p_R), ((gamma - 1) / 2 / gamma));

                double S_D = v_D - a_D * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_D + (gamma - 1) / 2 / gamma);
                double S_U = v_U + a_U * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_U + (gamma - 1) / 2 / gamma);

                double S_HD = v_D - a_D;
                double S_TD = v_star - pow(a_D * (p_star_DU / p_U), ((gamma - 1) / 2 / gamma));

                double S_HU = v_U + a_U;
                double S_TU = v_star + pow(a_U * (p_star_DU / p_U), ((gamma - 1) / 2 / gamma));
                if (u_star > 0)
                {
                    if (p_star_LR > p_L)
                    {
                        if (S_L > 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_L * u_L;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_L * u_L * u_L + p_L;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_L * u_L * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_L(i, j, 3)] + p_L) * u_L;
                        }
                        else
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rhoL_star * u_star;
                            U_TMP[Index_F_OLD(i, j, 1)] = rhoL_star * u_star * u_star + p_star_LR;
                            U_TMP[Index_F_OLD(i, j, 2)] = rhoL_star * u_star * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoL_star * (u_star * u_star + v_L * v_L)) * u_star;
                        }
                    }
                    else
                    {
                        if (S_HL > 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_L * u_L;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_L * u_L * u_L + p_L;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_L * u_L * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_L(i, j, 3)] + p_L) * u_L;
                        }
                        else
                        {
                            if (S_TL < 0)
                            {
                                U_TMP[Index_F_OLD(i, j, 0)] = rhoL_star * u_star;
                                U_TMP[Index_F_OLD(i, j, 1)] = rhoL_star * u_star * u_star + p_star_LR;
                                U_TMP[Index_F_OLD(i, j, 2)] = rhoL_star * u_star * v_L;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoL_star * (u_star * u_star + v_L * v_L)) * u_star;
                            }
                            else
                            {
                                double rho_L_fan = rho_L * pow(2 / (gamma + 1) + (gamma - 1) * (rho_L) / (gamma + 1) / a_L, 2 / (gamma - 1));
                                double u_L_fan = 2 * (a_L + (gamma - 1) * u_L / 2) / (gamma + 1);
                                double p_L_fan = p_L * pow((2 / (gamma + 1) + (gamma - 1) * (rho_L) / (gamma + 1) / a_L), 2 * gamma / (gamma - 1));
                                U_TMP[Index_F_OLD(i, j, 0)] = rho_L_fan * u_L_fan;
                                U_TMP[Index_F_OLD(i, j, 1)] = rho_L_fan * u_L_fan * u_L_fan + p_L_fan;
                                U_TMP[Index_F_OLD(i, j, 2)] = rho_L_fan * u_L_fan * v_L;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_L_fan / (gamma - 1) + 0.5 * rho_L_fan * (u_L_fan * u_L_fan + v_L * v_L)) * u_L_fan;
                            }
                        }
                    }
                }
                else
                {
                    if (p_star_LR > p_R)
                    {
                        if (S_R < 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_R * u_R;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_R * u_R * u_R + p_R;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_R * u_R * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_R(i, j, 3)] + p_R) * u_R;
                        }
                        else
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rhoR_star * u_star;
                            U_TMP[Index_F_OLD(i, j, 1)] = rhoR_star * u_star * u_star + p_star_LR;
                            U_TMP[Index_F_OLD(i, j, 2)] = rhoR_star * u_star * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoR_star * (u_star * u_star + v_R * v_R)) * u_star;
                        }
                    }
                    else
                    {
                        if (S_HR < 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_R * u_R;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_R * u_R * u_R + p_R;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_R * u_R * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_R(i, j, 3)] + p_R) * u_R;
                        }
                        else
                        {
                            if (S_TR > 0)
                            {
                                U_TMP[Index_F_OLD(i, j, 0)] = rhoR_star * u_star;
                                U_TMP[Index_F_OLD(i, j, 1)] = rhoR_star * u_star * u_star + p_star_LR;
                                U_TMP[Index_F_OLD(i, j, 2)] = rhoR_star * u_star * v_R;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoR_star * (u_star * u_star + v_R * v_R)) * u_star;
                            }
                            else
                            {
                                double rho_R_fan = rho_R * pow(2 / (gamma + 1) + (gamma - 1) * (rho_R) / (gamma + 1) / a_R, 2 / (gamma - 1));
                                double u_R_fan = 2 * (a_R + (gamma - 1) * u_R / 2) / (gamma + 1);
                                double p_R_fan = p_R * pow((2 / (gamma + 1) + (gamma - 1) * (rho_R) / (gamma + 1) / a_R), 2 * gamma / (gamma - 1));
                                U_TMP[Index_F_OLD(i, j, 0)] = rho_R_fan * u_R_fan;
                                U_TMP[Index_F_OLD(i, j, 1)] = rho_R_fan * u_R_fan * u_R_fan + p_R_fan;
                                U_TMP[Index_F_OLD(i, j, 2)] = rho_R_fan * u_R_fan * v_R;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_R_fan / (gamma - 1) + 0.5 * rho_R_fan * (u_R_fan * u_R_fan + v_R * v_R)) * u_R_fan;
                            }
                        }
                    }
                }
                if (v_star > 0)
                {
                    if (p_star_DU > p_D)
                    {
                        if (S_D > 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_D * u_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_D * v_D * v_D + p_D;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_D(i, j, 3)] + p_D) * v_D;
                        }
                        else
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rhoD_star * v_star;
                            U_TMP[Index_G_OLD(i, j, 1)] = rhoD_star * v_star * u_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rhoD_star * v_star * v_star + p_star_DU;
                            U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoD_star * (v_star * v_star + u_D * u_D)) * v_star;
                        }
                    }
                    else
                    {
                        if (S_HD > 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_D * u_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_D * v_D * v_D + p_D;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_D(i, j, 3)] + p_D) * v_D;
                        }
                        else
                        {
                            if (S_TD < 0)
                            {
                                U_TMP[Index_G_OLD(i, j, 0)] = rhoD_star * v_star;
                                U_TMP[Index_G_OLD(i, j, 1)] = rhoD_star * v_star * u_D;
                                U_TMP[Index_G_OLD(i, j, 2)] = rhoD_star * v_star * v_star + p_star_DU;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoD_star * (v_star * v_star + u_D * u_D)) * v_star;
                            }
                            else
                            {
                                double rho_D_fan = rho_D * pow(2 / (gamma + 1) + (gamma - 1) * (rho_D) / (gamma + 1) / a_D, 2 / (gamma - 1));
                                double v_D_fan = 2 * (a_D + (gamma - 1) * v_D / 2) / (gamma + 1);
                                double p_D_fan = p_D * pow((2 / (gamma + 1) + (gamma - 1) * (rho_D) / (gamma + 1) / a_D), 2 * gamma / (gamma - 1));
                                U_TMP[Index_G_OLD(i, j, 0)] = rho_D_fan * v_D_fan;
                                U_TMP[Index_G_OLD(i, j, 1)] = rho_D_fan * v_D_fan * u_D;
                                U_TMP[Index_G_OLD(i, j, 2)] = rho_D_fan * v_D_fan * v_D_fan + p_D_fan;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_D_fan / (gamma - 1) + 0.5 * rho_D_fan * (v_D_fan * v_D_fan + u_D * u_D)) * v_D_fan;
                            }
                        }
                    }
                }
                else
                {
                    if (p_star_DU > p_U)
                    {
                        if (S_U < 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_U * u_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_U * v_U * v_U + p_U;
                            U_TMP[Index(i, j, 3)] = (U_TMP[Index_U_U(i, j, 3)] + p_U) * v_U;
                        }
                        else
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rhoU_star * v_star;
                            U_TMP[Index_G_OLD(i, j, 1)] = rhoU_star * v_star * u_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rhoU_star * v_star * v_star + p_star_DU;
                            U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoU_star * (v_star * v_star + u_U * u_U)) * v_star;
                        }
                    }
                    else
                    {
                        if (S_HU < 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_U * u_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_U * v_U * v_U + p_U;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_U(i, j, 3)] + p_U) * v_U;
                        }
                        else
                        {
                            if (S_TU > 0)
                            {
                                U_TMP[Index_G_OLD(i, j, 0)] = rhoU_star * v_star;
                                U_TMP[Index_G_OLD(i, j, 1)] = rhoU_star * v_star * u_U;
                                U_TMP[Index_G_OLD(i, j, 2)] = rhoU_star * v_star * v_star + p_star_DU;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoU_star * (v_star * v_star + u_U * u_U)) * v_star;
                            }
                            else
                            {
                                double rho_U_fan = rho_U * pow(2 / (gamma + 1) + (gamma - 1) * (rho_U) / (gamma + 1) / a_U, 2 / (gamma - 1));
                                double v_U_fan = 2 * (a_U + (gamma - 1) * v_U / 2) / (gamma + 1);
                                double p_U_fan = p_U * pow((2 / (gamma + 1) + (gamma - 1) * (rho_U) / (gamma + 1) / a_U), 2 * gamma / (gamma - 1));
                                U_TMP[Index_G_OLD(i, j, 0)] = rho_U_fan * v_U_fan;
                                U_TMP[Index_G_OLD(i, j, 1)] = rho_U_fan * v_U_fan * u_U;
                                U_TMP[Index_G_OLD(i, j, 2)] = rho_U_fan * v_U_fan * v_U_fan + p_U_fan;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_U_fan / (gamma - 1) + 0.5 * rho_U_fan * (v_U_fan * v_U_fan + u_U * u_U)) * v_U_fan;
                            }
                        }
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
    long int lower = LOWER;
    long int upper = UPPER;
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size], \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
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
                double u_star = 0.5 * (u_L + u_R) + 0.5 * (p_L - p_R) / (rho_bar_LR * a_bar_LR);
                double rhoL_star = rho_L + (p_star_LR - p_L) / a_L / a_L;
                double rhoR_star = rho_R + (p_star_LR - p_R) / a_R / a_R;
                double S_L = u_L - a_L * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_L + (gamma - 1) / 2 / gamma);
                double S_R = u_R + a_R * sqrt((gamma + 1) * p_star_LR / 2 / gamma / p_R + (gamma - 1) / 2 / gamma);

                double S_HL = u_L - a_L;
                double S_TL = u_star - pow(a_L * (p_star_LR / p_R), ((gamma - 1) / 2 / gamma));

                double S_HR = u_R + a_R;
                double S_TR = u_star + pow(a_R * (p_star_LR / p_R), ((gamma - 1) / 2 / gamma));

                if (u_star > 0)
                {
                    if (p_star_LR > p_L)
                    {
                        if (S_L > 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_L * u_L;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_L * u_L * u_L + p_L;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_L * u_L * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_L(i, j, 3)] + p_L) * u_L;
                        }
                        else
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rhoL_star * u_star;
                            U_TMP[Index_F_OLD(i, j, 1)] = rhoL_star * u_star * u_star + p_star_LR;
                            U_TMP[Index_F_OLD(i, j, 2)] = rhoL_star * u_star * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoL_star * (u_star * u_star + v_L * v_L)) * u_star;
                        }
                    }
                    else
                    {
                        if (S_HL > 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_L * u_L;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_L * u_L * u_L + p_L;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_L * u_L * v_L;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_L(i, j, 3)] + p_L) * u_L;
                        }
                        else
                        {
                            if (S_TL < 0)
                            {
                                U_TMP[Index_F_OLD(i, j, 0)] = rhoL_star * u_star;
                                U_TMP[Index_F_OLD(i, j, 1)] = rhoL_star * u_star * u_star + p_star_LR;
                                U_TMP[Index_F_OLD(i, j, 2)] = rhoL_star * u_star * v_L;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoL_star * (u_star * u_star + v_L * v_L)) * u_star;
                            }
                            else
                            {
                                double rho_L_fan = rho_L * pow(2 / (gamma + 1) + (gamma - 1) * (rho_L) / (gamma + 1) / a_L, 2 / (gamma - 1));
                                double u_L_fan = 2 * (a_L + (gamma - 1) * u_L / 2) / (gamma + 1);
                                double p_L_fan = p_L * pow((2 / (gamma + 1) + (gamma - 1) * (rho_L) / (gamma + 1) / a_L), 2 * gamma / (gamma - 1));
                                U_TMP[Index_F_OLD(i, j, 0)] = rho_L_fan * u_L_fan;
                                U_TMP[Index_F_OLD(i, j, 1)] = rho_L_fan * u_L_fan * u_L_fan + p_L_fan;
                                U_TMP[Index_F_OLD(i, j, 2)] = rho_L_fan * u_L_fan * v_L;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_L_fan / (gamma - 1) + 0.5 * rho_L_fan * (u_L_fan * u_L_fan + v_L * v_L)) * u_L_fan;
                            }
                        }
                    }
                }
                else
                {
                    if (p_star_LR > p_R)
                    {
                        if (S_R < 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_R * u_R;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_R * u_R * u_R + p_R;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_R * u_R * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_R(i, j, 3)] + p_R) * u_R;
                        }
                        else
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rhoR_star * u_star;
                            U_TMP[Index_F_OLD(i, j, 1)] = rhoR_star * u_star * u_star + p_star_LR;
                            U_TMP[Index_F_OLD(i, j, 2)] = rhoR_star * u_star * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoR_star * (u_star * u_star + v_R * v_R)) * u_star;
                        }
                    }
                    else
                    {
                        if (S_HR < 0)
                        {
                            U_TMP[Index_F_OLD(i, j, 0)] = rho_R * u_R;
                            U_TMP[Index_F_OLD(i, j, 1)] = rho_R * u_R * u_R + p_R;
                            U_TMP[Index_F_OLD(i, j, 2)] = rho_R * u_R * v_R;
                            U_TMP[Index_F_OLD(i, j, 3)] = (U_TMP[Index_U_R(i, j, 3)] + p_R) * u_R;
                        }
                        else
                        {
                            if (S_TR > 0)
                            {
                                U_TMP[Index_F_OLD(i, j, 0)] = rhoR_star * u_star;
                                U_TMP[Index_F_OLD(i, j, 1)] = rhoR_star * u_star * u_star + p_star_LR;
                                U_TMP[Index_F_OLD(i, j, 2)] = rhoR_star * u_star * v_R;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_star_LR / (gamma - 1) + 0.5 * rhoR_star * (u_star * u_star + v_R * v_R)) * u_star;
                            }
                            else
                            {
                                double rho_R_fan = rho_R * pow(2 / (gamma + 1) + (gamma - 1) * (rho_R) / (gamma + 1) / a_R, 2 / (gamma - 1));
                                double u_R_fan = 2 * (a_R + (gamma - 1) * u_R / 2) / (gamma + 1);
                                double p_R_fan = p_R * pow((2 / (gamma + 1) + (gamma - 1) * (rho_R) / (gamma + 1) / a_R), 2 * gamma / (gamma - 1));
                                U_TMP[Index_F_OLD(i, j, 0)] = rho_R_fan * u_R_fan;
                                U_TMP[Index_F_OLD(i, j, 1)] = rho_R_fan * u_R_fan * u_R_fan + p_R_fan;
                                U_TMP[Index_F_OLD(i, j, 2)] = rho_R_fan * u_R_fan * v_R;
                                U_TMP[Index_F_OLD(i, j, 3)] = (gamma * p_R_fan / (gamma - 1) + 0.5 * rho_R_fan * (u_R_fan * u_R_fan + v_R * v_R)) * u_R_fan;
                            }
                        }
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
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size], \
                         U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
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
                double v_star = 0.5 * (v_D + v_U) + 0.5 * (p_D - p_U) / (rho_bar_DU * a_bar_DU);
                double rhoD_star = rho_D + (p_star_DU - p_D) / a_D / a_D;
                double rhoU_star = rho_U + (p_star_DU - p_U) / a_U / a_U;

                double S_D = v_D - a_D * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_D + (gamma - 1) / 2 / gamma);
                double S_U = v_U + a_U * sqrt((gamma + 1) * p_star_DU / 2 / gamma / p_U + (gamma - 1) / 2 / gamma);

                double S_HD = v_D - a_D;
                double S_TD = v_star - pow(a_D * (p_star_DU / p_U), ((gamma - 1) / 2 / gamma));

                double S_HU = v_U + a_U;
                double S_TU = v_star + pow(a_U * (p_star_DU / p_U), ((gamma - 1) / 2 / gamma));
                if (v_star > 0)
                {
                    if (p_star_DU > p_D)
                    {
                        if (S_D > 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_D * u_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_D * v_D * v_D + p_D;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_D(i, j, 3)] + p_D) * v_D;
                        }
                        else
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rhoD_star * v_star;
                            U_TMP[Index_G_OLD(i, j, 1)] = rhoD_star * v_star * u_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rhoD_star * v_star * v_star + p_star_DU;
                            U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoD_star * (v_star * v_star + u_D * u_D)) * v_star;
                        }
                    }
                    else
                    {
                        if (S_HD > 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_D * u_D * v_D;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_D * v_D * v_D + p_D;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_D(i, j, 3)] + p_D) * v_D;
                        }
                        else
                        {
                            if (S_TD < 0)
                            {
                                U_TMP[Index_G_OLD(i, j, 0)] = rhoD_star * v_star;
                                U_TMP[Index_G_OLD(i, j, 1)] = rhoD_star * v_star * u_D;
                                U_TMP[Index_G_OLD(i, j, 2)] = rhoD_star * v_star * v_star + p_star_DU;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoD_star * (v_star * v_star + u_D * u_D)) * v_star;
                            }
                            else
                            {
                                double rho_D_fan = rho_D * pow(2 / (gamma + 1) + (gamma - 1) * (rho_D) / (gamma + 1) / a_D, 2 / (gamma - 1));
                                double v_D_fan = 2 * (a_D + (gamma - 1) * v_D / 2) / (gamma + 1);
                                double p_D_fan = p_D * pow((2 / (gamma + 1) + (gamma - 1) * (rho_D) / (gamma + 1) / a_D), 2 * gamma / (gamma - 1));
                                U_TMP[Index_G_OLD(i, j, 0)] = rho_D_fan * v_D_fan;
                                U_TMP[Index_G_OLD(i, j, 1)] = rho_D_fan * v_D_fan * u_D;
                                U_TMP[Index_G_OLD(i, j, 2)] = rho_D_fan * v_D_fan * v_D_fan + p_D_fan;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_D_fan / (gamma - 1) + 0.5 * rho_D_fan * (v_D_fan * v_D_fan + u_D * u_D)) * v_D_fan;
                            }
                        }
                    }
                }
                else
                {
                    if (p_star_DU > p_U)
                    {
                        if (S_U < 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_U * u_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_U * v_U * v_U + p_U;
                            U_TMP[Index(i, j, 3)] = (U_TMP[Index_U_U(i, j, 3)] + p_U) * v_U;
                        }
                        else
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rhoU_star * v_star;
                            U_TMP[Index_G_OLD(i, j, 1)] = rhoU_star * v_star * u_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rhoU_star * v_star * v_star + p_star_DU;
                            U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoU_star * (v_star * v_star + u_U * u_U)) * v_star;
                        }
                    }
                    else
                    {
                        if (S_HU < 0)
                        {
                            U_TMP[Index_G_OLD(i, j, 0)] = rho_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 1)] = rho_U * u_U * v_U;
                            U_TMP[Index_G_OLD(i, j, 2)] = rho_U * v_U * v_U + p_U;
                            U_TMP[Index_G_OLD(i, j, 3)] = (U_TMP[Index_U_U(i, j, 3)] + p_U) * v_U;
                        }
                        else
                        {
                            if (S_TU > 0)
                            {
                                U_TMP[Index_G_OLD(i, j, 0)] = rhoU_star * v_star;
                                U_TMP[Index_G_OLD(i, j, 1)] = rhoU_star * v_star * u_U;
                                U_TMP[Index_G_OLD(i, j, 2)] = rhoU_star * v_star * v_star + p_star_DU;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_star_DU / (gamma - 1) + 0.5 * rhoU_star * (v_star * v_star + u_U * u_U)) * v_star;
                            }
                            else
                            {
                                double rho_U_fan = rho_U * pow(2 / (gamma + 1) + (gamma - 1) * (rho_U) / (gamma + 1) / a_U, 2 / (gamma - 1));
                                double v_U_fan = 2 * (a_U + (gamma - 1) * v_U / 2) / (gamma + 1);
                                double p_U_fan = p_U * pow((2 / (gamma + 1) + (gamma - 1) * (rho_U) / (gamma + 1) / a_U), 2 * gamma / (gamma - 1));
                                U_TMP[Index_G_OLD(i, j, 0)] = rho_U_fan * v_U_fan;
                                U_TMP[Index_G_OLD(i, j, 1)] = rho_U_fan * v_U_fan * u_U;
                                U_TMP[Index_G_OLD(i, j, 2)] = rho_U_fan * v_U_fan * v_U_fan + p_U_fan;
                                U_TMP[Index_G_OLD(i, j, 3)] = (gamma * p_U_fan / (gamma - 1) + 0.5 * rho_U_fan * (v_U_fan * v_U_fan + u_U * u_U)) * v_U_fan;
                            }
                        }
                    }
                }
            }
        }
    }
}