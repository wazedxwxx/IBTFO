#include "Conserve2Flux.H"
#include <math.h>
#include "EQDefine.H"
#include "CoordDefine.H"
#include <algorithm>
using namespace std;
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

                double e_L = p_L / (gamma - 1) / rho_L;
                double e_R = p_R / (gamma - 1) / rho_R;
                double e_U = p_U / (gamma - 1) / rho_U;
                double e_D = p_D / (gamma - 1) / rho_D;

                double h_L = e_L + u_L * u_L / 2 + v_L * v_L / 2 + p_L / rho_L;
                double h_R = e_R + u_R * u_R / 2 + v_R * v_R / 2 + p_R / rho_R;
                double h_U = e_U + u_U * u_U / 2 + v_U * v_U / 2 + p_U / rho_U;
                double h_D = e_D + u_D * u_D / 2 + v_D * v_D / 2 + p_D / rho_D;

                double delta = 1e-6;

                double lamta1_L_p = 0.5 * (u_L + sqrt(u_L * u_L + delta * delta));                         // lamta_L+
                double lamta1_L_n = 0.5 * (u_L - sqrt(u_L * u_L + delta * delta));                         // lamta_L-
                double lamta2_L_p = 0.5 * (u_L + sqrt(u_L * u_L + delta * delta));                         // lamta_L+
                double lamta2_L_n = 0.5 * (u_L - sqrt(u_L * u_L + delta * delta));                         // lamta_L-
                double lamta3_L_p = 0.5 * ((u_L - a_L) + sqrt((u_L - a_L) * (u_L - a_L) + delta * delta)); // lamta_L+
                double lamta3_L_n = 0.5 * ((u_L - a_L) - sqrt((u_L - a_L) * (u_L - a_L) + delta * delta)); // lamta_L-
                double lamta4_L_p = 0.5 * ((u_L + a_L) + sqrt((u_L + a_L) * (u_L + a_L) + delta * delta)); // lamta_L+
                double lamta4_L_n = 0.5 * ((u_L + a_L) - sqrt((u_L + a_L) * (u_L + a_L) + delta * delta)); // lamta_L-

                double lamta1_R_p = 0.5 * (u_R + sqrt(u_R * u_R + delta * delta));                         // lamta_R+
                double lamta1_R_n = 0.5 * (u_R - sqrt(u_R * u_R + delta * delta));                         // lamta_R-
                double lamta2_R_p = 0.5 * (u_R + sqrt(u_R * u_R + delta * delta));                         // lamta_R+
                double lamta2_R_n = 0.5 * (u_R - sqrt(u_R * u_R + delta * delta));                         // lamta_R-
                double lamta3_R_p = 0.5 * ((u_R - a_R) + sqrt((u_R - a_R) * (u_R - a_R) + delta * delta)); // lamta_R+
                double lamta3_R_n = 0.5 * ((u_R - a_R) - sqrt((u_R - a_R) * (u_R - a_R) + delta * delta)); // lamta_R-
                double lamta4_R_p = 0.5 * ((u_R + a_R) + sqrt((u_R + a_R) * (u_R + a_R) + delta * delta)); // lamta_R+
                double lamta4_R_n = 0.5 * ((u_R + a_R) - sqrt((u_R + a_R) * (u_R + a_R) + delta * delta)); // lamta_R-

                double mu1_D_p = 0.5 * (v_D + sqrt(v_D * v_D + delta * delta));                         // mu_D+
                double mu1_D_n = 0.5 * (v_D - sqrt(v_D * v_D + delta * delta));                         // mu_D-
                double mu2_D_p = 0.5 * (v_D + sqrt(v_D * v_D + delta * delta));                         // mu_D+
                double mu2_D_n = 0.5 * (v_D - sqrt(v_D * v_D + delta * delta));                         // mu_D-
                double mu3_D_p = 0.5 * ((v_D - a_D) + sqrt((v_D - a_D) * (v_D - a_D) + delta * delta)); // mu_D+
                double mu3_D_n = 0.5 * ((v_D - a_D) - sqrt((v_D - a_D) * (v_D - a_D) + delta * delta)); // mu_D-
                double mu4_D_p = 0.5 * ((v_D + a_D) + sqrt((v_D + a_D) * (v_D + a_D) + delta * delta)); // mu_D+
                double mu4_D_n = 0.5 * ((v_D + a_D) - sqrt((v_D + a_D) * (v_D + a_D) + delta * delta)); // mu_D-

                double mu1_U_p = 0.5 * (v_U + sqrt(v_U * v_U + delta * delta));                         // mu_U+
                double mu1_U_n = 0.5 * (v_U - sqrt(v_U * v_U + delta * delta));                         // mu_U-
                double mu2_U_p = 0.5 * (v_U + sqrt(v_U * v_U + delta * delta));                         // mu_U+
                double mu2_U_n = 0.5 * (v_U - sqrt(v_U * v_U + delta * delta));                         // mu_U-
                double mu3_U_p = 0.5 * ((v_U - a_U) + sqrt((v_U - a_U) * (v_U - a_U) + delta * delta)); // mu_U+
                double mu3_U_n = 0.5 * ((v_U - a_U) - sqrt((v_U - a_U) * (v_U - a_U) + delta * delta)); // mu_U-
                double mu4_U_p = 0.5 * ((v_U + a_U) + sqrt((v_U + a_U) * (v_U + a_U) + delta * delta)); // mu_U+
                double mu4_U_n = 0.5 * ((v_U + a_U) - sqrt((v_U + a_U) * (v_U + a_U) + delta * delta)); // mu_U-

                double f1_L_p = rho_L * (2 * (gamma - 1) * lamta1_L_p + lamta3_L_p + lamta4_L_p) / 2 / gamma;
                double f1_L_n = rho_L * (2 * (gamma - 1) * lamta1_L_n + lamta3_L_n + lamta4_L_n) / 2 / gamma;
                double f2_L_p = rho_L * (2 * (gamma - 1) * u_L * lamta1_L_p + (u_L - a_L) * lamta3_L_p + (u_L + a_L) * lamta4_L_p) / 2 / gamma;
                double f2_L_n = rho_L * (2 * (gamma - 1) * u_L * lamta1_L_n + (u_L - a_L) * lamta3_L_n + (u_L + a_L) * lamta4_L_n) / 2 / gamma;
                double f3_L_p = rho_L * (2 * (gamma - 1) * v_L * lamta1_L_p + lamta3_L_p * v_L + lamta4_L_p * v_L) / 2 / gamma;
                double f3_L_n = rho_L * (2 * (gamma - 1) * v_L * lamta1_L_n + lamta3_L_n * v_L + lamta4_L_n * v_L) / 2 / gamma;
                double f4_L_p = rho_L * ((gamma - 1) * (u_L * u_L + v_L * v_L) * lamta1_L_p + (h_L - a_L * u_L) * lamta3_L_p + (h_L + a_L * u_L) * lamta4_L_p) / 2 / gamma;
                double f4_L_n = rho_L * ((gamma - 1) * (u_L * u_L + v_L * v_L) * lamta1_L_n + (h_L - a_L * u_L) * lamta3_L_n + (h_L + a_L * u_L) * lamta4_L_n) / 2 / gamma;

                double f1_R_p = rho_R * (2 * (gamma - 1) * lamta1_R_p + lamta3_R_p + lamta4_R_p) / 2 / gamma;
                double f1_R_n = rho_R * (2 * (gamma - 1) * lamta1_R_n + lamta3_R_n + lamta4_R_n) / 2 / gamma;
                double f2_R_p = rho_R * (2 * (gamma - 1) * u_R * lamta1_R_p + (u_R - a_R) * lamta3_R_p + (u_R + a_R) * lamta4_R_p) / 2 / gamma;
                double f2_R_n = rho_R * (2 * (gamma - 1) * u_R * lamta1_R_n + (u_R - a_R) * lamta3_R_n + (u_R + a_R) * lamta4_R_n) / 2 / gamma;
                double f3_R_p = rho_R * (2 * (gamma - 1) * v_R * lamta1_R_p + lamta3_R_p * v_R + lamta4_R_p * v_R) / 2 / gamma;
                double f3_R_n = rho_R * (2 * (gamma - 1) * v_R * lamta1_R_n + lamta3_R_n * v_R + lamta4_R_n * v_R) / 2 / gamma;
                double f4_R_p = rho_R * ((gamma - 1) * (u_R * u_R + v_R * v_R) * lamta1_R_p + (h_R - a_R * u_R) * lamta3_R_p + (h_R + a_R * u_R) * lamta4_R_p) / 2 / gamma;
                double f4_R_n = rho_R * ((gamma - 1) * (u_R * u_R + v_R * v_R) * lamta1_R_n + (h_R - a_R * u_R) * lamta3_R_n + (h_R + a_R * u_R) * lamta4_R_n) / 2 / gamma;

                double g1_D_p = rho_D * (2 * (gamma - 1) * mu1_D_p + mu3_D_p + mu4_D_p) / 2 / gamma;
                double g1_D_n = rho_D * (2 * (gamma - 1) * mu1_D_n + mu3_D_n + mu4_D_n) / 2 / gamma;
                double g2_D_p = rho_D * (2 * (gamma - 1) * u_D * mu1_D_p + mu3_D_p * u_D + mu4_D_p * u_D) / 2 / gamma;
                double g2_D_n = rho_D * (2 * (gamma - 1) * u_D * mu1_D_n + mu3_D_n * u_D + mu4_D_n * u_D) / 2 / gamma;
                double g3_D_p = rho_D * (2 * (gamma - 1) * v_D * mu1_D_p + (v_D - a_D) * mu3_D_p + (v_D + a_D) * mu4_D_p) / 2 / gamma;
                double g3_D_n = rho_D * (2 * (gamma - 1) * v_D * mu1_D_n + (v_D - a_D) * mu3_D_n + (v_D + a_D) * mu4_D_n) / 2 / gamma;
                double g4_D_p = rho_D * ((gamma - 1) * (u_D * u_D + v_D * v_D) * mu1_D_p + (h_D - a_D * v_D) * mu3_D_p + (h_D + a_D * v_D) * mu4_D_p) / 2 / gamma;
                double g4_D_n = rho_D * ((gamma - 1) * (u_D * u_D + v_D * v_D) * mu1_D_n + (h_D - a_D * v_D) * mu3_D_n + (h_D + a_D * v_D) * mu4_D_n) / 2 / gamma;

                double g1_U_p = rho_U * (2 * (gamma - 1) * mu1_U_p + mu3_U_p + mu4_U_p) / 2 / gamma;
                double g1_U_n = rho_U * (2 * (gamma - 1) * mu1_U_n + mu3_U_n + mu4_U_n) / 2 / gamma;
                double g2_U_p = rho_U * (2 * (gamma - 1) * u_U * mu1_U_p + mu3_U_p * u_U + mu4_U_p * u_U) / 2 / gamma;
                double g2_U_n = rho_U * (2 * (gamma - 1) * u_U * mu1_U_n + mu3_U_n * u_U + mu4_U_n * u_U) / 2 / gamma;
                double g3_U_p = rho_U * (2 * (gamma - 1) * v_U * mu1_U_p + (v_U - a_U) * mu3_U_p + (v_U + a_U) * mu4_U_p) / 2 / gamma;
                double g3_U_n = rho_U * (2 * (gamma - 1) * v_U * mu1_U_n + (v_U - a_U) * mu3_U_n + (v_U + a_U) * mu4_U_n) / 2 / gamma;
                double g4_U_p = rho_U * ((gamma - 1) * (u_U * u_U + v_U * v_U) * mu1_U_p + (h_U - a_U * v_U) * mu3_U_p + (h_U + a_U * v_U) * mu4_U_p) / 2 / gamma;
                double g4_U_n = rho_U * ((gamma - 1) * (u_U * u_U + v_U * v_U) * mu1_U_n + (h_U - a_U * v_U) * mu3_U_n + (h_U + a_U * v_U) * mu4_U_n) / 2 / gamma;

                U_TMP[Index_F_OLD(i, j, 0)] = f1_L_p + f1_R_n;
                U_TMP[Index_F_OLD(i, j, 1)] = f2_L_p + f2_R_n;
                U_TMP[Index_F_OLD(i, j, 2)] = f3_L_p + f3_R_n;
                U_TMP[Index_F_OLD(i, j, 3)] = f4_L_p + f4_R_n;

                U_TMP[Index_G_OLD(i, j, 0)] = g1_D_p + g1_U_n;
                U_TMP[Index_G_OLD(i, j, 1)] = g2_D_p + g2_U_n;
                U_TMP[Index_G_OLD(i, j, 2)] = g3_D_p + g3_U_n;
                U_TMP[Index_G_OLD(i, j, 3)] = g4_D_p + g4_U_n;
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
#pragma acc data present(U_TMP [(M)*lower * num_tmp_size:(M) * (upper - lower) * num_tmp_size],   \
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

                double e_L = p_L / (gamma - 1) / rho_L;
                double e_R = p_R / (gamma - 1) / rho_R;

                double h_L = e_L + u_L * u_L / 2 + v_L * v_L / 2 + p_L / rho_L;
                double h_R = e_R + u_R * u_R / 2 + v_R * v_R / 2 + p_R / rho_R;

                double delta = 1e-6;

                double lamta1_L_p = 0.5 * (u_L + sqrt(u_L * u_L + delta * delta));                         // lamta_L+
                double lamta1_L_n = 0.5 * (u_L - sqrt(u_L * u_L + delta * delta));                         // lamta_L-
                double lamta2_L_p = 0.5 * (u_L + sqrt(u_L * u_L + delta * delta));                         // lamta_L+
                double lamta2_L_n = 0.5 * (u_L - sqrt(u_L * u_L + delta * delta));                         // lamta_L-
                double lamta3_L_p = 0.5 * ((u_L - a_L) + sqrt((u_L - a_L) * (u_L - a_L) + delta * delta)); // lamta_L+
                double lamta3_L_n = 0.5 * ((u_L - a_L) - sqrt((u_L - a_L) * (u_L - a_L) + delta * delta)); // lamta_L-
                double lamta4_L_p = 0.5 * ((u_L + a_L) + sqrt((u_L + a_L) * (u_L + a_L) + delta * delta)); // lamta_L+
                double lamta4_L_n = 0.5 * ((u_L + a_L) - sqrt((u_L + a_L) * (u_L + a_L) + delta * delta)); // lamta_L-

                double lamta1_R_p = 0.5 * (u_R + sqrt(u_R * u_R + delta * delta));                         // lamta_R+
                double lamta1_R_n = 0.5 * (u_R - sqrt(u_R * u_R + delta * delta));                         // lamta_R-
                double lamta2_R_p = 0.5 * (u_R + sqrt(u_R * u_R + delta * delta));                         // lamta_R+
                double lamta2_R_n = 0.5 * (u_R - sqrt(u_R * u_R + delta * delta));                         // lamta_R-
                double lamta3_R_p = 0.5 * ((u_R - a_R) + sqrt((u_R - a_R) * (u_R - a_R) + delta * delta)); // lamta_R+
                double lamta3_R_n = 0.5 * ((u_R - a_R) - sqrt((u_R - a_R) * (u_R - a_R) + delta * delta)); // lamta_R-
                double lamta4_R_p = 0.5 * ((u_R + a_R) + sqrt((u_R + a_R) * (u_R + a_R) + delta * delta)); // lamta_R+
                double lamta4_R_n = 0.5 * ((u_R + a_R) - sqrt((u_R + a_R) * (u_R + a_R) + delta * delta)); // lamta_R-

                double f1_L_p = rho_L * (2 * (gamma - 1) * lamta1_L_p + lamta3_L_p + lamta4_L_p) / 2 / gamma;
                double f1_L_n = rho_L * (2 * (gamma - 1) * lamta1_L_n + lamta3_L_n + lamta4_L_n) / 2 / gamma;
                double f2_L_p = rho_L * (2 * (gamma - 1) * u_L * lamta1_L_p + (u_L - a_L) * lamta3_L_p + (u_L + a_L) * lamta4_L_p) / 2 / gamma;
                double f2_L_n = rho_L * (2 * (gamma - 1) * u_L * lamta1_L_n + (u_L - a_L) * lamta3_L_n + (u_L + a_L) * lamta4_L_n) / 2 / gamma;
                double f3_L_p = rho_L * (2 * (gamma - 1) * v_L * lamta1_L_p + lamta3_L_p * v_L + lamta4_L_p * v_L) / 2 / gamma;
                double f3_L_n = rho_L * (2 * (gamma - 1) * v_L * lamta1_L_n + lamta3_L_n * v_L + lamta4_L_n * v_L) / 2 / gamma;
                double f4_L_p = rho_L * ((gamma - 1) * (u_L * u_L + v_L * v_L) * lamta1_L_p + (h_L - a_L * u_L) * lamta3_L_p + (h_L + a_L * u_L) * lamta4_L_p) / 2 / gamma;
                double f4_L_n = rho_L * ((gamma - 1) * (u_L * u_L + v_L * v_L) * lamta1_L_n + (h_L - a_L * u_L) * lamta3_L_n + (h_L + a_L * u_L) * lamta4_L_n) / 2 / gamma;

                double f1_R_p = rho_R * (2 * (gamma - 1) * lamta1_R_p + lamta3_R_p + lamta4_R_p) / 2 / gamma;
                double f1_R_n = rho_R * (2 * (gamma - 1) * lamta1_R_n + lamta3_R_n + lamta4_R_n) / 2 / gamma;
                double f2_R_p = rho_R * (2 * (gamma - 1) * u_R * lamta1_R_p + (u_R - a_R) * lamta3_R_p + (u_R + a_R) * lamta4_R_p) / 2 / gamma;
                double f2_R_n = rho_R * (2 * (gamma - 1) * u_R * lamta1_R_n + (u_R - a_R) * lamta3_R_n + (u_R + a_R) * lamta4_R_n) / 2 / gamma;
                double f3_R_p = rho_R * (2 * (gamma - 1) * v_R * lamta1_R_p + lamta3_R_p * v_R + lamta4_R_p * v_R) / 2 / gamma;
                double f3_R_n = rho_R * (2 * (gamma - 1) * v_R * lamta1_R_n + lamta3_R_n * v_R + lamta4_R_n * v_R) / 2 / gamma;
                double f4_R_p = rho_R * ((gamma - 1) * (u_R * u_R + v_R * v_R) * lamta1_R_p + (h_R - a_R * u_R) * lamta3_R_p + (h_R + a_R * u_R) * lamta4_R_p) / 2 / gamma;
                double f4_R_n = rho_R * ((gamma - 1) * (u_R * u_R + v_R * v_R) * lamta1_R_n + (h_R - a_R * u_R) * lamta3_R_n + (h_R + a_R * u_R) * lamta4_R_n) / 2 / gamma;

                U_TMP[Index_F_OLD(i, j, 0)] = f1_L_p + f1_R_n;
                U_TMP[Index_F_OLD(i, j, 1)] = f2_L_p + f2_R_n;
                U_TMP[Index_F_OLD(i, j, 2)] = f3_L_p + f3_R_n;
                U_TMP[Index_F_OLD(i, j, 3)] = f4_L_p + f4_R_n;
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

                double e_U = p_U / (gamma - 1) / rho_U;
                double e_D = p_D / (gamma - 1) / rho_D;

                double h_U = e_U + u_U * u_U / 2 + v_U * v_U / 2 + p_U / rho_U;
                double h_D = e_D + u_D * u_D / 2 + v_D * v_D / 2 + p_D / rho_D;

                double delta = 1e-6;

                double mu1_D_p = 0.5 * (v_D + sqrt(v_D * v_D + delta * delta));                         // mu_D+
                double mu1_D_n = 0.5 * (v_D - sqrt(v_D * v_D + delta * delta));                         // mu_D-
                double mu2_D_p = 0.5 * (v_D + sqrt(v_D * v_D + delta * delta));                         // mu_D+
                double mu2_D_n = 0.5 * (v_D - sqrt(v_D * v_D + delta * delta));                         // mu_D-
                double mu3_D_p = 0.5 * ((v_D - a_D) + sqrt((v_D - a_D) * (v_D - a_D) + delta * delta)); // mu_D+
                double mu3_D_n = 0.5 * ((v_D - a_D) - sqrt((v_D - a_D) * (v_D - a_D) + delta * delta)); // mu_D-
                double mu4_D_p = 0.5 * ((v_D + a_D) + sqrt((v_D + a_D) * (v_D + a_D) + delta * delta)); // mu_D+
                double mu4_D_n = 0.5 * ((v_D + a_D) - sqrt((v_D + a_D) * (v_D + a_D) + delta * delta)); // mu_D-

                double mu1_U_p = 0.5 * (v_U + sqrt(v_U * v_U + delta * delta));                         // mu_U+
                double mu1_U_n = 0.5 * (v_U - sqrt(v_U * v_U + delta * delta));                         // mu_U-
                double mu2_U_p = 0.5 * (v_U + sqrt(v_U * v_U + delta * delta));                         // mu_U+
                double mu2_U_n = 0.5 * (v_U - sqrt(v_U * v_U + delta * delta));                         // mu_U-
                double mu3_U_p = 0.5 * ((v_U - a_U) + sqrt((v_U - a_U) * (v_U - a_U) + delta * delta)); // mu_U+
                double mu3_U_n = 0.5 * ((v_U - a_U) - sqrt((v_U - a_U) * (v_U - a_U) + delta * delta)); // mu_U-
                double mu4_U_p = 0.5 * ((v_U + a_U) + sqrt((v_U + a_U) * (v_U + a_U) + delta * delta)); // mu_U+
                double mu4_U_n = 0.5 * ((v_U + a_U) - sqrt((v_U + a_U) * (v_U + a_U) + delta * delta)); // mu_U-


                double g1_D_p = rho_D * (2 * (gamma - 1) * mu1_D_p + mu3_D_p + mu4_D_p) / 2 / gamma;
                double g1_D_n = rho_D * (2 * (gamma - 1) * mu1_D_n + mu3_D_n + mu4_D_n) / 2 / gamma;
                double g2_D_p = rho_D * (2 * (gamma - 1) * u_D * mu1_D_p + mu3_D_p * u_D + mu4_D_p * u_D) / 2 / gamma;
                double g2_D_n = rho_D * (2 * (gamma - 1) * u_D * mu1_D_n + mu3_D_n * u_D + mu4_D_n * u_D) / 2 / gamma;
                double g3_D_p = rho_D * (2 * (gamma - 1) * v_D * mu1_D_p + (v_D - a_D) * mu3_D_p + (v_D + a_D) * mu4_D_p) / 2 / gamma;
                double g3_D_n = rho_D * (2 * (gamma - 1) * v_D * mu1_D_n + (v_D - a_D) * mu3_D_n + (v_D + a_D) * mu4_D_n) / 2 / gamma;
                double g4_D_p = rho_D * ((gamma - 1) * (u_D * u_D + v_D * v_D) * mu1_D_p + (h_D - a_D * v_D) * mu3_D_p + (h_D + a_D * v_D) * mu4_D_p) / 2 / gamma;
                double g4_D_n = rho_D * ((gamma - 1) * (u_D * u_D + v_D * v_D) * mu1_D_n + (h_D - a_D * v_D) * mu3_D_n + (h_D + a_D * v_D) * mu4_D_n) / 2 / gamma;

                double g1_U_p = rho_U * (2 * (gamma - 1) * mu1_U_p + mu3_U_p + mu4_U_p) / 2 / gamma;
                double g1_U_n = rho_U * (2 * (gamma - 1) * mu1_U_n + mu3_U_n + mu4_U_n) / 2 / gamma;
                double g2_U_p = rho_U * (2 * (gamma - 1) * u_U * mu1_U_p + mu3_U_p * u_U + mu4_U_p * u_U) / 2 / gamma;
                double g2_U_n = rho_U * (2 * (gamma - 1) * u_U * mu1_U_n + mu3_U_n * u_U + mu4_U_n * u_U) / 2 / gamma;
                double g3_U_p = rho_U * (2 * (gamma - 1) * v_U * mu1_U_p + (v_U - a_U) * mu3_U_p + (v_U + a_U) * mu4_U_p) / 2 / gamma;
                double g3_U_n = rho_U * (2 * (gamma - 1) * v_U * mu1_U_n + (v_U - a_U) * mu3_U_n + (v_U + a_U) * mu4_U_n) / 2 / gamma;
                double g4_U_p = rho_U * ((gamma - 1) * (u_U * u_U + v_U * v_U) * mu1_U_p + (h_U - a_U * v_U) * mu3_U_p + (h_U + a_U * v_U) * mu4_U_p) / 2 / gamma;
                double g4_U_n = rho_U * ((gamma - 1) * (u_U * u_U + v_U * v_U) * mu1_U_n + (h_U - a_U * v_U) * mu3_U_n + (h_U + a_U * v_U) * mu4_U_n) / 2 / gamma;

                U_TMP[Index_G_OLD(i, j, 0)] = g1_D_p + g1_U_n;
                U_TMP[Index_G_OLD(i, j, 1)] = g2_D_p + g2_U_n;
                U_TMP[Index_G_OLD(i, j, 2)] = g3_D_p + g3_U_n;
                U_TMP[Index_G_OLD(i, j, 3)] = g4_D_p + g4_U_n;
            }
        }
    }
}