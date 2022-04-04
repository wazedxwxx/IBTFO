#include "Boundary.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include <iostream>
#include <math.h>
void Boundary(const int N_x,
              const int N_y,
              const int num_ghost_cell,
              const double gamma,
              double *U_OLD,
              double *U_NEW,
              double *XYCOORD,
              double *GFM_Index)
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
                if (XYCOORD[Index_Coord(i, j, 5 )] == 0)
                    U_OLD[Index(i, j, k)] = U_NEW[Index(i, j, k)];
            }
        }
    }

#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
        {
#pragma acc loop
            for (int k = 0; k < num_eq; k++)
            {
                U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)];
            }
        }
    }

    // Ghost-cell
    //            i,j+1     i+1,j+1             a2=1/d2      a3=1/d3
    //                                 =>
    //            i,j        i+1,j              a1=1/d1      a4=1/d4
    //
    //

    /*
    Index Define
    0:  Image point x coord
    1:  Image point y coord
    2:  Image point neighbor point x index
    3:  Image point neighbor point y index
    4:  Image point inverse distance coffcient
    5:  Image point inverse distance coffcient
    6:  Image point inverse distance coffcient
    7:  Image point inverse distance coffcient
    8:  Extra point x coord
    9:  Extra point y coord
    10: Extra point neighbor point x index
    11: Extra point neighbor point y index
    12: Extra point inverse distance coffcient
    13: Extra point inverse distance coffcient
    14: Extra point inverse distance coffcient
    15: Extra point inverse distance coffcient

    */
    const double dx = XYCOORD[Index_Coord(1, 0, 0 )] - XYCOORD[Index_Coord(0, 0, 0 )];
    const double dy = XYCOORD[Index_Coord(0, 1, 0 )] - XYCOORD[Index_Coord(0, 0, 0 )];

    const double delta = sqrt(dx * dx + dy * dy) / 2;
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {

            if (XYCOORD[Index_Coord(i, j, 5 )] > 1)
            {
                double R_E = 1.5 * delta;
                double R_G = std::abs(XYCOORD[Index_Coord(i, j, 2 )]);

                double n_x = XYCOORD[Index_Coord(i, j, 3 )];
                double n_y = XYCOORD[Index_Coord(i, j, 4 )];

                double a1_image = GFM_Index[Index_sch(i, j, 4 )];
                double a2_image = GFM_Index[Index_sch(i, j, 5 )];
                double a3_image = GFM_Index[Index_sch(i, j, 6 )];
                double a4_image = GFM_Index[Index_sch(i, j, 7 )];

                int IDX = GFM_Index[Index_sch(i, j, 2 )];
                int IDY = GFM_Index[Index_sch(i, j, 3 )];
                double rho1_image = U_NEW[Index(IDX, IDY, 0)];
                double u1_image = U_NEW[Index(IDX, IDY, 1)] / rho1_image;
                double v1_image = U_NEW[Index(IDX, IDY, 2)] / rho1_image;
                double p1_image = (gamma - 1) * (U_NEW[Index(IDX, IDY, 3)] - 0.5 * rho1_image * (u1_image * u1_image + v1_image * v1_image));
                double rho2_image = U_NEW[Index(IDX, IDY + 1, 0)];
                double u2_image = U_NEW[Index(IDX, IDY + 1, 1)] / rho2_image;
                double v2_image = U_NEW[Index(IDX, IDY + 1, 2)] / rho2_image;
                double p2_image = (gamma - 1) * (U_NEW[Index(IDX, IDY + 1, 3)] - 0.5 * rho2_image * (u2_image * u2_image + v2_image * v2_image));
                double rho3_image = U_NEW[Index(IDX + 1, IDY + 1, 0)];
                double u3_image = U_NEW[Index(IDX + 1, IDY + 1, 1)] / rho3_image;
                double v3_image = U_NEW[Index(IDX + 1, IDY + 1, 2)] / rho3_image;
                double p3_image = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY + 1, 3)] - 0.5 * rho3_image * (u3_image * u3_image + v3_image * v3_image));
                double rho4_image = U_NEW[Index(IDX + 1, IDY, 0)];
                double u4_image = U_NEW[Index(IDX + 1, IDY, 1)] / rho4_image;
                double v4_image = U_NEW[Index(IDX + 1, IDY, 2)] / rho4_image;
                double p4_image = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY, 3)] - 0.5 * rho4_image * (u4_image * u4_image + v4_image * v4_image));

                double rho_image = (a1_image * rho1_image + a2_image * rho2_image + a3_image * rho3_image + a4_image * rho4_image) / (a1_image + a2_image + a3_image + a4_image);
                double u_image = (a1_image * u1_image + a2_image * u2_image + a3_image * u3_image + a4_image * u4_image) / (a1_image + a2_image + a3_image + a4_image);
                double v_image = (a1_image * v1_image + a2_image * v2_image + a3_image * v3_image + a4_image * v4_image) / (a1_image + a2_image + a3_image + a4_image);
                double p_image = (a1_image * p1_image + a2_image * p2_image + a3_image * p3_image + a4_image * p4_image) / (a1_image + a2_image + a3_image + a4_image);

                double u_n_image = u_image * n_x + v_image * n_y;
                double u_t_image = u_image * n_y - v_image * n_x;

                double a1_extra = GFM_Index[Index_sch(i, j, 12 )];
                double a2_extra = GFM_Index[Index_sch(i, j, 13 )];
                double a3_extra = GFM_Index[Index_sch(i, j, 14 )];
                double a4_extra = GFM_Index[Index_sch(i, j, 15 )];

                IDX = GFM_Index[Index_sch(i, j, 10 )];
                IDY = GFM_Index[Index_sch(i, j, 11 )];
                double rho1_extra = U_NEW[Index(IDX, IDY, 0)];
                double u1_extra = U_NEW[Index(IDX, IDY, 1)] / rho1_extra;
                double v1_extra = U_NEW[Index(IDX, IDY, 2)] / rho1_extra;
                double p1_extra = (gamma - 1) * (U_NEW[Index(IDX, IDY, 3)] - 0.5 * rho1_extra * (u1_extra * u1_extra + v1_extra * v1_extra));
                double rho2_extra = U_NEW[Index(IDX, IDY + 1, 0)];
                double u2_extra = U_NEW[Index(IDX, IDY + 1, 1)] / rho2_extra;
                double v2_extra = U_NEW[Index(IDX, IDY + 1, 2)] / rho2_extra;
                double p2_extra = (gamma - 1) * (U_NEW[Index(IDX, IDY + 1, 3)] - 0.5 * rho2_extra * (u2_extra * u2_extra + v2_extra * v2_extra));
                double rho3_extra = U_NEW[Index(IDX + 1, IDY + 1, 0)];
                double u3_extra = U_NEW[Index(IDX + 1, IDY + 1, 1)] / rho3_extra;
                double v3_extra = U_NEW[Index(IDX + 1, IDY + 1, 2)] / rho3_extra;
                double p3_extra = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY + 1, 3)] - 0.5 * rho3_extra * (u3_extra * u3_extra + v3_extra * v3_extra));
                double rho4_extra = U_NEW[Index(IDX + 1, IDY, 0)];
                double u4_extra = U_NEW[Index(IDX + 1, IDY, 1)] / rho4_extra;
                double v4_extra = U_NEW[Index(IDX + 1, IDY, 2)] / rho4_extra;
                double p4_extra = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY, 3)] - 0.5 * rho4_extra * (u4_extra * u4_extra + v4_extra * v4_extra));

                double rho_extra = (a1_extra * rho1_extra + a2_extra * rho2_extra + a3_extra * rho3_extra + a4_extra * rho4_extra) / (a1_extra + a2_extra + a3_extra + a4_extra);
                double u_extra = (a1_extra * u1_extra + a2_extra * u2_extra + a3_extra * u3_extra + a4_extra * u4_extra) / (a1_extra + a2_extra + a3_extra + a4_extra);
                double v_extra = (a1_extra * v1_extra + a2_extra * v2_extra + a3_extra * v3_extra + a4_extra * v4_extra) / (a1_extra + a2_extra + a3_extra + a4_extra);
                double p_extra = (a1_extra * p1_extra + a2_extra * p2_extra + a3_extra * p3_extra + a4_extra * p4_extra) / (a1_extra + a2_extra + a3_extra + a4_extra);

                double u_n_extra = u_extra * n_x + v_extra * n_y;
                double u_t_extra = u_extra * n_y - v_extra * n_x;

                double rho_mirror, p_mirror, u_t_mirror, u_n_mirror;
                if (std::abs(XYCOORD[Index_Coord(i, j, 2 )]) > delta)
                {
                    rho_mirror = rho_image;
                    p_mirror = p_image;
                    u_t_mirror = u_t_image;
                    u_n_mirror = -u_n_image;
                }
                else
                {
                    rho_mirror = (rho_image * ((delta + R_E) * (delta + R_E) - R_G * R_G) - rho_extra * (delta * delta - R_G * R_G)) / (R_E * R_E + 2 * delta * R_E);
                    p_mirror = (p_image * ((delta + R_E) * (delta + R_E) - R_G * R_G) - p_extra * (delta * delta - R_G * R_G)) / (R_E * R_E + 2 * delta * R_E);
                    u_t_mirror = (u_t_image * ((delta + R_E) * (delta + R_E) - R_G * R_G) - u_t_extra * (delta * delta - R_G * R_G)) / (R_E * R_E + 2 * delta * R_E);
                    u_n_mirror = -(u_n_image * ((delta + R_E) * (delta + R_E) - R_G * R_G) - u_n_extra * (delta * delta - R_G * R_G)) / (R_E * R_E + 2 * delta * R_E);
                }

                /* u_n = -u_n; */
                double u_mirror = u_n_mirror * n_x + u_t_mirror * n_y;
                double v_mirror = u_n_mirror * n_y - u_t_mirror * n_x;

                U_OLD[Index(i, j, 0)] = rho_mirror;
                U_OLD[Index(i, j, 1)] = rho_mirror * u_mirror;
                U_OLD[Index(i, j, 2)] = rho_mirror * v_mirror;
                U_OLD[Index(i, j, 3)] = p_mirror / (gamma - 1) + 0.5 * rho_mirror * (u_mirror * u_mirror + v_mirror * v_mirror);
            }
        }
    }

}