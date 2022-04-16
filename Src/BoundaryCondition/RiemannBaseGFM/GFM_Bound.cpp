// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Boundary.H"
void Boundary(const int N_x,
              const int N_y,
              const int num_ghost_cell,
              const double gamma,
              double *U_OLD,
              double *U_NEW,
              double *XYCOORD,
              double *GFM_Index,
              const int ndevices,
              const int device)
{
    int lower = LOWER;
    int upper = UPPER;
#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord]) \
    present(GFM_Index [(M)*lower * num_sch:(M) * (upper - lower) * num_sch])                \
        present(U_OLD [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])                  \
            present(U_NEW [(M)*lower * num_eq:(M) * (upper - lower) * num_eq])
    {
#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
#pragma acc loop
                for (int k = 0; k < num_eq; k++)
                {
                    if (XYCOORD[Index_Coord(i, j, 5)] == 0)
                        U_OLD[Index(i, j, k)] = U_NEW[Index(i, j, k)];
                }
            }
        }

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
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
#pragma acc parallel loop async
        for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {

                if (XYCOORD[Index_Coord(i, j, 5)] > 1)
                {

                    double a1 = 0;
                    double a2 = 0;
                    double a3 = 0;
                    double a4 = 0;
                    double n_x = XYCOORD[Index_Coord(i, j, 3)];
                    double n_y = XYCOORD[Index_Coord(i, j, 4)];
                    int IDX = GFM_Index[Index_sch(i, j, 2)];
                    int IDY = GFM_Index[Index_sch(i, j, 3)];
                    double dx = XYCOORD[Index_Coord(1, lower + 1, 0)] - XYCOORD[Index_Coord(0, lower + 1, 0)];
                    double dy = XYCOORD[Index_Coord(1, lower + 2, 1)] - XYCOORD[Index_Coord(1, lower + 1, 1)];
                    double rho_weight = 0;
                    double u_weight = 0;
                    double v_weight = 0;
                    double p_weight = 0;

                    double rho1 = U_NEW[Index(IDX, IDY, 0)];
                    double u1 = U_NEW[Index(IDX, IDY, 1)] / rho1;
                    double v1 = U_NEW[Index(IDX, IDY, 2)] / rho1;
                    double p1 = (gamma - 1) * (U_NEW[Index(IDX, IDY, 3)] - 0.5 * rho1 * (u1 * u1 + v1 * v1));
                    double rho2 = U_NEW[Index(IDX, IDY + 1, 0)];
                    double u2 = U_NEW[Index(IDX, IDY + 1, 1)] / rho2;
                    double v2 = U_NEW[Index(IDX, IDY + 1, 2)] / rho2;
                    double p2 = (gamma - 1) * (U_NEW[Index(IDX, IDY + 1, 3)] - 0.5 * rho2 * (u2 * u2 + v2 * v2));
                    double rho3 = U_NEW[Index(IDX + 1, IDY + 1, 0)];
                    double u3 = U_NEW[Index(IDX + 1, IDY + 1, 1)] / rho3;
                    double v3 = U_NEW[Index(IDX + 1, IDY + 1, 2)] / rho3;
                    double p3 = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY + 1, 3)] - 0.5 * rho3 * (u3 * u3 + v3 * v3));
                    double rho4 = U_NEW[Index(IDX + 1, IDY, 0)];
                    double u4 = U_NEW[Index(IDX + 1, IDY, 1)] / rho4;
                    double v4 = U_NEW[Index(IDX + 1, IDY, 2)] / rho4;
                    double p4 = (gamma - 1) * (U_NEW[Index(IDX + 1, IDY, 3)] - 0.5 * rho4 * (u4 * u4 + v4 * v4));

                    double x_d = GFM_Index[Index_sch(i, j, 0)] - XYCOORD[Index_Coord(IDX, IDY, 0)];
                    double y_d = GFM_Index[Index_sch(i, j, 0)] - XYCOORD[Index_Coord(IDX, IDY, 0)];
                    double rhox_1 = 0;
                    double ux_1 = 0;
                    double vx_1 = 0;
                    double px_1 = 0;
                    double rhox_2 = 0;
                    double ux_2 = 0;
                    double vx_2 = 0;
                    double px_2 = 0;

                    // case 1
                    //              *2     *3
                    //
                    //              *1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {

                        if (x_d > dx / 2)
                        {
                            rhox_1 = rho4;
                            ux_1 = u4;
                            vx_1 = v4;
                            px_1 = p4;
                            rhox_2 = rho3;
                            ux_2 = u3;
                            vx_2 = v3;
                            px_2 = p3;
                        }
                        else
                        {
                            rhox_1 = rho1;
                            ux_1 = u1;
                            vx_1 = v1;
                            px_1 = p1;
                            rhox_2 = rho2;
                            ux_2 = u2;
                            vx_2 = v2;
                            px_2 = p2;
                        }
                        if (y_d > dy / 2)
                        {

                            rho_weight = rhox_2;
                            u_weight = ux_2;
                            v_weight = vx_2;
                            p_weight = px_2;
                        }
                        else
                        {
                            rho_weight = rhox_1;
                            u_weight = ux_1;
                            v_weight = vx_1;
                            p_weight = px_1;
                        }
                    }

                    // case 2
                    //              *2     *3
                    //
                    //              o1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {
                        rhox_1 = rho4;
                        ux_1 = u4;
                        vx_1 = v4;
                        px_1 = p4;
                        if (x_d > dx / 2)
                        {

                            rhox_2 = rho3;
                            ux_2 = u3;
                            vx_2 = v3;
                            px_2 = p3;
                        }
                        else
                        {
                            rhox_2 = rho2;
                            ux_2 = u2;
                            vx_2 = v2;
                            px_2 = p2;
                        }
                        if (y_d > dy / 2)
                        {

                            rho_weight = rhox_2;
                            u_weight = ux_2;
                            v_weight = vx_2;
                            p_weight = px_2;
                        }
                        else
                        {
                            rho_weight = rhox_1;
                            u_weight = ux_1;
                            v_weight = vx_1;
                            p_weight = px_1;
                        }
                    }

                    // case 3
                    //              *2     *3
                    //
                    //              *1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {
                        rhox_1 = rho1;
                        ux_1 = u1;
                        vx_1 = v1;
                        px_1 = p1;
                        if (x_d > dx / 2)
                        {

                            rhox_2 = rho3;
                            ux_2 = u3;
                            vx_2 = v3;
                            px_2 = p3;
                        }
                        else
                        {
                            rhox_2 = rho2;
                            ux_2 = u2;
                            vx_2 = v2;
                            px_2 = p2;
                        }
                        if (y_d > dy / 2)
                        {

                            rho_weight = rhox_2;
                            u_weight = ux_2;
                            v_weight = vx_2;
                            p_weight = px_2;
                        }
                        else
                        {
                            rho_weight = rhox_1;
                            u_weight = ux_1;
                            v_weight = vx_1;
                            p_weight = px_1;
                        }
                    }

                    // case 4
                    //              o2     *3
                    //
                    //              *1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {
                        rhox_2 = rho3;
                        ux_2 = u3;
                        vx_2 = v3;
                        px_2 = p3;
                        if (x_d > dx / 2)
                        {
                            rhox_1 = rho4;
                            ux_1 = u4;
                            vx_1 = v4;
                            px_1 = p4;
                        }
                        else
                        {
                            rhox_1 = rho1;
                            ux_1 = u1;
                            vx_1 = v1;
                            px_1 = p1;
                        }
                        if (y_d > dy / 2)
                        {

                            rho_weight = rhox_2;
                            u_weight = ux_2;
                            v_weight = vx_2;
                            p_weight = px_2;
                        }
                        else
                        {
                            rho_weight = rhox_1;
                            u_weight = ux_1;
                            v_weight = vx_1;
                            p_weight = px_1;
                        }
                    }

                    // case 5
                    //              *2     o3
                    //
                    //              *1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {
                        rhox_2 = rho2;
                        ux_2 = u2;
                        vx_2 = v2;
                        px_2 = p2;
                        if (x_d > dx / 2)
                        {
                            rhox_1 = rho4;
                            ux_1 = u4;
                            vx_1 = v4;
                            px_1 = p4;
                        }
                        else
                        {
                            rhox_1 = rho1;
                            ux_1 = u1;
                            vx_1 = v1;
                            px_1 = p1;
                        }
                        if (y_d > dy / 2)
                        {

                            rho_weight = rhox_2;
                            u_weight = ux_2;
                            v_weight = vx_2;
                            p_weight = px_2;
                        }
                        else
                        {
                            rho_weight = rhox_1;
                            u_weight = ux_1;
                            v_weight = vx_1;
                            p_weight = px_1;
                        }
                    }

                    // case 6
                    //              o2     o3
                    //
                    //              *1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        if (x_d > dx / 2)
                        {
                            rho_weight = rho4;
                            u_weight = u4;
                            v_weight = v4;
                            p_weight = p4;
                        }
                        else
                        {
                            rho_weight = rho1;
                            u_weight = u1;
                            v_weight = v1;
                            p_weight = p1;
                        }
                    }

                    // case 7
                    //              *2     *3
                    //
                    //              o1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {

                        if (x_d > dx / 2)
                        {
                            rho_weight = rho3;
                            u_weight = u3;
                            v_weight = v3;
                            p_weight = p3;
                        }
                        else
                        {
                            rho_weight = rho2;
                            u_weight = u2;
                            v_weight = v2;
                            p_weight = p2;
                        }
                    }

                    // case 8
                    //              *2     o3
                    //
                    //              *1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        if (y_d > dy / 2)
                        {

                            rho_weight = rho2;
                            u_weight = u2;
                            v_weight = v2;
                            p_weight = p2;
                        }
                        else
                        {
                            rho_weight = rho1;
                            u_weight = u1;
                            v_weight = v1;
                            p_weight = p1;
                        }
                    }

                    // case 9
                    //              o2     *3
                    //
                    //              o1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {

                        if (y_d > dy / 2)
                        {

                            rho_weight = rho3;
                            u_weight = u3;
                            v_weight = v3;
                            p_weight = p3;
                        }
                        else
                        {
                            rho_weight = rho4;
                            u_weight = u4;
                            v_weight = v4;
                            p_weight = p4;
                        }
                    }

                    // case 10
                    //              o2     *3
                    //
                    //              *1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {

                        if (y_d > dy / 2)
                        {

                            rho_weight = rho3;
                            u_weight = u3;
                            v_weight = v3;
                            p_weight = p3;
                        }
                        else
                        {
                            rho_weight = rho1;
                            u_weight = u1;
                            v_weight = v1;
                            p_weight = p1;
                        }
                    }

                    // case 11
                    //              *2     o3
                    //
                    //              o1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        if (y_d > dy / 2)
                        {

                            rho_weight = rho2;
                            u_weight = u2;
                            v_weight = v2;
                            p_weight = p2;
                        }
                        else
                        {
                            rho_weight = rho4;
                            u_weight = u4;
                            v_weight = v4;
                            p_weight = p4;
                        }
                    }

                    // case 12
                    //              o2     o3
                    //
                    //              *1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        rho_weight = rho1;
                        u_weight = u1;
                        v_weight = v1;
                        p_weight = p1;
                    }

                    // case 13
                    //              *2     o3
                    //
                    //              o1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        rho_weight = rho2;
                        u_weight = u2;
                        v_weight = v2;
                        p_weight = p2;
                    }

                    // case 14
                    //              o2     *3
                    //
                    //              o1     o4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] < 1)
                    {

                        rho_weight = rho3;
                        u_weight = u3;
                        v_weight = v3;
                        p_weight = p3;
                    }

                    // case 15
                    //              o2     o3
                    //
                    //              o1     *4

                    if (XYCOORD[Index_Coord(IDX, IDY, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY, 5)] < 1 &&
                        XYCOORD[Index_Coord(IDX, IDY + 1, 5)] > 0 &&
                        XYCOORD[Index_Coord(IDX + 1, IDY + 1, 5)] > 0)
                    {

                        rho_weight = rho4;
                        u_weight = u4;
                        v_weight = v4;
                        p_weight = p4;
                    }

                    double u_n = u_weight * n_x + v_weight * n_y;
                    double u_t = u_weight * n_y - v_weight * n_x;
                    u_n = -u_n;
                    u_weight = u_n * n_x + u_t * n_y;
                    v_weight = u_n * n_y - u_t * n_x;

                    U_OLD[Index(i, j, 0)] = rho_weight;
                    U_OLD[Index(i, j, 1)] = rho_weight * u_weight;
                    U_OLD[Index(i, j, 2)] = rho_weight * v_weight;
                    U_OLD[Index(i, j, 3)] = p_weight / (gamma - 1) + 0.5 * rho_weight * (u_weight * u_weight + v_weight * v_weight);

                    /*                                       std::cout << " i " << i << " j " << j
                                                                  << " x " << XYCOORD[Index_Coord(i, j, 0 )]
                                                                  << " y " << XYCOORD[Index_Coord(i, j, 1 )]
                                                                  << " Phi " << std::abs(XYCOORD[Index_Coord(i, j, 2 )])
                                                                  << " nx " << XYCOORD[Index_Coord(i, j, 3 )]
                                                                  << " ny " << XYCOORD[Index_Coord(i, j, 4 )]
                                                                  << " MirrorIDX " << GFM_Index[Index_sch(i, j, 2 )]
                                                                  << " MirrorIDY " << GFM_Index[Index_sch(i, j, 3 )]
                                                                  << " A1 " << a1
                                                                  << " rho1 " << rho1
                                                                  << " u1 " << u1
                                                                  << " v1 " << v1
                                                                  << " p1 " << p1
                                                                  << " A2 " << a2
                                                                  << " rho2 " << rho2
                                                                  << " u2 " << u2
                                                                  << " v2 " << v2
                                                                  << " p2 " << p2
                                                                  << " A3 " << a3
                                                                  << " rho3 " << rho3
                                                                  << " u3 " << u3
                                                                  << " v3 " << v3
                                                                  << " p3 " << p3
                                                                  << " A4 " << a4
                                                                  << " rho4 " << rho4
                                                                  << " u4 " << u4
                                                                  << " v4 " << v4
                                                                  << " p4 " << p4
                                                                  << " u_n " << u_n
                                                                  << " u_t " << u_t
                                                                  << " rhow " << rho_weight
                                                                  << " uw " << u_weight
                                                                  << " vw " << v_weight
                                                                  << " pw " << p_weight
                                                                  << " " << std::endl;   */
                }
            }
        }
    }
}
