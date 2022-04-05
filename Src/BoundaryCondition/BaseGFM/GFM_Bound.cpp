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
                if (XYCOORD[Index_Coord(i, j, 5)] == 0)
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
#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {

            if (XYCOORD[Index_Coord(i, j, 5)] > 1)
            {
                double a1 = GFM_Index[Index_sch(i, j, 4)];
                double a2 = GFM_Index[Index_sch(i, j, 5)];
                double a3 = GFM_Index[Index_sch(i, j, 6)];
                double a4 = GFM_Index[Index_sch(i, j, 7)];
                double n_x = XYCOORD[Index_Coord(i, j, 3)];
                double n_y = XYCOORD[Index_Coord(i, j, 4)];
                int IDX = GFM_Index[Index_sch(i, j, 2)];
                int IDY = GFM_Index[Index_sch(i, j, 3)];
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

                double rho_weight = (a1 * rho1 + a2 * rho2 + a3 * rho3 + a4 * rho4) / (a1 + a2 + a3 + a4);
                double u_weight = (a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4) / (a1 + a2 + a3 + a4);
                double v_weight = (a1 * v1 + a2 * v2 + a3 * v3 + a4 * v4) / (a1 + a2 + a3 + a4);
                double p_weight = (a1 * p1 + a2 * p2 + a3 * p3 + a4 * p4) / (a1 + a2 + a3 + a4);

                double u_n = u_weight * n_x + v_weight * n_y;
                double u_t = u_weight * n_y - v_weight * n_x;
                u_n = -u_n;
                u_weight = u_n * n_x + u_t * n_y;
                v_weight = u_n * n_y - u_t * n_x;

                U_OLD[Index(i, j, 0)] = rho_weight;
                U_OLD[Index(i, j, 1)] = rho_weight * u_weight;
                U_OLD[Index(i, j, 2)] = rho_weight * v_weight;
                U_OLD[Index(i, j, 3)] = p_weight / (gamma - 1) + 0.5 * rho_weight * (u_weight * u_weight + v_weight * v_weight);

/*                  std::cout << " i " << i << " j " << j
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
                          << " " << std::endl;  */
            }
        }
    }
}


