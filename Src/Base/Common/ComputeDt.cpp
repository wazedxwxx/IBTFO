// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn
#include "ComputeDt.H"
void ComputeDt(const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               const double gamma,
               const double CFL_number,
               double *U_OLD,
               double *XYCOORD,
               double *TMP_DT,
               const int ndevices,
               const int device)
{
        int h = (N_y + 2 * num_ghost_cell) / ndevices;
        int lower = max(h * device + 1 - 2 * num_ghost_cell, 0);
        int upper = min(h * (device + 1) + 2 * num_ghost_cell, N);
#pragma acc data present(XYCOORD[(M) * lower * num_coord:(M) * (upper-lower) * num_coord], \
                 U_OLD[(M) * lower * num_eq:(M) * (upper-lower) * num_eq],TMP_DT[device:1])
        {
            double a_max = 0.0;
            double u_max = 0.0;
            double v_max = 0.0;
            const double dx = Psy_L / N_x;
            const double dy = Psy_H / N_y;
#pragma acc parallel loop reduction(max: a_max, u_max, v_max)
            for (int j = lower + num_ghost_cell; j < upper - num_ghost_cell; j++)
            {
#pragma acc loop
                for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
                {
                    if (XYCOORD[Index_Coord(i, j, 5)] < 1)
                    {

                        double a = 0.0;
                        double rho = U_OLD[Index(i, j, 0)];
                        double u = U_OLD[Index(i, j, 1)] / U_OLD[Index(i, j, 0)];
                        double v = U_OLD[Index(i, j, 2)] / U_OLD[Index(i, j, 0)];
                        double p = (gamma - 1) * (U_OLD[Index(i, j, 3)] - 0.5 * rho * (u * u + v * v));
                        a = std::pow((gamma * p / rho), 0.5);
                        a_max = max(a_max, a);
                        u_max = max(u_max, abs(u));
                        v_max = max(v_max, abs(v));
                    }
                }
            }
            TMP_DT[device] = CFL_number * dx / (a_max + max(u_max, v_max));
        }
}
