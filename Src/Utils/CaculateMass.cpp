// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "CaculateMass.H"
void CaculateMass(const double Psy_L,
                  const double Psy_H,
                  const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *U_OLD,
                  double *XYCOORD,
                  double *MASS_DEVICE,
                  const int ndevices,
                  const int device)
{

    int real_lower = REAL_LOWER;
    int real_upper = REAL_UPPER;
    acc_set_device_num(device, acc_device_default);
#pragma acc data present(U_OLD [(M)*real_lower * num_eq:(M) * (real_upper - real_lower) * num_eq], \
                         XYCOORD [(M)*real_lower * num_coord:(M) * (real_upper - real_lower) * num_coord])
    {
        double mass = 0;
        double dx = Psy_L / N_x;
        double dy = Psy_H / N_y;
#pragma acc update self(U_OLD [(M)*real_lower * num_eq:(M) * (real_upper - real_lower) * num_eq])
#pragma acc parallel loop
        for (int j = real_lower; j < real_upper; j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                {

                    mass += U_OLD[Index(i, j, 0)] * dx * dy;
                }
            }
        }
        MASS_DEVICE[device] = mass;
    }
}
