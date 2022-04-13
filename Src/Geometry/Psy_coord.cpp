// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#include "Psy_coord.H"
void Psy_coord(const double lo_x,
               const double lo_y,
               const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *XYCOORD,
               const int ndevices,
               const int device)
{
    int lower = LOWER;
    int upper = UPPER;

#pragma acc data present(XYCOORD[(N_x + 2 * num_ghost_cell) * lower * num_coord:(N_x + 2 * num_ghost_cell) * (upper-lower) * num_coord]) 
    {
        const double dx = Psy_L / N_x;
        const double dy = Psy_L / N_x;
        double Li = lo_x;
        double Hi = lo_y;
        Li = -dx * num_ghost_cell + lo_x + dx / 2;
        Hi = -dy * num_ghost_cell + lo_y + dy / 2;

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop 
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                XYCOORD[Index_Coord(i, j, 0)] = Li + i * dx;
                XYCOORD[Index_Coord(i, j, 1)] = Hi + j * dx;
            }
        }
    }
}
