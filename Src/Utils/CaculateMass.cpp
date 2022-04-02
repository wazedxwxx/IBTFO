#include "EQDefine.H"
#include "CaculateMass.H"
#include "CoordDefine.H"
double CaculateMass(const double Psy_L,
                    const double Psy_H,
                    const int N_x,
                    const int N_y,
                    const int num_ghost_cell,
                    double *U_OLD,
                    double *XYCOORD)
{
    double mass = 0;
    double dx = Psy_L / N_x;
    double dy = Psy_H / N_y;
#pragma acc parallel loop
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
#pragma acc loop
        for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5)
            {

                mass += U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell, N_y + 2 * num_ghost_cell)] * dx * dy;
            }
        }
    }
    return mass;
}