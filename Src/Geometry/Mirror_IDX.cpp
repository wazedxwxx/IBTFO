#include "Mirror_IDX.H"
#include <math.h>
#include "CoordDefine.H"
void Mirror_IDX(const int N_x,
                const int N_y,
                const int num_ghost_cell,
                int *IDX,
                int *IDY,
                double *XYCOORD,
                const double x_mirror,
                const double y_mirror)
{
// x direction search
#pragma acc parallel loop
    for (int i = 0; i < N_x + 2 * num_ghost_cell - 1; i++)
    {
        if (x_mirror < XYCOORD[Index_Coord(i + 1, 0, 0)] && x_mirror >= XYCOORD[Index_Coord(i, 0, 0)])
        {
            *IDX = i;
        }
    }

#pragma acc parallel loop
    for (int j = 0; j < N_y + 2 * num_ghost_cell - 1; j++)
    {
        if (y_mirror < XYCOORD[Index_Coord(0, j + 1, 1)] && y_mirror >= XYCOORD[Index_Coord(0, j, 1)])
        {
            *IDY = j;
        }
    }
}