#include "Scheme_Index.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
void Scheme_Index(const int N_x,
                  const int N_y,
                  const int num_ghost_cell,
                  double *XYCOORD,
                  double *SCHEME_IDX)
{

#pragma acc parallel loop
    for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
    {
#pragma acc loop
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {

            SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] = 0; // x-direction index
            SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] = 0; // y-direction index
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] > 1)
            {

                for (int k = 0; k < num_ghost_cell + 1; k++)
                {
                    if (XYCOORD[Index_Coord(i - k, j, 5, N_x + 2 * num_ghost_cell)] < 1) // left search
                    {
                        SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] = 1 - 2 * k;
                        break;
                    }
                }
                if (SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] == 0)
                {
                    for (int k = 0; k < num_ghost_cell + 1; k++)
                    {
                        if (XYCOORD[Index_Coord(i + k, j, 5, N_x + 2 * num_ghost_cell)] < 1) // right search
                        {
                            SCHEME_IDX[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)] = 2 * k - 1;
                            break;
                        }
                    }
                }

                for (int k = 0; k < num_ghost_cell + 1; k++)
                {
                    if (XYCOORD[Index_Coord(i, j - k, 5, N_x + 2 * num_ghost_cell)] < 1) // down search
                    {
                        SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] = 1 - 2 * k;
                        break;
                    }
                }
                if (SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] == 0)
                {
                    for (int k = 0; k < num_ghost_cell + 1; k++)
                    {
                        if (XYCOORD[Index_Coord(i, j + k, 5, N_x + 2 * num_ghost_cell)] < 1) // right search
                        {
                            SCHEME_IDX[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)] = 2 * k - 1;
                            break;
                        }
                    }
                }
            }
        }
    }
};