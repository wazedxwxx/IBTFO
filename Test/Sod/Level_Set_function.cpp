#include <math.h>
#include "Level_Set_function.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
using namespace std;


void Level_Set_function(char *filename,
                          const int N_x,
                          const int N_y,
                          const int num_ghost_cell,
                          double *XYCOORD,
                          const int ndevices,
                          const int device)
{
    int lower = LOWER;
    int upper = UPPER;

    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    double angle = 3.1415926535 * para.get("angle", 30) / 180;
    double AxisD = para.get("AxisD", 0.4);

#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord])
    {
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {

                XYCOORD[Index_Coord(i, j, 2)] = 1;
            }
        }


    }
}