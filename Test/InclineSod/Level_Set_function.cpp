#include <math.h>
#include "Level_Set_function.H"
#include "EQDefine.H"
#include "CoordDefine.H"
#include <iostream>
using namespace std;

#pragma acc routine worker
inline double line(double x, double y, double k, double b, int dir) // dist for y= k x + b
{
    double phi;
    if (y > k * x + b)
    {
        phi = dir * abs(k * x - y + b) / pow(1 + k * k, 0.5);
    }
    else
    {
        phi = -dir * abs(k * x - y + b) / pow(1 + k * k, 0.5);
    }
    return phi;
}

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
                double phi, phi1, phi2;
                phi1 = line(XYCOORD[Index_Coord(i, j, 0)], XYCOORD[Index_Coord(i, j, 1)], tan(angle), 0, 1);
                phi2 = line(XYCOORD[Index_Coord(i, j, 0)], XYCOORD[Index_Coord(i, j, 1)], tan(angle), AxisD, -1);

                phi = ((phi1 < phi2) ? phi1 : phi2);
                XYCOORD[Index_Coord(i, j, 2)] = phi;
                // Level_Set_function(gemo_factor, XYCOORD[Index_Coord(i, j, 0)], XYCOORD[Index_Coord(i, j, 1)]); // Phi
            }
        }


    }
}
