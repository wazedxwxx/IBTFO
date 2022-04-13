#include "Level_Set_function.H"

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
    double f0 = para.get("f0", 0);
    double f1 = para.get("f1", 1);
    const double rvortex = 0.40;
    const double rl0 = f0 * rvortex;
    const double rl1 = f1 * rvortex;

#pragma acc data present(XYCOORD [(M)*lower * num_coord:(M) * (upper - lower) * num_coord])
    {
#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                double phi1, phi2;
                double xlc = 0.50;
                double ylc = 0.50;
                double xp = XYCOORD[Index_Coord(i, j, 0)];
                double yp = XYCOORD[Index_Coord(i, j, 1)];
                double r = 0;
                r = sqrt((xp - xlc) * (xp - xlc) + (yp - ylc) * (yp - ylc));

                phi1 = rl1 - r;
                phi2 = r - rl0;

                XYCOORD[Index_Coord(i, j, 2)] = ((phi1 < phi2) ? phi1 : phi2);
            }
        }
    }
}
