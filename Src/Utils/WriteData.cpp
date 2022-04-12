#include <math.h>
#include "WriteData.H"
#include "EQDefine.H"
#include <string>
#include "CoordDefine.H"
#include <iostream>
#include <stdio.h>
#include <string>
#include "openacc.h"
using namespace std;
void WriteData(const double lo_x,
               const double lo_y,
               const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               int iter,
               double now_t,
               const double gamma,
               double *U_OLD,
               double *XYCOORD,
               const int ndevices,
               const int device)
{
    FILE *header;
    header = fopen("result.visit", "a");
    if (iter == 0 && device == 0)
    {
        fprintf(header, "!NBLOCKS %d\n",ndevices);
    }
        fprintf(header, "dat_%d_%d.vtk\n",device,iter);
    fclose(header);

    int real_lower = REAL_LOWER;
    int real_upper = REAL_UPPER;
acc_set_device_num(device, acc_device_default);
#pragma acc data present(U_OLD [(M)*real_lower * num_eq:(M) * (real_upper - real_lower) * num_eq], \
                         XYCOORD [(M)*real_lower * num_coord:(M) * (real_upper - real_lower) * num_coord])
    {

        FILE *out;
        string filename = "dat_" + to_string(device) + "_" + to_string(iter) + ".vtk";
        out = fopen(filename.c_str(), "w");
        fprintf(out, "# vtk DataFile Version 2.0 \n");
        fprintf(out, "Volume example \n");
        fprintf(out, "ASCII \n");
        fprintf(out, "DATASET STRUCTURED_POINTS \n");
        fprintf(out, "DIMENSIONS %d %d %d  \n", N_x, min(real_upper + 1, N_x + num_ghost_cell) - max(real_lower, num_ghost_cell), 1);
        fprintf(out, "ASPECT_RATIO 1 1 1 \n");
        fprintf(out, "ORIGIN %f %f %f \n", lo_x + Psy_L / N_x / 2, lo_y + Psy_H / N_y / 2 + (max(real_lower, num_ghost_cell) - num_ghost_cell) * Psy_H / N_y, 0.0);
        fprintf(out, "SPACING %f %f %f \n", Psy_L / N_x, Psy_H / N_y, Psy_H / N_y);
        fprintf(out, "FIELD FieldData 1  \n");
        fprintf(out, "TIME 1 1 double  \n");
        fprintf(out, "%f  \n", now_t);
        fprintf(out, "POINT_DATA  %d\n", N_x * (min(real_upper + 1, N_x + num_ghost_cell) - max(real_lower, num_ghost_cell)));
        fprintf(out, "SCALARS density float  \n");
        fprintf(out, "LOOKUP_TABLE default  \n");

#pragma acc update self(U_OLD [(M)*real_lower * num_eq:(M) * (real_upper - real_lower) * num_eq])
#pragma acc loop seq
        for (int j = max(real_lower, num_ghost_cell); j < min(real_upper + 1, N_y + num_ghost_cell); j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                    if (U_OLD[Index(i, j, 0)] > 1e-6)
                        fprintf(out, "%f  \n", U_OLD[Index(i, j, 0)]);
                    else
                        fprintf(out, "%f  \n", 1e-6);
                else
                    fprintf(out, "%f  \n", -1.0);
            }
        }

        fprintf(out, "SCALARS x-vel float \n");
        fprintf(out, "LOOKUP_TABLE default  \n");

#pragma acc loop seq
        for (int j = max(real_lower, num_ghost_cell); j < min(real_upper + 1, N_y + num_ghost_cell); j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                    fprintf(out, "%f  \n", U_OLD[Index(i, j, 1)] / U_OLD[Index(i, j, 0)]);
                else
                    fprintf(out, "%f  \n", -100.0);
            }
        }

        fprintf(out, "SCALARS y-vel float  \n");
        fprintf(out, "LOOKUP_TABLE default  \n");

#pragma acc loop seq
        for (int j = max(real_lower, num_ghost_cell); j < min(real_upper + 1, N_y + num_ghost_cell); j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                    fprintf(out, "%f  \n", U_OLD[Index(i, j, 2)] / U_OLD[Index(i, j, 0)]);
                else
                    fprintf(out, "%f  \n", -100.0);
            }
        }

        fprintf(out, "SCALARS pressure float  \n");
        fprintf(out, "LOOKUP_TABLE default  \n");

#pragma acc loop seq
        for (int j = max(real_lower, num_ghost_cell); j < min(real_upper + 1, N_y + num_ghost_cell); j++)
        {
#pragma acc loop
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                {
                    double rho = U_OLD[Index(i, j, 0)];
                    double u = U_OLD[Index(i, j, 1)] / U_OLD[Index(i, j, 0)];
                    double v = U_OLD[Index(i, j, 2)] / U_OLD[Index(i, j, 0)];
                    double p = (gamma - 1) * (U_OLD[Index(i, j, 3)] - 0.5 * rho * (u * u + v * v));
                    fprintf(out, "%f  \n", p);
                }
                else
                {
                    fprintf(out, "%f  \n", -1.0);
                }
            }
        }

        fclose(out);
    }

    /*     outfile << "SCALARS x-vel float" << endl;
        outfile << "LOOKUP_TABLE default" << endl;
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                    outfile << U_OLD[Index(i, j, 1)] / U_OLD[Index(i, j, 0)] << endl;
                else
                    outfile << -100.0 << endl;
            }
        }

        outfile << "SCALARS y-vel float" << endl;
        outfile << "LOOKUP_TABLE default" << endl;
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                    outfile << U_OLD[Index(i, j, 2)] / U_OLD[Index(i, j, 0)] << endl;
                else
                    outfile << -100.0 << endl;
            }
        }

        outfile << "SCALARS pressure float" << endl;
        outfile << "LOOKUP_TABLE default" << endl;
        for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
        {
            for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
            {
                if (XYCOORD[Index_Coord(i, j, 5)] < 0.5)
                {
                    double rho = U_OLD[Index(i, j, 0)];
                    double u = U_OLD[Index(i, j, 1)] / U_OLD[Index(i, j, 0)];
                    double v = U_OLD[Index(i, j, 2)] / U_OLD[Index(i, j, 0)];
                    double p = (gamma - 1) * (U_OLD[Index(i, j, 3)] - 0.5 * rho * (u * u + v * v));
                    outfile << p << endl;
                }
                else
                    outfile << -1.0 << endl;
            }
        } */
}
