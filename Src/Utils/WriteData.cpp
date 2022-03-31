#include <fstream>
#include <iostream>
#include "WriteData.H"
#define num_eq 4
#define Index(a, b, c, N) ((N) * (b) + (a)) * num_eq + (c)
#define Index_Coord(a, b, c, N) ((N) * (b) + (a)) * 6 + (c)
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
               double *XYCOORD)
{
    ofstream outfile;
    outfile.open("dat_" + to_string(iter) + ".vtk", ios::out | ios::trunc);
    outfile << "# vtk DataFile Version 2.0" << endl;
    outfile << "Volume example" << endl;
    outfile << "ASCII" << endl;
    outfile << "DATASET STRUCTURED_POINTS" << endl;
    outfile << "DIMENSIONS"
            << " " << N_x << " "
            << " " << N_y << " " << 1 << endl;
    outfile << "ASPECT_RATIO 1 1 1" << endl;
    outfile << "ORIGIN " << lo_x << " " << lo_y << " 0 " << endl;
    outfile << "SPACING " << Psy_L / N_x << " " << Psy_H / N_y << " " << Psy_H / N_y << " " << endl;
    outfile << "FIELD FieldData 1" << endl;
    outfile << "TIME 1 1 double " << endl;
    outfile << now_t << endl;
    outfile << "POINT_DATA " << N_x * N_y << endl;
    outfile << "SCALARS density float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
        for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5)
                if (U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] > 1e-6)
                    outfile << U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] << endl;
                else
                    outfile << 1e-6 << endl;
            else
                outfile << -1.0 << endl;
        }
    }

    outfile << "SCALARS x-vel float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
        for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
        {
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5)
                outfile << U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] << endl;
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
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5)
                outfile << U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)] << endl;
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
            if (XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] < 0.5)
            {
                double rho = U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double u = U_OLD[Index(i, j, 1, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double v = U_OLD[Index(i, j, 2, N_x + 2 * num_ghost_cell)] / U_OLD[Index(i, j, 0, N_x + 2 * num_ghost_cell)];
                double p = (gamma - 1) * (U_OLD[Index(i, j, 3, N_x + 2 * num_ghost_cell)] - 0.5 * rho * (u * u + v * v));
                outfile << p << endl;
            }
            else
                outfile << -1.0 << endl;
        }
    }

    outfile.close();
}