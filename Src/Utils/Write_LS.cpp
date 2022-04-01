#include <fstream>
#include <iostream>
#include "Write_LS.H"
#include "CoordDefine.H"
#define Index_sch(a, b, c, N) ((N) * (b) + (a)) * 2 + (c)
using namespace std;
void Write_LS(const double Psy_L,
              const double Psy_H,
              const int N_x,
              const int N_y,
              const int num_ghost_cell,
              double *XYCOORD)
{
    double Ori = -Psy_L * num_ghost_cell / N_x;
    ofstream outfile;
    outfile.open("gemo_ls.vtk", ios::out | ios::trunc);
    outfile << "# vtk DataFile Version 2.0" << endl;
    outfile << "Volume example" << endl;
    outfile << "ASCII" << endl;
    outfile << "DATASET STRUCTURED_POINTS" << endl;
    outfile << "DIMENSIONS"
            << " " << N_x + 2 * num_ghost_cell << " "
            << " " << N_y + 2 * num_ghost_cell << " " << 1 << endl;
    outfile << "ASPECT_RATIO 1 1 1" << endl;
    outfile << "ORIGIN " << Ori << " " << Ori << " " << Ori << endl;
    outfile << "SPACING " << Psy_L / N_x << " " << Psy_H / N_y << " " << Psy_H / N_y << " " << endl;
    outfile << "POINT_DATA " << (N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) << endl;
    outfile << "SCALARS Phi float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << XYCOORD[Index_Coord(i, j, 2, N_x + 2 * num_ghost_cell)] << endl;
        }
    }

    outfile << "SCALARS nx float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << XYCOORD[Index_Coord(i, j, 3, N_x + 2 * num_ghost_cell)] << endl;
        }
    }

    outfile << "SCALARS ny float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << XYCOORD[Index_Coord(i, j, 4, N_x + 2 * num_ghost_cell)] << endl;
        }
    }

    outfile << "SCALARS cell_type float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << XYCOORD[Index_Coord(i, j, 5, N_x + 2 * num_ghost_cell)] << endl;
        }
    }

    outfile.close();
}