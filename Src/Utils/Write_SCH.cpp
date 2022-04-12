#include <fstream>
#include <iostream>
#include "Write_SCH.H"
#include "SchDefine.H"
using namespace std;
void Write_SCH(const double Psy_L,
               const double Psy_H,
               const int N_x,
               const int N_y,
               const int num_ghost_cell,
               double *SCHEME_IDX)
{
        double Ori = -Psy_L * num_ghost_cell / N_x;
    ofstream outfile;
    outfile.open("gemo_sch.vtk", ios::out | ios::trunc);
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

    outfile << "SCALARS Iamge_X_Coord float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 0)] << endl;
        }
    }

    outfile << "SCALARS Iamge_Y_Coord float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 1)] << endl;
        }
    }


    outfile << "SCALARS Iamge_MIRROR_IDX float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 2)] << endl;
        }
    }

    outfile << "SCALARS Iamge_MIRROR_IDY float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 3)] << endl;
        }
    }

    outfile << "SCALARS Extra_X_Coord float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 8)] << endl;
        }
    }

    outfile << "SCALARS Extra_Y_Coord float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 9)] << endl;
        }
    }


    outfile << "SCALARS EXTRA_MIRROR_IDX float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 10)] << endl;
        }
    }

    outfile << "SCALARS Extra_MIRROR_IDY float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 11)] << endl;
        }
    } 

    outfile << "SCALARS Image_A1 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 4)] << endl;
        }
    }


    outfile << "SCALARS Image_A2 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 5)] << endl;
        }
    }

    outfile << "SCALARS Image_A3 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 6)] << endl;
        }
    }

    outfile << "SCALARS Image_A4 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 7)] << endl;
        }
    }


    outfile << "SCALARS Extra_A1 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 12)] << endl;
        }
    }


    outfile << "SCALARS Extra_A2 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 13)] << endl;
        }
    }


    outfile << "SCALARS Extra_A3 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 14)] << endl;
        }
    }


    outfile << "SCALARS Extra_A4 float" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int j = 0; j < N_y + 2 * num_ghost_cell; j++)
    {
        for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
        {
            outfile << SCHEME_IDX[Index_sch(i, j, 15)] << endl;
        }
    } 



    outfile.close();
}