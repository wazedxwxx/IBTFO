#include <fstream>
#include <iostream>
#include "WriteIDX2TXT.H"
#define Index_sch(a, b, c, N) ((N) * (b) + (a)) * 2 + (c)
using namespace std;
void WriteIDX2TXT(const int N_x,
               const int N_y,
               const int num_ghost_cell,
               int *U)
               {
    ofstream outfile1;
    outfile1.open("X_index.txt", ios::out | ios::trunc);
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
        for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
        {
            outfile1 << U[Index_sch(i, j, 0, N_x + 2 * num_ghost_cell)]<<" " ;
        }
        outfile1 << endl;
    }
    outfile1.close();

    ofstream outfile2;
    outfile2.open("Y_index.txt", ios::out | ios::trunc);
    for (int j = num_ghost_cell; j < N_y + num_ghost_cell; j++)
    {
        for (int i = num_ghost_cell; i < N_x + num_ghost_cell; i++)
        {
            outfile2 << U[Index_sch(i, j, 1, N_x + 2 * num_ghost_cell)]<<" " ;
        }
        outfile2 << endl;
    }
    outfile2.close();


               }