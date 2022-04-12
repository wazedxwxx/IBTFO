#include "ParamReader.H"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>

struct ReadP
{
    int max_iter = 0;
    double lo_x = 0.0; // x direction grid number in computational domain
    double lo_y = 0.0;

    double hi_x = 0.0;       // x direction grid number in computational domain
    double hi_y = 0.0;       // y direction grid number in computational domain
    int N_x = 0;             // x direction grid number in computational domain
    int N_y = 0;             // y direction grid number in computational domain
    double Psy_time = 0.0;   // physical time
    double CFL_number = 0.0; // CFL_number
    int num_ghost_cell = 0;
    int plot_per = 0;
    double gamma = 0.0; // Gas parameters
    double angle = 0.0;
    double yaxisD = 0;
    int plot_int; // Outputs NUM
};

ReadP ReadPara(char *filename);
