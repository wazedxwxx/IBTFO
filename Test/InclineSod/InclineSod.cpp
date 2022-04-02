#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>

#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include "Initialize.H"
#include "Psy_coord.H"
#include "WriteData.H"
#include "Write_LS.H"
#include "CaculateMass.H"
#include "Scheme_Index.H"
#include "Boundary.H"
#include "Level_Set.H"
#include "Conserve2Flux.H"
#include "Scheme_Index.H"
#include "Advance.H"
#include "TimeAdvance.H"
#include "ComputeDt.H"
#include "ParamReader.H"

using namespace std;

int main(int argc,char** argv)
{


    ParamReader DetectParams;
    //parameter settings model
    char* filename = argv[1];
    Params<double> para(DetectParams.open("sod.inp").numbers());
    const int max_iter = int(para.get("max_iter", 100000));//max time steps

    const double lo_x = para.get("lo_x", 0);//x direction grid number in computational domain
    const double lo_y = para.get("lo_y", 0);//y direction grid number in computational domain

    const double hi_x = para.get("hi_x", 1);//x direction grid number in computational domain
    const double hi_y = para.get("hi_y", 1);//y direction grid number in computational domain

    const int N_x = int(para.get("N_x", 500));//x direction grid number in computational domain
    const int N_y = int(para.get("N_y", 200));//y direction grid number in computational domain
    const double Psy_time = para.get("Psy_time", 1.0); // physical time
    const double CFL_number = para.get("CFL_number", 0.3);//CFL_number
    const int num_ghost_cell = int(para.get("num_ghost_cell", 2));
    const int plot_per = int(para.get("plot_per", 100));
    const double gamma = para.get("gamma", 1.4);//Gas parameters
    const int plot_int = para.get("plot_int", 0);//Outputs NUM
    int output_int = 1;


    cout <<" ==== parameters are read ===="<<endl;

    const double h = 1.0 / N_x;        //cell size
    const double delta_t = 1 / 1000;   //Preset time step
    double now_t = 0;                  //Preset time
    int iter = 0;
    double dt = 0.0;
    const double Psy_L = hi_x - lo_x;
    const double Psy_H = hi_y - lo_y;

    double *U_OLD = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq];
    double *XYCOORD = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_coord];
    double *U_TMP = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_NEW = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_TMPRK = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_L = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_R = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_D = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *U_U = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *F_L = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *F_R = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *G_D = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *G_U = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *F_OLD = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *G_OLD = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_eq];
    double *SCHEME_IDX = new double[(N_x + 2 *num_ghost_cell) * (N_y + 2 *num_ghost_cell) * num_sch];
    double Mass_start, Mass_now, Mass_loss;

    cout <<" ====  memory allocation complete ===="<<endl;

    Psy_coord(lo_x,lo_y,Psy_L, Psy_H, N_x, N_y, num_ghost_cell, XYCOORD);
    Level_Set(filename, N_x, N_y, num_ghost_cell, XYCOORD);
    Write_LS(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, XYCOORD);
    Scheme_Index(N_x, N_y, num_ghost_cell, XYCOORD, SCHEME_IDX);
    // WriteIDX2TXT(N_x, N_y, num_ghost_cell, SCHEME_IDX);
    Initialize(filename, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD);
    Mass_start = CaculateMass(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, U_OLD, XYCOORD);
    cout << " ====  Initial Mass is " <<Mass_start <<" ===="<< endl;
    cout << " ====  Geometry Initialize complete ====" << endl;
    Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, SCHEME_IDX);
    //Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, SCHEME_IDX);
    WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD,XYCOORD);

#pragma acc data copy(U_OLD[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq]) \
                 copy(F_L[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(F_R[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(G_D[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(G_U[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_L[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_R[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_U[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_D[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(F_OLD[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(G_OLD[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_TMP[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(U_NEW[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq])\
                 copy(SCHEME_IDX[:(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_sch])\

    while (now_t < Psy_time && iter < max_iter)
    {    
        ComputeDt(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, CFL_number, U_OLD, XYCOORD, &dt);
        TimeAdvance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, F_OLD, G_OLD, U_TMP, U_NEW, U_TMPRK, F_L, F_R, G_D, G_U, U_L, U_R, U_D, U_U,XYCOORD);
        Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW,XYCOORD,SCHEME_IDX);

        now_t = now_t + dt;
        iter++;   

if (plot_int==0){
    if (iter % plot_per == 0)
    {
        cout <<" ==== Writing Data Wait ===="<<endl;
        WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD, XYCOORD);
    }
}
else{
    if (now_t > output_int * Psy_time / plot_int)
    {
        cout <<" ==== Writing Data Wait ===="<<endl;
        WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD, XYCOORD);
        output_int++;
    }
}
       Mass_now = CaculateMass(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, U_OLD, XYCOORD);
       Mass_loss = abs(Mass_start - Mass_now);
       cout.precision(6);
       cout << "iter : " << iter << ".   dt :" << dt << ".   time :" << now_t << ".   mass loss :" << Mass_loss << endl;
    }
    return 0;
}