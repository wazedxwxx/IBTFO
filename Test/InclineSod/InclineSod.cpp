#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <math.h>
#include <algorithm>

#include "openacc.h"
#include "EQDefine.H"
#include "CoordDefine.H"
#include "SchDefine.H"
#include "Initialize.H"
#include "Psy_coord.H"
#include "WriteData.H"
#include "Write_LS.H"
#include "Write_SCH.H"
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

int main(int argc, char **argv)
{

    ParamReader DetectParams;
    // parameter settings model
    char *filename = argv[1];
    Params<double> para(DetectParams.open(filename).numbers());
    const int max_iter = int(para.get("max_iter", 100000)); // max time steps
    const double lo_x = para.get("lo_x", 0);                // x direction grid number in computational domain
    const double lo_y = para.get("lo_y", 0);                // y direction grid number in computational domain

    const double hi_x = para.get("hi_x", 1); // x direction grid number in computational domain
    const double hi_y = para.get("hi_y", 1); // y direction grid number in computational domain

    const int N_x = int(para.get("N_x", 500));             // x direction grid number in computational domain
    const int N_y = int(para.get("N_y", 200));             // y direction grid number in computational domain
    const double Psy_time = para.get("Psy_time", 1.0);     // physical time
    const double CFL_number = para.get("CFL_number", 0.3); // CFL_number
    const int num_ghost_cell = int(para.get("num_ghost_cell", 2));
    const int plot_per = int(para.get("plot_per", 100));
    const double gamma = para.get("gamma", 1.4);  // Gas parameters
    const int plot_int = para.get("plot_int", 0); // Outputs NUM
    int output_int = 1;

    cout << " ==== parameters are read ====" << endl;

    const double h = 1.0 / N_x;      // cell size
    const double delta_t = 1 / 1000; // Preset time step
    double now_t = 0;                // Preset time
    int iter = 0;
    double dt = 0.0;
    const double Psy_L = hi_x - lo_x;
    const double Psy_H = hi_y - lo_y;

    int ndevices = acc_get_num_devices(acc_get_device_type());
    // double Mass_start, Mass_now, Mass_loss;

    double *U_OLD = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq];
    double *XYCOORD = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_coord];
    double *U_TMP = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_tmp_size];
    double *U_NEW = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_eq];
    double *SCHEME_IDX = new double[(N_x + 2 * num_ghost_cell) * (N_y + 2 * num_ghost_cell) * num_sch];
    double TMP_DT[ndevices];
    for (int device = 0; device < ndevices; device++)
    {
        acc_set_device_num(device, acc_device_default);
#pragma acc enter data create(U_OLD [(M)*LOWER * num_eq:(M) * (UPPER - LOWER) * num_eq],             \
                              U_NEW [(M)*LOWER * num_eq:(M) * (UPPER - LOWER) * num_eq],             \
                              U_TMP [(M)*LOWER * num_tmp_size:(M) * (UPPER - LOWER) * num_tmp_size], \
                              XYCOORD [(M)*LOWER * num_coord:(M) * (UPPER - LOWER) * num_coord],     \
                              SCHEME_IDX [(M)*LOWER * num_sch:(M) * (UPPER - LOWER) * num_sch], TMP_DT [device:1])
    }

    cout << " ====  memory allocation complete ====" << endl;

    acc_wait_all();
    for (int device = 0; device < ndevices; device++)
    {
        acc_set_device_num(device, acc_device_default);
        printf("Launching device %d\n", device);
#pragma acc update device(U_OLD [(M)*LOWER * num_eq:(M) * (UPPER - LOWER) * num_eq],             \
                          U_NEW [(M)*LOWER * num_eq:(M) * (UPPER - LOWER) * num_eq],             \
                          U_TMP [(M)*LOWER * num_tmp_size:(M) * (UPPER - LOWER) * num_tmp_size], \
                          XYCOORD [(M)*LOWER * num_coord:(M) * (UPPER - LOWER) * num_coord],     \
                          SCHEME_IDX [(M)*LOWER * num_sch:(M) * (UPPER - LOWER) * num_sch]) async

        Psy_coord(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, XYCOORD, ndevices, device);

        Level_Set(filename, N_x, N_y, num_ghost_cell, XYCOORD, ndevices, device);

        Scheme_Index(N_x, N_y, num_ghost_cell, XYCOORD, SCHEME_IDX, ndevices, device);

        Initialize(filename, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, ndevices, device);

        Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, SCHEME_IDX, ndevices, device);

        // WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD,XYCOORD);

        // Write_LS(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, XYCOORD);
        // WriteIDX2TXT(N_x, N_y, num_ghost_cell, SCHEME_IDX);
        // Mass_start = CaculateMass(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, U_OLD, XYCOORD);
        // cout << " ====  Initial Mass is " <<Mass_start <<" ===="<< endl;

#pragma acc update self(U_OLD [(M)*REAL_LOWER * num_eq:(M) * (REAL_UPPER - REAL_LOWER) * num_eq],             \
                        U_NEW [(M)*REAL_LOWER * num_eq:(M) * (REAL_UPPER - REAL_LOWER) * num_eq],             \
                        U_TMP [(M)*REAL_LOWER * num_tmp_size:(M) * (REAL_UPPER - REAL_LOWER) * num_tmp_size], \
                        XYCOORD [(M)*REAL_LOWER * num_coord:(M) * (REAL_UPPER - REAL_LOWER) * num_coord],     \
                        SCHEME_IDX [(M)*REAL_LOWER * num_sch:(M) * (REAL_UPPER - REAL_LOWER) * num_sch]) async
    }
    acc_wait_all();

    for (int device = 0; device < ndevices; device++)
    {
        acc_set_device_num(device, acc_device_default);

        ComputeDt(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, CFL_number, U_OLD, XYCOORD, TMP_DT, ndevices, device);

        if (device + 2 > ndevices)
            dt = *min_element(TMP_DT, TMP_DT + ndevices);

        //        cout << "Device " << device << " REAL_LOWER " << REAL_LOWER << " REAL_UPPER " << REAL_UPPER << " LOWER " << LOWER << " UPPER " << UPPER << endl;
    }

    while (now_t < Psy_time && iter < max_iter)
    {
        acc_wait_all();
        for (int device = 0; device < ndevices; device++)
        {
            acc_set_device_num(device, acc_device_default);
#pragma acc update device(U_OLD [(M)*LOWER * num_eq:(M) * (2 * num_ghost_cell) * num_eq],                          \
                          U_OLD [(M) * (UPPER - 2 * num_ghost_cell) * num_eq:(M) * (2 * num_ghost_cell) * num_eq], \
                          U_NEW [(M)*LOWER * num_eq:(M) * (2 * num_ghost_cell) * num_eq],                          \
                          U_NEW [(M) * (UPPER - 2 * num_ghost_cell) * num_eq:(M) * (2 * num_ghost_cell) * num_eq]) async

            if (plot_int == 0)
            {
                if (iter % plot_per == 0)
                {
                    cout << " ==== Writing Data Wait ====" << endl;
                    WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD, XYCOORD, ndevices, device);
                }
            }
            else
            {
                if (now_t > output_int * Psy_time / plot_int || iter == 0)
                {
                    cout << " ==== Writing Data Wait ====" << endl;
                    WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD, XYCOORD, ndevices, device);
                    if (device == ndevices - 1)
                        output_int = output_int + 1 * (iter > 0);
                }
            }

            TimeAdvance(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, dt, U_OLD, U_TMP, U_NEW, XYCOORD, SCHEME_IDX, ndevices, device);

            Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, XYCOORD, SCHEME_IDX, ndevices, device);

            Global_Boundary(N_x, N_y, num_ghost_cell, gamma, U_OLD, U_NEW, ndevices, device);

            ComputeDt(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, gamma, CFL_number, U_OLD, XYCOORD, TMP_DT, ndevices, device);

            if (device + 2 > ndevices)
                dt = *min_element(TMP_DT, TMP_DT + ndevices);

#pragma acc update self(U_OLD [(M)*REAL_LOWER * num_eq:(M) * (2 * num_ghost_cell) * num_eq],                          \
                        U_OLD [(M) * (REAL_UPPER - 2 * num_ghost_cell) * num_eq:(M) * (2 * num_ghost_cell) * num_eq], \
                        U_NEW [(M)*REAL_LOWER * num_eq:(M) * (2 * num_ghost_cell) * num_eq],                          \
                        U_NEW [(M) * (REAL_UPPER - 2 * num_ghost_cell) * num_eq:(M) * (2 * num_ghost_cell) * num_eq]) async
        }
        acc_wait_all();

        now_t = now_t + dt;
        iter++;

        /* Mass_now = CaculateMass(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, U_OLD, XYCOORD);
       Mass_loss = abs(Mass_start - Mass_now); */

        cout.precision(6);
        cout << "iter : " << iter << ".   dt :" << dt << ".   time :" << now_t << endl;
        // cout << "iter : " << iter << ".   dt :" << dt << ".   time :" << now_t << ".   mass loss :" << Mass_loss << endl;
    }

    for (int device = 0; device < ndevices; device++)
    {
        acc_set_device_num(device, acc_device_default);
        WriteData(lo_x, lo_y, Psy_L, Psy_H, N_x, N_y, num_ghost_cell, iter, now_t, gamma, U_OLD, XYCOORD, ndevices, device);

#pragma acc exit data copyout(U_OLD [(M)*REAL_LOWER * num_eq:(M) * (REAL_UPPER - REAL_LOWER) * num_eq],             \
                              U_NEW [(M)*REAL_LOWER * num_eq:(M) * (REAL_UPPER - REAL_LOWER) * num_eq],             \
                              U_TMP [(M)*REAL_LOWER * num_tmp_size:(M) * (REAL_UPPER - REAL_LOWER) * num_tmp_size], \
                              XYCOORD [(M)*REAL_LOWER * num_coord:(M) * (REAL_UPPER - REAL_LOWER) * num_coord],     \
                              SCHEME_IDX [(M)*REAL_LOWER * num_sch:(M) * (REAL_UPPER - REAL_LOWER) * num_sch])
    }
    // Write_SCH(Psy_L, Psy_H, N_x, N_y, num_ghost_cell, SCHEME_IDX);
    return 0;
}