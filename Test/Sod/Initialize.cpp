#include "Initialize.H"

void Initialize(char *filename,
                const double Psy_L,
                const double Psy_H,
                const int N_x,
                const int N_y,
                const int num_ghost_cell,
                const double gamma,
                double *U_OLD,
                double *U_NEW,
                double *XYCOORD,
                const int ndevices,
                const int device)
{
    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    const double rho_L = para.get("rho_L", 1);     // left side density
    const double u_L = para.get("u_L", 0);         // left side x-vel
    const double v_L = para.get("v_L", 0);         // left side y-vel
    const double p_L = para.get("p_L", 1);         // left side pressure
    const double rho_R = para.get("rho_R", 0.125); // right side density
    const double u_R = para.get("u_R", 0);         // right side x-vel
    const double v_R = para.get("v_R", 0);         // right side y-vel
    const double p_R = para.get("p_R", 0.1);       // right side pressure
    
    int lower = LOWER;
    int upper = UPPER;
    acc_set_device_num(device, acc_device_default);
#pragma acc data present(U_OLD [(M) * lower * num_eq:(M) * (upper - lower) * num_eq], \
                         U_NEW [(M) * lower * num_eq:(M) * (upper - lower) * num_eq], \
                         XYCOORD [(M) * lower * num_coord:(M) * (upper - lower) * num_coord])
    {

#pragma acc parallel loop async
        for (int j = lower; j < upper; j++)
        {
#pragma acc loop 
            for (int i = 0; i < N_x + 2 * num_ghost_cell; i++)
            {
                if (j < (N_x + 2 * num_ghost_cell) / 2)
                //   if (XYCOORD[Index_Coord(i, j, 1, N_x + 2 * num_ghost_cell)] <
                //       -pow(3, 0.5) * XYCOORD[Index_Coord(i, j, 0, N_x + 2 * num_ghost_cell)] + 0.5 * pow(3, 0.5) + 0.5)
                {
                    U_OLD[Index(i, j, 0)] = rho_L;
                    U_OLD[Index(i, j, 1)] = rho_L * u_L;
                    U_OLD[Index(i, j, 2)] = rho_L * v_L;
                    U_OLD[Index(i, j, 3)] = p_L / (gamma - 1) + 0.5 * rho_L * (u_L * u_L + v_L * v_L);
                }
                else
                {
                    U_OLD[Index(i, j, 0)] = rho_R;
                    U_OLD[Index(i, j, 1)] = rho_R * u_R;
                    U_OLD[Index(i, j, 2)] = rho_R * v_R;
                    U_OLD[Index(i, j, 3)] = p_R / (gamma - 1) + 0.5 * rho_R * (u_R * u_R + v_R * v_R);
                }
                for (int k = 0; k < num_eq; k++)
                {
                    U_NEW[Index(i, j, k)] = U_OLD[Index(i, j, k)];
                }
            }
        }
    }
}
