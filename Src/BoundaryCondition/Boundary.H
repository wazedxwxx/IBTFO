void Boundary(const int N_x,
              const int N_y,
              const int num_ghost_cell,
              const double gamma,
              double *U_OLD,
              double *U_NEW,
              double *XYCOORD,
              double *SCHEME_IDX,
              const int ndevices,
              const int device);

void Global_Boundary(const int N_x,
                     const int N_y,
                     const int num_ghost_cell,
                     const double gamma,
                     double *U_OLD,
                     double *U_NEW,
                     const int ndevices,
                     const int device);
