void ComputeDt(const double Psy_L,
                 const double Psy_H,
                 const int N_x,
                 const int N_y,
                 const int num_ghost_cell,
                 const double gamma,
                 const double CFL_number,
                 double *U_OLD,
                 double *dt);