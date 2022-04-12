void Conserve2Flux(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *F,
                   double *G,
                   const int ndevices,
                   const int device);

void Conserve2FluxX(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *G,
                   const int ndevices,
                   const int device);

void Conserve2FluxY(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *F,
                   const int ndevices,
                   const int device);

void Conserve2FluxZ(const int N_x,
                   const int N_y,
                   const int num_ghost_cell,
                   const double gamma,
                   double *U,
                   double *H,
                   const int ndevices,
                   const int device);
