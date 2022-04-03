#define num_sch 8
#define M N_x + 2 * num_ghost_cell
#define N N_y + 2 * num_ghost_cell
#define Index_sch(a, b, c) ((N) * (b) + (a)) * num_sch + (c)