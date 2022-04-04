#define num_sch 2
#define M N_x + 2 * num_ghost_cell
#define N N_y + 2 * num_ghost_cell
#define Index_sch(a, b, c) ((M) * (b) + (a)) * num_sch + (c)