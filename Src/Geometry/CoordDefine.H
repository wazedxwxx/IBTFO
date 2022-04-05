#define num_coord 6
#define M N_x + 2 * num_ghost_cell
#define N N_y + 2 * num_ghost_cell
#define Index_Coord(a, b, c) ((M) * (b) + (a)) * 6 + (c)

