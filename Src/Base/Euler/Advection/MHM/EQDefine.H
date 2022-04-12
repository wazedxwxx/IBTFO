#define num_eq 4
#define num_tmp_hight 7
#define num_tmp_size num_eq * num_tmp_hight
#define M N_x + 2 * num_ghost_cell
#define N N_y + 2 * num_ghost_cell
#define H_SIZE ((N_y + 2 * num_ghost_cell) / ndevices)
#define LOWER (max(H_SIZE * device - 4 * num_ghost_cell, 0))
#define UPPER (min(H_SIZE * (device + 1) + 4 * num_ghost_cell, N))

#define REAL_LOWER (max(H_SIZE * device, 0))
#define REAL_UPPER (min(H_SIZE * (device + 1), N))


/*
U_TMP
U_TMPRK
F_OLD
G_OLD
F_L
F_R
G_D
G_U
U_L
U_R
U_D
U_U
*/

#define Index(a, b, c) ((M) * (b) + (a)) * num_eq + (c)

#define Index_U_TMP(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c)
#define Index_U_TMPRK(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + num_eq
#define Index_F_OLD(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 2 * num_eq
#define Index_G_OLD(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 2 * num_eq
#define Index_F_L(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 3 * num_eq
#define Index_F_R(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 4 * num_eq
#define Index_G_D(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 3 * num_eq
#define Index_G_U(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 4 * num_eq
#define Index_U_L(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 5 * num_eq
#define Index_U_R(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 6 * num_eq
#define Index_U_D(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 5 * num_eq
#define Index_U_U(a, b, c) ((M) * (b) + (a)) * num_tmp_size + (c) + 6 * num_eq