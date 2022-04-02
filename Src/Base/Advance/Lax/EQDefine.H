#define num_eq 4
#define num_tmp_size num_eq * 12
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

#define Index(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c)
#define Index_U_TMP(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c)
#define Index_U_TMPRK(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + M * N * num_eq
#define Index_F_OLD(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 2 * M * N * num_eq
#define Index_G_OLD(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 3 * M * N * num_eq
#define Index_F_L(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 4 * M * N * num_eq
#define Index_F_R(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 5 * M * N * num_eq
#define Index_G_D(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 6 * M * N * num_eq
#define Index_G_U(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 7 * M * N * num_eq
#define Index_U_L(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 8 * M * N * num_eq
#define Index_U_R(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 9 * M * N * num_eq
#define Index_U_D(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 10 * M * N * num_eq
#define Index_U_U(a, b, c, M, N) ((M) * (b) + (a)) * num_eq + (c) + 11 * M * N * num_eq
