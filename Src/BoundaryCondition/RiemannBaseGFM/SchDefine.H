// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#define num_sch 4
#define M N_x + 2 * num_ghost_cell
#define N N_y + 2 * num_ghost_cell
#define Index_sch(a, b, c) ((M) * (b) + (a)) * num_sch + (c)


