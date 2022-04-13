// Copyright (C) 2022 , National University of Defense Technology
// Xinxin Wang , wxx@nudt.edu.cn

#define limiter(a,b)   ((((a + 1e-7) / (b + 1e-7))>0)?((((a + 1e-7) / (b + 1e-7)) * ((a + 1e-7) / (b + 1e-7)) + ((a + 1e-7) / (b + 1e-7))) * b / (((a + 1e-7) / (b + 1e-7)) * ((a + 1e-7) / (b + 1e-7)) + 1)):0)


