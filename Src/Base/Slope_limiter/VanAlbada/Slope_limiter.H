#define limiter(a,b)   ((((a + 1e-7) / (b + 1e-7))>0)?((r + std::abs(((a + 1e-7) / (b + 1e-7)))) * b / (1.0 + std::abs(((a + 1e-7) / (b + 1e-7))))):0)


