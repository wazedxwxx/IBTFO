#include "Slope_limiter.H"
#include <math.h>
using namespace std;
double Slope_limiter(double a,
                     double b)
{
    double r = (a + 1e-7) / (b + 1e-7);
    double c;
    if (r <= 0)
    {
        c = 0;
    }
    else
    {

        c = (r + std::abs(r)) * b / (1.0 + std::abs(r));
    }

    return c;
}