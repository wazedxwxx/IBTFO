#include "Slope_limiter.H"
#include <math.h>
using namespace std;
double Slope_limiter(double a,
                     double b)
{
    double r = (a + 1e-7) / (b + 1e-7);
    double c;
    if (r > 0)
    {
        c = (r * r + r) * b / (r * r + 1);
    }
    else
    {
        c = 0;
    }

    return c;
}