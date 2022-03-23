#include "Slope_limiter.H"
#include <math.h>
using namespace std;
double Slope_limiter(double a,
                     double b)
{

    double c;
    if (a * b > 0)
    {
        double r = a / b;

        if (r > 1)
        {
            c = b;
        }
        else
        {
            c = r * b;
        }
    }
    else
    {
        c = 0;
    }

    return c;
}