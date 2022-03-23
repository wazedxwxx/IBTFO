#include "Slope_limiter.H"
#include <math.h>
using namespace std;
double Slope_limiter(double a,
                     double b)
{

    double c;
    if (a * b <= 0)
    {
        c  = 0;
    }
    else
    {
        double r = a  / b ;

        if (r > 1)
        {
            c = b;
        }
        else
        {
            c = r * b;
        }
    }

    return c;
}