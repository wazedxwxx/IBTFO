#include <math.h>
#include "Level_Set_function.H"
#include <iostream>
using namespace std;

inline double line(double x, double y, double k, double b, int dir) // dist for y= k x + b
{
    double phi;
    if (y > k * x + b)
    {
        phi = dir * abs(k * x - y + b) / pow(1 + k * k, 0.5);
    }
    else
    {
        phi = -dir * abs(k * x - y + b) / pow(1 + k * k, 0.5);
    }
    return phi;
}

double Level_Set_function(double a, double b)
{
    double phi, phi1, phi2;
    phi1 = line(a, b, pow(3, 0.5) / 3, 0, 1);
    phi2 = line(a, b, pow(3, 0.5) / 3, 0.4, -1);
    //cout << " x " << a << " y " << b << " phi1 " << phi1 << " phi2 " << phi2 << endl;
    phi = ((phi1 < phi2) ? phi1 : phi2);
    return phi;
}
