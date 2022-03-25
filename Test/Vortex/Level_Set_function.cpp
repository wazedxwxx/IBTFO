#include <math.h>
#include "Level_Set_function.H"
#include "ParamReader.H"
#include <iostream>
using namespace std;
double Level_Set_function(char*filename,double a, double b)
{
    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    const double f0 = para.get("f0", 0);
    const double f1 = para.get("f1", 1);
    const double rvortex = 0.40;
    const double rl0 = f0 * rvortex;
    const double rl1 = f1 * rvortex;

    const double xlc = 0.50;
    const double ylc = 0.50;
    const double xp = a;
    const double yp = b;
    double r = 0;

    double phi, phi1, phi2;

    r = sqrt((xp - xlc) * (xp - xlc) + (yp - ylc) * (yp - ylc));

    phi1 = rl1 - r;
    phi2 = r - rl0;

    phi = ((phi1 < phi2) ? phi1 : phi2);

    //cout << " x " << a << " y "<<b<<" r "<<r <<endl;
    return phi;
}
