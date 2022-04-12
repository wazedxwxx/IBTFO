#include "ReadParam.H"



struct ReadP ReadPara(char *filename)
{
    struct ReadP prob;
    ParamReader DetectParams;
    Params<double> para(DetectParams.open(filename).numbers());
    prob.max_iter= int(para.get("max_iter", 100000)); // max time steps
    prob.lo_x = para.get("lo_x", 0);                // x direction grid number in computational domain
    prob.lo_y = para.get("lo_y", 0);                // y direction grid number in computational domain

    prob.hi_x = para.get("hi_x", 1); // x direction grid number in computational domain
    prob.hi_y = para.get("hi_y", 1); // y direction grid number in computational domain

    prob.N_x = int(para.get("N_x", 500));             // x direction grid number in computational domain
    prob.N_y = int(para.get("N_y", 200));             // y direction grid number in computational domain
    prob.Psy_time = para.get("Psy_time", 1.0);     // physical time
    prob.CFL_number = para.get("CFL_number", 0.3); // CFL_number
    prob.num_ghost_cell = int(para.get("num_ghost_cell", 2));
    prob.plot_per = int(para.get("plot_per", 100));
    prob.gamma = para.get("gamma", 1.4);  // Gas parameters
    prob.plot_int = para.get("plot_int", 0); // Outputs NUM
    prob.angle = 3.1415926535 * para.get("angle", 30) / 180; // tube angle
    prob.yaxisD = para.get("yaxisD", 0.4);                   // The y-intercept of the shock tube >
    return prob;
}
