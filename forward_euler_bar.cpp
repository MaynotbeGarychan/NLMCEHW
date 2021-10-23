#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

struct Init_param
{
    float k_bar;
    float df;
    float Fmax;
    int num_iter;
};

struct Step_output
{
    float Kt;
    float ui;
    float Fi;
};

float formula_Fbar(float k_bar, float u)
{
    float F_bar;
    F_bar = (u + 1.5*pow(u,2) + 05*pow(u,3)) + k_bar;
    return F_bar;

}

float formula_Kt(float k_bar, float u)
{
    float Kt;
    Kt = 1+3*u+1.5*pow(u,2)+k_bar;
    return Kt;
}


int main()
{
    Init_param init_param;
    init_param.k_bar = 0.8;
    init_param.Fmax = -2;
    init_param.df = -0.05;
    init_param.num_iter = init_param.Fmax/init_param.df;

    int size = init_param.num_iter + 1;
    std::vector<Step_output> step_output[size];
    
    // give first value
    Step_output step_output0;
    step_output0.Kt = 0;
    step_output0.ui = 0;
    step_output0.Fi = 0;
    step_output->push_back(step_output0);
    
    // loop 
    for (int i = 1; i < init_param.num_iter; i++)
    {
        Step_output temp;
        float du;

        // Calculate the value
        temp.Kt = formula_Kt(init_param.k_bar,step_output->back().ui);
        du = init_param.df/temp.Kt;

        // Update the displacement and force
        temp.ui = step_output->back().ui + du;
        temp.Fi = step_output->back().Fi + init_param.df;
        
        step_output->push_back(temp);
    }

    return 0;
}
