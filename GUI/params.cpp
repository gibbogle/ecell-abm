#include <qstring.h>
#include "params.h"
#include "log.h"

LOG_USE();

Params::Params()
{
    PARAM_SET params[] = {

{"BLOCK_SIZE", 3, 0, 0,
"Dimension of initial block",
"Dimension N of initial block of cells (NxNxN)"},

{"GROWTH", 1, 0, 1,
"Cells grow?",
"Simulate cell growth and division."},

{"DIVIDE_TIME_MEDIAN", 24, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"HOURS", 100, 0, 0,
"Number of hours",
"Length of the simulation.\n\
[hours]"},

{"DELTA_T", 100.0, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[sec]"},

{"SOLVER",0,0,0,
"Solver",
"Choice of solvers: Explicit, Backward Euler, Runge-Kutta 45"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 6, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps)."},

{"F_DRAG", 20, 0, 0,
"Drag factor",
"Velocity v = F/Fdrag, angular velocity w = M/Mdrag\n\
    F = total force vector, M = total moment vector, Fdrag = force drag coefficient, Mdrag = moment drag coefficient."},

{"M_DRAG", 200, 0, 0,
"Moment drag factor",
"Moment drag factor."},

{"F_ALPHA", 0.5, 0, 1,
"Smoothing factor",
"To control oscillations, cell motion is determined by force F and moment M computed as a weighted sum of current contacts and previous F and M: \n\
 F = (1-Falpha)*(sum of forces) + Falpha*Fprev, where Falpha is the force smoothing factor and Fprev is the value of F in the previous time step \n\
 M = (1-Malpha)*(sum of moments) + Malpha*Mprev, where Malpha is the moment smoothing factor and Mprev is the value of M in the previous time step."},

{"M_ALPHA", 0.5, 0, 1,
"Moment smoothing factor",
"Moment smoothing factor."},

{"F_JIGGLE", 0.01, 0, 0,
"Jiggle factor",
"A random perturbation can be applied to cell motion, independently to linear and angular velocities:\n\
 Each component of velocity is perturbed by dv = R*Fjiggle/Fdrag\n\
 and each component of angular velocity is perturbed by dw = R*Mjiggle/Mdrag \n\
 where a value of Gaussian(0,1) distributed random variate R is generated for each velocity component."},

{"M_JIGGLE", 0.01, 0, 0,
"Moment jiggle factor",
"Moment jiggle factor."},

{"CELLML_FILE", 0, 0, 0,
"",
"CellML file for the cell growth model"}

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}
}

PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
    workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}
