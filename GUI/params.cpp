#include <qstring.h>
#include "params.h"

Params::Params()
{
    PARAM_SET params[] = {

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

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

{"DELTA_T", 120.0, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[sec]"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NT_ANIMATION", 10, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps)."},

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
