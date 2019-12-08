#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <string>

const double PI = 3.141592653589793238462;
const double deg2rad = PI / 180.0;
const double rad2deg = 180.0 / PI;
const double F = 13.99624; // betta/h  [GHz/T]
const double fwhm2sd = 1.0 / 2.35482; // FWHM = 2.35482 * SD

const int angular_indices[20] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

const int meanval_indices[13] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 25};

char const * const param_names[] = {
	"distance mean 1 (nm)",
	"distance width 1 (nm)",
	"xi mean 1 (deg)",
	"xi width 1 (deg)",
	"phi mean 1 (deg)",
	"phi width 1 (deg)",
	"alpha mean 1 (deg)",
	"alpha width 1 (deg)",
	"beta mean 1 (deg)",
	"beta width 1 (deg)",
	"gamma mean 1 (deg)",
	"gamma width 1 (deg)",
	"distance mean 2 (nm)",
	"distance width 2 (nm)",
	"xi mean 2 (deg)",
	"xi width 2 (deg)",
	"phi mean 2 (deg)",
	"phi width 2 (deg)",
	"alpha mean 2 (deg)",
	"alpha width 2 (deg)",
	"beta mean 2 (deg)",
	"beta width 2 (deg)",
	"gamma mean 2 (deg)",
	"gamma width 2 (deg)",
	"ratio 1 to 2",
	"J mean (MHz)",
	"J width (MHz)",
	"pump efficiency" 
};

char const * const angle_names[] = {
	"xi mean 1 (deg)",
	"phi mean 1 (deg)",
	"alpha mean 1 (deg)",
	"beta mean 1 (deg)",
	"gamma mean 1 (deg)",
	"xi mean 2 (deg)",
	"phi mean 2 (deg)",
	"alpha mean 2 (deg)",
	"beta mean 2 (deg)",
	"gamma mean 2 (deg)",
};

char const * const transform_names[] = {
	"A:I_B:I",
	"A:I_B:Rx",
	"A:I_B:Ry",
	"A:I_B:Rz",
	"A:Rx_B:I",
	"A:Rx_B:Rx",
	"A:Rx_B:Ry",
	"A:Rx_B:Rz",
	"A:Ry_B:I",
	"A:Ry_B:Rx",
	"A:Ry_B:Ry",
	"A:Ry_B:Rz",
	"A:Rz_B:I",
	"A:Rz_B:Rx",
	"A:Rz_B:Ry",
	"A:Rz_B:Rz"
};

const int esr_weights[3][6][19] = {
{
	{ 1, 1 },
	{ 1, 2, 1 },
	{ 1, 3, 3, 1 },
	{ 1, 4, 6, 4, 1 },
	{ 1, 5, 10, 10, 5, 1 },
	{ 1, 6, 15, 20, 15, 6, 1 }
},
{
	{ 1, 1, 1 },
	{ 1, 2, 3, 2, 1 },
	{ 1, 3, 6, 7, 6, 3, 1 },
	{ 1, 4, 10, 16, 19, 16, 10, 4, 1 },
	{ 1, 5, 15, 30, 45, 51, 45, 30, 15, 5, 1 },
	{ 1, 6, 21, 50, 90, 126, 141, 126, 90, 50, 21, 6, 1 }
},
{
	{ 1, 1, 1, 1 },
	{ 1, 2, 3, 4, 3, 2, 1 },
	{ 1, 3, 6, 10, 12, 12, 10, 6, 3, 1 },
	{ 1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1 },
	{ 1, 5, 15, 35, 65, 101, 135, 155, 155, 135, 101, 65, 35, 15, 5, 1 },
	{ 1, 6, 21, 56, 120, 216, 336, 456, 546, 580, 546, 456, 336, 216, 120, 56, 21, 6, 1 }
} 
};

struct experiment
{
	std::vector<double> timeValues;
	std::vector<double> signalValues;
	double detPiLength;
	double detPiHalfLength;
	double pumpPiLength;
	double detFreq;
	double pumpFreq;
	double magnField;
};

struct bound
{
	double lower;
	double upper;
};

struct optimization_parameters
{
	std::vector<int> opt_param_numbers;
	std::vector<int> fixed_param_numbers;
	std::vector<int> param_modes;
	std::vector<bound> opt_param_bounds;
	std::vector<double> fixed_param_values;
	int multiSpin;
};

struct genetic_parameters
{
	int num_generations_max;
	int size_generation;
	double prob_crossover;
    double prob_mutation;
	int merit_function;
	int num_avg;
};

struct output_parameters
{
	std::string directory;
	int record_spectrum;
    int record_score;
    int record_parameters;
    int record_fit;
    int record_symmetric_solutions;
    int record_error_plot;
    std::vector<std::vector<int>> error_plot_variables;
	int error_plot_size;
};

struct errorplot_parameters
{
	int enable;
	std::vector<std::vector<int>> error_plot_variables;
	int error_plot_size;
	std::string input_directory;
	std::string output_directory;
	std::vector<double> optimized_genes;
};
 
#endif