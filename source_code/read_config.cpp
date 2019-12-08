#include "read_config.h"
#ifdef _WIN32
#pragma comment(lib, "libconfig")
#endif
#include "libconfig.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

void read_config(const char* path_config, 
	             std::vector<experiment>& experiments,
	             std::vector<Spin>& spin_system,
                 optimization_parameters& opt_param,
				 genetic_parameters& genetic_param,
				 output_parameters& output_param,
				 errorplot_parameters& errorplot_param)
{
	// Read a configuration file
	config_t cfg, *cf;
	cf = &cfg;
	config_init(cf);
	if (!config_read_file(cf, path_config)) {
		std::cerr << config_error_file(cf) << ": " << config_error_line(cf) << " - " << config_error_text(cf) << std::endl;
		config_destroy(cf);
		exit(EXIT_FAILURE);
	}
	// Enable the auto convertion of data types
	config_set_auto_convert(cf, true);
	////////////////////////////////////////////////////////////////////////////
	//                        EXPERIMENTAL PARAMETERS                         //
	////////////////////////////////////////////////////////////////////////////
	const config_setting_t *experimentals = config_lookup(cf, "experimentals");
	size_t numOfSignals = config_setting_length(experimentals);
	experiments.reserve(numOfSignals);
	std::cout << "  The number of PELDOR time traces is " << numOfSignals << std::endl;
	const config_setting_t *set;
	const char *path = NULL;
	std::string dataline;
	double timeValue, signalValue;
	std::vector<double> timeValues, signalValues;
	experiment single_experiment;
	double detPiLength, detPiHalfLength, pumpPiLength, detFreq, pumpFreq, magnField;
	
	for (size_t i = 0; i < numOfSignals; ++i) {
		set = config_setting_get_elem(experimentals, i);
		path = config_setting_get_string(config_setting_get_member(set, "filename"));
		std::ifstream datafile(path, std::ios::in);
		if (!datafile.is_open()) {
			std::cerr << " Error: Could not open the data file: " << std::endl << path << std::endl;
			exit(EXIT_FAILURE);
		}
		else {
			while (getline(datafile, dataline)) {
				std::stringstream datastream(dataline);
				datastream >> timeValue >> signalValue;
				timeValues.push_back(timeValue);
				signalValues.push_back(signalValue);
			}
			single_experiment.timeValues = timeValues;
			single_experiment.signalValues = signalValues;
			timeValues.clear();
			signalValues.clear();
		}
		detPiLength     = config_setting_get_float(config_setting_get_member(set, "detPiLength"));
		detPiHalfLength = config_setting_get_float(config_setting_get_member(set, "detPiHalfLength"));
		pumpPiLength    = config_setting_get_float(config_setting_get_member(set, "pumpPiLength"));
		detFreq         = config_setting_get_float(config_setting_get_member(set, "detFreq"));
		pumpFreq        = config_setting_get_float(config_setting_get_member(set, "pumpFreq"));
		magnField       = config_setting_get_float(config_setting_get_member(set, "magnField"));
		single_experiment.detPiHalfLength = detPiHalfLength;
		single_experiment.detPiLength     = detPiLength;
		single_experiment.pumpPiLength    = pumpPiLength;
		single_experiment.detFreq         = detFreq;
		single_experiment.pumpFreq        = pumpFreq;
		single_experiment.magnField        = magnField;
		experiments.push_back(single_experiment);
	}
	////////////////////////////////////////////////////////////////////////////
	//								 SPIN PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	// The number of spins
	int nSpins(2);
	try {
		nSpins = config_setting_get_int(config_lookup(cf, "nSpins"));
	} catch (const std::exception& e) {}
	if ((nSpins != 2) && (nSpins != 3)) {
		std::cerr << "Incorrect number of spins!" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Check if the multi-spin mode is swiched on
	if (nSpins == 3) opt_param.multiSpin = config_setting_get_int(config_lookup(cf, "multiSpinEffects"));
	else             opt_param.multiSpin = 0;

	// Read the spectroscopic parameters of spin A
	Spin spinA;
	const config_setting_t *pspinA = config_lookup(cf, "spinA");
	for (size_t i = 0; i < 3; ++i) {
		spinA.g[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "g"), i);
	}
	size_t sizeA = config_setting_length(config_setting_get_member(pspinA, "n"));
	for (size_t i = 0; i < sizeA; ++i) {
		spinA.n.push_back(config_setting_get_int_elem(config_setting_get_member(pspinA, "n"), i));
		spinA.I.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "I"), i));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 0));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 1));
		spinA.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinA, "A"), 3 * i + 2));
	}
	if (config_setting_length(config_setting_get_member(pspinA, "gStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i) 
			spinA.gStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "gStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinA.gStrain[i] = 0.0;
	}
	if (config_setting_length(config_setting_get_member(pspinA, "AStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i) 
			spinA.AStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinA, "AStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinA.AStrain[i] = 0.0;
	}
	spinA.lwpp = config_setting_get_float(config_setting_get_member(pspinA, "lwpp"));
	spinA.initialize();
	spin_system.push_back(spinA);

	// Read the spectroscopic parameters of spin B
	Spin spinB;
	const config_setting_t *pspinB = config_lookup(cf, "spinB");
	for (size_t i = 0; i < 3; ++i) {
		spinB.g[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "g"), i);
	}
	size_t sizeB = config_setting_length(config_setting_get_member(pspinB, "n"));
	for (size_t i = 0; i < sizeB; ++i) {
		spinB.n.push_back(config_setting_get_int_elem(config_setting_get_member(pspinB, "n"), i));
		spinB.I.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "I"), i));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 0));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 1));
		spinB.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinB, "A"), 3 * i + 2));
	}
	if (config_setting_length(config_setting_get_member(pspinB, "gStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i)
			spinB.gStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "gStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinB.gStrain[i] = 0.0;
	}
	if (config_setting_length(config_setting_get_member(pspinB, "AStrain")) != 0) {
		for (size_t i = 0; i < 3; ++i)
			spinB.AStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinB, "AStrain"), i);
	}
	else {
		for (size_t i = 0; i < 3; ++i)
			spinB.AStrain[i] = 0.0;
	}
	spinB.lwpp = config_setting_get_float(config_setting_get_member(pspinB, "lwpp"));
	spinB.initialize();
	spin_system.push_back(spinB);

	// Read the spectroscopic parameters of spin C
	if (nSpins == 3) {
		Spin spinC;
		const config_setting_t *pspinC = config_lookup(cf, "spinC");
		for (size_t i = 0; i < 3; ++i) {
			spinC.g[i] = config_setting_get_float_elem(config_setting_get_member(pspinC, "g"), i);
		}
		size_t sizeC = config_setting_length(config_setting_get_member(pspinC, "n"));
		for (size_t i = 0; i < sizeC; ++i) {
			spinC.n.push_back(config_setting_get_int_elem(config_setting_get_member(pspinC, "n"), i));
			spinC.I.push_back(config_setting_get_float_elem(config_setting_get_member(pspinC, "I"), i));
			spinC.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinC, "A"), 3 * i + 0));
			spinC.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinC, "A"), 3 * i + 1));
			spinC.A.push_back(config_setting_get_float_elem(config_setting_get_member(pspinC, "A"), 3 * i + 2));
		}
		if (config_setting_length(config_setting_get_member(pspinC, "gStrain")) != 0) {
			for (size_t i = 0; i < 3; ++i)
				spinC.gStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinC, "gStrain"), i);
		}
		else {
			for (size_t i = 0; i < 3; ++i)
				spinC.gStrain[i] = 0.0;
		}
		if (config_setting_length(config_setting_get_member(pspinC, "AStrain")) != 0) {
			for (size_t i = 0; i < 3; ++i)
				spinC.AStrain[i] = config_setting_get_float_elem(config_setting_get_member(pspinC, "AStrain"), i);
		}
		else {
			for (size_t i = 0; i < 3; ++i)
				spinC.AStrain[i] = 0.0;
		}
		spinC.lwpp = config_setting_get_float(config_setting_get_member(pspinC, "lwpp"));
		spinC.initialize();
		spin_system.push_back(spinC);
		spin_system.reserve(3);
	}
	else {
		spin_system.reserve(2);
	}

	////////////////////////////////////////////////////////////////////////////
	//							  FITTING PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	opt_param.opt_param_numbers.reserve(28);
	opt_param.fixed_param_numbers.reserve(28);
	opt_param.param_modes.reserve(28);
	const config_setting_t *fit_param = config_lookup(cf, "parameters");
	const config_setting_t *single_param;
	int opt(0);
	int mode(0);
	bound single_bound;
	double value(0);
	int count(0);
	int count2(0);
	for (int i = 0; i < 28; ++i) {
		// Read parameter after parameter
		single_param = config_setting_get_elem(fit_param, i);
		// Check if this parameter is enabled
		opt = config_setting_get_int(config_setting_get_member(single_param, "opt"));
		if (opt == 0) {
			opt_param.opt_param_numbers.push_back(-1);
			opt_param.fixed_param_numbers.push_back(-1);
			opt_param.param_modes.push_back(0);
		}
		else if (opt == 1) {
			opt_param.fixed_param_numbers.push_back(-1);
			opt_param.opt_param_numbers.push_back(count);
			++count;
			mode = config_setting_get_int(config_setting_get_member(single_param, "mode"));
			opt_param.param_modes.push_back(mode);
			single_bound.lower = config_setting_get_float_elem(config_setting_get_member(single_param, "range"), 0);
			single_bound.upper = config_setting_get_float_elem(config_setting_get_member(single_param, "range"), 1);
			if (std::count(angular_indices, angular_indices + 20, i)) {
				single_bound.lower *= deg2rad;
				single_bound.upper *= deg2rad;
			}
			opt_param.opt_param_bounds.push_back(single_bound);
			if ((i == 27) && (mode == 1)) {
				for (size_t j = 0; j < (numOfSignals - 1); ++j) {
					opt_param.opt_param_bounds.push_back(single_bound);
					++count;
				}
			}
		}
		else if (opt == 2) {
			opt_param.opt_param_numbers.push_back(-1);
			opt_param.fixed_param_numbers.push_back(count2);
			++count2;
			mode = config_setting_get_int(config_setting_get_member(single_param, "mode"));
			opt_param.param_modes.push_back(mode);
			value = config_setting_get_float(config_setting_get_member(single_param, "value"));
			if (std::count(angular_indices, angular_indices + 20, i)) {
				value *= deg2rad;
			}
			opt_param.fixed_param_values.push_back(value);
		}
		else {
			std::cerr << "Incorrect value of opt!" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	//opt_param.opt_param_bounds.reserve(count);
	//opt_param.fixed_param_values.reserve(count2);

	////////////////////////////////////////////////////////////////////////////
	//			   	   PARAMETERS OF THE GENETIC ALGORITHM                    //
	////////////////////////////////////////////////////////////////////////////
	config_setting_t *genetic = config_lookup(cf, "genetic");
	genetic_param.num_generations_max = config_setting_get_int(config_setting_get_member(genetic, "num_generations_max"));
	genetic_param.size_generation = config_setting_get_int(config_setting_get_member(genetic, "size_generation"));
	genetic_param.prob_crossover = config_setting_get_float(config_setting_get_member(genetic, "prob_crossover"));
    genetic_param.prob_mutation = config_setting_get_float(config_setting_get_member(genetic, "prob_mutation"));
	genetic_param.merit_function = config_setting_get_int(config_setting_get_member(genetic, "merit_function"));
	genetic_param.num_avg = config_setting_get_int(config_setting_get_member(genetic, "num_averages"));
	////////////////////////////////////////////////////////////////////////////
	//				               OUTPUT PARAMETERS                          //
	////////////////////////////////////////////////////////////////////////////
	config_setting_t *output = config_lookup(cf, "output");
	const char *pfileroot = NULL;
	pfileroot = config_setting_get_string(config_setting_get_member(output, "directory"));
	std::string output_dir(pfileroot);
	output_param.directory = output_dir;
    output_param.record_spectrum = config_setting_get_int(config_setting_get_member(output, "record_spectrum"));
    output_param.record_score = config_setting_get_int(config_setting_get_member(output, "record_score"));
    output_param.record_parameters = config_setting_get_int(config_setting_get_member(output, "record_parameters"));
    output_param.record_fit = config_setting_get_int(config_setting_get_member(output, "record_fit"));
    output_param.record_symmetric_solutions = config_setting_get_int(config_setting_get_member(output, "record_symmetric_solutions"));
    output_param.record_error_plot = config_setting_get_int(config_setting_get_member(output, "record_error_plot"));
    config_setting_t *ep_var = config_setting_get_member(output, "error_plot_variables");
	int nErrorPlots = config_setting_length(ep_var);
	int nVarParam(0);
	if (nErrorPlots == 0) {
		output_param.record_error_plot = 0;
	}
	else {
		output_param.error_plot_variables.reserve(nErrorPlots);
		for (int i = 0; i < nErrorPlots; ++i) {
			nVarParam = config_setting_length(config_setting_get_elem(ep_var, i));
			std::vector<int> varParam;  varParam.reserve(nVarParam);
			for (int j = 0; j < nVarParam; ++j) varParam.push_back(config_setting_get_int_elem(config_setting_get_elem(ep_var,i),j));
			output_param.error_plot_variables.push_back(varParam);
			varParam.clear();
		}
	}
	output_param.error_plot_size = config_setting_get_int(config_setting_get_member(output, "error_plot_size"));
	////////////////////////////////////////////////////////////////////////////
	//				            'ERROR PLOT ONLY' MODE                        //
	////////////////////////////////////////////////////////////////////////////
	const config_setting_t *errorplot = config_lookup(cf, "error_plot_only");
	errorplot_param.enable = config_setting_get_int(config_setting_get_member(errorplot, "enable"));
	if (errorplot_param.enable) {
		config_setting_t *ep_var = config_setting_get_member(errorplot, "error_plot_variables");
		int nErrorPlots = config_setting_length(ep_var);
		int nVarParam(0);
		if (nErrorPlots != 0) {
			errorplot_param.error_plot_variables.reserve(nErrorPlots);
			for (int i = 0; i < nErrorPlots; ++i) {
				nVarParam = config_setting_length(config_setting_get_elem(ep_var, i));
				std::vector<int> varParam;  varParam.reserve(nVarParam);
				for (int j = 0; j < nVarParam; ++j) {
					varParam.push_back(config_setting_get_int_elem(config_setting_get_elem(ep_var,i),j));
				}
				errorplot_param.error_plot_variables.push_back(varParam);
				varParam.clear();
			}
		}
		errorplot_param.error_plot_size = config_setting_get_int(config_setting_get_member(errorplot, "error_plot_size"));

		// Read out the file with the optimized parameters
		errorplot_param.optimized_genes.reserve(opt_param.opt_param_bounds.size());
		// Set the number of lines to read
		int n_datalines(28);
		if ((opt_param.opt_param_numbers[27] != -1) && (opt_param.param_modes[27] == 1)) n_datalines += (numOfSignals - 1);
		// Set the current number of line (-1 line corresponds to the capitel)
		int count(-1); 
		// Read the file line by line
		double datavalue;
		std::string dataline;
		const char *pInputFile = NULL;
		pInputFile = config_setting_get_string(config_setting_get_member(errorplot, "input_directory"));
		std::ifstream inputFile(pInputFile, std::ios::in);
		if (!inputFile.is_open()) {
			std::cerr << " Error: Could not open the data file!" << std::endl << pInputFile << std::endl;
			exit(EXIT_FAILURE);
		}
		else {
			while (getline(inputFile, dataline)) {
				if ((count >= 0) && (count < n_datalines) && (opt_param.opt_param_numbers[count] != -1)) {
					std::string subline = dataline.substr(32);
					std::stringstream datastream(subline);
					datastream >> datavalue;
					if (std::count(angular_indices, angular_indices + 20, count)) {
						errorplot_param.optimized_genes.push_back(datavalue * deg2rad);
					}
					else {
						errorplot_param.optimized_genes.push_back(datavalue);
					}
				}
				++count;
			}
		}
		const char *pOutputFile = NULL;
		pOutputFile = config_setting_get_string(config_setting_get_member(errorplot, "output_directory"));
		std::string outputFile(pOutputFile);
		errorplot_param.output_directory = outputFile;
	}
	// Close a configuration file
	config_destroy(cf);
};