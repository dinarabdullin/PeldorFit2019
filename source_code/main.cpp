#include "definitions.h"
#include "read_config.h"
#include "genetic_algorithm.h"
#include "timer.h"
#include <iostream>
#include <vector>
//#include <memory>

int main(int argc, char** argv)
{  
	// Check that the number of arguments equals 2	
	if ( argc != 2 ) {
		std::cerr << "Incorrect number of arguments!" << std::endl;
        return EXIT_FAILURE;
	}	
	
	// Define input variables
	std::vector<experiment> experiments;
	std::vector<Spin> spin_system;
	optimization_parameters opt_param;
	genetic_parameters genetic_param;
	output_parameters output_param;
	errorplot_parameters errorplot_param;

	// Read a configuration file
	std::cout << std::endl << "Loading input data from the config file..." << std::endl;
	read_config(argv[1], experiments, spin_system, opt_param, genetic_param, output_param, errorplot_param);
	std::cout << "Input data is loaded!" << std::endl;

	// Optimization
	if (errorplot_param.enable == 0) { // Calculate the fit to the PELDOR signals
		double wall_time0(0), wall_time1(0);
		wall_time0 = get_wall_time();
		std::cout << "Optimization of the spin geometry via genetic algorithm... Please be patient!" << std::endl;
		//std::shared_ptr<GeneticAlgorithm> optimizer(new GeneticAlgorithm());
		GeneticAlgorithm* optimizer = new GeneticAlgorithm();
		optimizer->run_optimization(experiments, spin_system, opt_param, genetic_param, output_param);
		wall_time1 = get_wall_time();
		std::cout << "Finished! The optimization took " << (wall_time1 - wall_time0) / 3600 << " hours." << std::endl;
		delete optimizer;
	}
	else if (errorplot_param.enable) { // Calculate the error plot
		double wall_time0(0), wall_time1(0);
		wall_time0 = get_wall_time();
		std::cout << "Recording the error plot... Please be patient!" << std::endl;
		//std::shared_ptr<GeneticAlgorithm> optimizer(new GeneticAlgorithm());
		GeneticAlgorithm* optimizer = new GeneticAlgorithm();
		optimizer->run_error_plot(experiments, spin_system, opt_param, genetic_param, output_param, errorplot_param);
		wall_time1 = get_wall_time();
		std::cout << "Finished! The calculation took " << (wall_time1 - wall_time0) / 3600 << " hours." << std::endl;
		delete optimizer;
	}

	return EXIT_SUCCESS;
}