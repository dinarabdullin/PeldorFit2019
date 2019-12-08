#ifndef GENETIC_ALGORITH_H
#define GENETIC_ALGORITH_H

#include "definitions.h"
#include "generation.h"
#include <vector>

class GeneticAlgorithm
{
public:
	GeneticAlgorithm();

	// Run optimization using the genetic algorithm
	void run_optimization(std::vector<experiment> const& experiments, std::vector<Spin>& spin_system, optimization_parameters const& opt_param,
		genetic_parameters const& genetic_param, output_parameters const& output_param) const;

	// Save the spectrum of the spin system
	void record_spectrum(PeldorCalculator const& peldorCalculator, experiment const& exp, std::vector<Spin> const& spin_system, output_parameters const& output_param) const;

	// Save the goodness-of-fit (the fitness of the best chromosome) vs optimization step
	void record_score(std::vector<double> const& best_fitness, genetic_parameters const& genetic_param, output_parameters const& output_param) const;

	// Save the optimized parameters (the genes of the best chromosome)
	void record_parameters(Chromosome const& chromosome, std::vector<experiment> const& experiments, optimization_parameters const& opt_param, output_parameters const& output_param) const;

	// Save the fits to the PELDOR time traces
	void record_fit(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin>& spin_system,
		optimization_parameters const& opt_param, output_parameters const& output_param) const;

	// Save the symmetry-related sets of fitting parameters
	void record_symmetric_parameters(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system,
		optimization_parameters const& opt_param, genetic_parameters const& genetic_param, output_parameters const& output_param) const;

	// Record the error plot
	void record_error_plot(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system,
		optimization_parameters const& opt_param, genetic_parameters const& genetic_param, output_parameters const& output_param) const;

	// Record the error plot using previous fitting results
	void run_error_plot(std::vector<experiment> const& experiments, std::vector<Spin>& spin_system, optimization_parameters const& opt_param,
		genetic_parameters const& genetic_param, output_parameters& output_param, errorplot_parameters const& errorplot_param) const;
};

#endif