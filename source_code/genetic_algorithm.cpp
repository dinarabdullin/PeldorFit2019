#include "genetic_algorithm.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
//#include <memory>
#include <tuple>
#include <random>
#include <array>
#include <algorithm>

GeneticAlgorithm::GeneticAlgorithm() {}

void GeneticAlgorithm::run_optimization(std::vector<experiment> const& experiments, std::vector<Spin>& spin_system, 
	optimization_parameters const& opt_param, genetic_parameters const& genetic_param, output_parameters const& output_param) const
{
	// Initialize the calculator of PELDOR signals
	//std::shared_ptr<PeldorCalculator> peldorCalculator(new PeldorCalculator(genetic_param.num_avg));
	PeldorCalculator* peldorCalculator = new PeldorCalculator(genetic_param.num_avg);
	// Calculate the spectrum of the spin system
	if (output_param.record_spectrum) {
		std::cout << "  Recording the spectrum of the spin system... ";
		record_spectrum(*peldorCalculator, experiments[0], spin_system, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Start optimization
	// Pre-compute the excitation probability of spinA by the detection and pump pulses 
	peldorCalculator->exitation_probabilities_spinA(experiments, spin_system[0]);
	// 1st optimization step
	int nopt(1);
	std::cout << "  Optimization step " << nopt << " / " << genetic_param.num_generations_max << std::endl;
	// Create a first generation
	//std::shared_ptr<Generation> generation(new Generation(genetic_param));
	Generation* generation = new Generation(genetic_param);
	generation->create_initial_generation(opt_param.opt_param_bounds, genetic_param);
	// Score the first generation
	generation->score_chromosomes(*peldorCalculator, experiments, spin_system, opt_param, genetic_param);
	// Sort chromosomes according to their fitness
	generation->sort_chromosomes();
	// Save the best fitness
	std::vector<double> best_fitness; best_fitness.reserve(genetic_param.num_generations_max);
	best_fitness.push_back(generation->chromosomes[0].fitness);
	// Proceed with next generations
	for (int nopt = 2; nopt <= genetic_param.num_generations_max; ++nopt) {
		std::cout << "  Optimization step " << nopt << " / " << genetic_param.num_generations_max << std::endl;
		// Create the next generation
		generation->produce_offspring(opt_param.opt_param_bounds, genetic_param);
		// Score the generation
		generation->score_chromosomes(*peldorCalculator, experiments, spin_system, opt_param, genetic_param);
		// Sort chromosomes according to their fitness
		generation->sort_chromosomes();
		// Save the best fitness
		best_fitness.push_back(generation->chromosomes[0].fitness);
	}
	// Save the goodness-of-fit (the fitness of the best chromosome) vs optimization step
	if (output_param.record_score) {
		std::cout << "  Recording the goodness-of-fit vs optimization step... ";
		record_score(best_fitness, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Save the optimized parameters (the genes of the best chromosome)
	if (output_param.record_parameters) {
		std::cout << "  Recording the optimized values of fitting parameters... ";
		record_parameters(generation->chromosomes[0], experiments, opt_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Save the fits to the PELDOR time traces
	if (output_param.record_fit) {
		std::cout << "  Recording the fits to the PELDOR signals... ";
		record_fit(*peldorCalculator, generation->chromosomes[0], experiments, spin_system, opt_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Record the symmetry-related sets of fitting parameters
	if (output_param.record_symmetric_solutions) {
		std::cout << "  Recording the symmetry-related sets of fitting parameters... ";
		record_symmetric_parameters(*peldorCalculator, generation->chromosomes[0], experiments, spin_system, opt_param, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	// Record the error plots
	if (output_param.record_error_plot) {
		std::cout << "  Recording the error plot... ";
		record_error_plot(*peldorCalculator, generation->chromosomes[0], experiments, spin_system, opt_param, genetic_param, output_param);
		std::cout << "Done!" << std::endl;
	}
	delete generation;
	delete peldorCalculator;
}

void GeneticAlgorithm::record_spectrum(PeldorCalculator const& peldorCalculator, experiment const& exp, std::vector<Spin> const& spin_system, output_parameters const& output_param) const
{
	// Calculate the spectrum
	std::vector<std::vector<double>> spectrum; spectrum.reserve(2);
	spectrum = peldorCalculator.calculate_spectrum(exp, spin_system);

	// Create a file to which the spectrum will be saved
	std::ostringstream filename;
	filename << output_param.directory << "spectrum.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);

	// Write the spectrum to the file
	file << std::left << std::setw(20) << "Frequency(GHz)" << std::setw(20) << "Amplitude(a.u.)" << std::endl;
	for (size_t i = 0; i < spectrum[0].size(); i++) {
		file << std::left << std::setprecision(5) << std::setw(20) << spectrum[0][i];
		file << std::left << std::setprecision(5) << std::setw(20) << spectrum[1][i];
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_score(std::vector<double> const& best_fitness, genetic_parameters const& genetic_param, output_parameters const& output_param) const
{
	// Create a file
	std::ostringstream filename;
	filename << output_param.directory << "score.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the column names
	size_t const col_width = 20;
	file << std::left << std::setw(col_width) << "Generation";
	if (genetic_param.merit_function == 1)      file << std::left << std::setw(col_width) << "RMSD" << std::endl;
	else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width) << "RMSD/PCC" << std::endl;
	else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width) << "PCC" << std::endl;
	// Write the data
	for (size_t n = 0; n < best_fitness.size(); ++n) {
		file << std::left << std::setw(col_width) << (n + 1);
		file << std::left << std::setw(col_width) << std::fixed << std::setprecision(5) << best_fitness[n];
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_parameters(Chromosome const& chromosome, std::vector<experiment> const& experiments, optimization_parameters const& opt_param, output_parameters const& output_param) const
{
	// Create a file
	std::ostringstream filename;
	filename << output_param.directory << "parameters.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the column names
	size_t const col_width = 32;
	file << std::left << std::setw(col_width) << "Parameter";
	file << std::left << std::setw(col_width) << "Value";
	file << std::endl;
	// Fill in the rows with the names of fitting parameters and their values
	const size_t numOfExp = experiments.size();
	for (size_t i = 0; i < 27; ++i) {
		file << std::left << std::setw(col_width) << param_names[i];
		if (opt_param.opt_param_numbers[i] != -1) {
			if (std::count(angular_indices, angular_indices + 20, i)) {
				file << std::left << std::setw(col_width) << std::fixed << std::setprecision(0) << chromosome.genes[opt_param.opt_param_numbers[i]] * rad2deg;
			}
			else {
				file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << chromosome.genes[opt_param.opt_param_numbers[i]];
			}
		}
		else {
			if (opt_param.fixed_param_numbers[i] != -1) {
				std::ostringstream row_name;
				if (std::count(angular_indices, angular_indices + 20, i)) {
					row_name << std::fixed << std::setprecision(0) << opt_param.fixed_param_values[opt_param.fixed_param_numbers[i]] * rad2deg << " (fixed)";
					file << std::left << std::setw(col_width) << row_name.str();
				}
				else {
					row_name << std::fixed << std::setprecision(2) << opt_param.fixed_param_values[opt_param.fixed_param_numbers[i]] << " (fixed)";
					file << std::left << std::setw(col_width) << row_name.str();
				}
				row_name.str(std::string()); row_name.clear();
			}
			else {
				file << std::setw(col_width) << " ";
			}	
		}
		file << std::endl;
	}
	// Write the pump efficiency
	if (opt_param.opt_param_numbers[27] != -1) {
		if (opt_param.param_modes[27] == 1) {
			std::ostringstream row_name;
			for (size_t i = 0; i < numOfExp; ++i) { 
				row_name << param_names[27] << " (dataset " << (i + 1) << ")";
				file << std::left << std::setw(col_width) << row_name.str();
				file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << chromosome.genes[opt_param.opt_param_numbers[27] + i];
				file << std::endl;		
				row_name.str(std::string()); row_name.clear();		
			}
		}
		else {
			file << std::left << std::setw(col_width) << param_names[27];
			file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << chromosome.genes[opt_param.opt_param_numbers[27]];
			file << std::endl;
		}	
	}
	else {
		file << std::left << std::setw(col_width) << param_names[27];
		if (opt_param.fixed_param_numbers[27] != -1) {
			file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << opt_param.fixed_param_values[opt_param.fixed_param_numbers[27]];
		}
		else {
			file << std::setw(col_width) << " ";
		}
		file << std::endl;
	}

	file.close();
}

void GeneticAlgorithm::record_fit(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin>& spin_system,
	optimization_parameters const& opt_param, output_parameters const& output_param) const
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Calculate PELDOR signal for all frequency offsets
	const size_t numOfExp = experiments.size();
	std::vector<std::vector<double>> signals; signals.reserve(numOfExp);
	for (size_t j = 0; j < numOfExp; ++j) {
		std::vector<double> signal; signal.reserve(experiments[j].signalValues.size());
		signal = peldorCalculator.calculate_peldor_signal(chromosome.genes, experiments[j], spin_system, opt_param, urng, j);
		signals.push_back(signal);
	}

	// Calculate the maximal number of time/signal points
	size_t Nmax(0);
	for (size_t i = 0; i < numOfExp; ++i) {
		if (experiments[i].timeValues.size() > Nmax) {
			Nmax = experiments[i].timeValues.size();
		}
	}

	// Create a file
	std::ostringstream filename;
	filename << output_param.directory << "fit.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);

	// Write the column names
	size_t const col_width = 10;
	std::ostringstream exp_name;
	std::ostringstream fit_name;
	for (size_t i = 0; i < numOfExp; ++i) {
		exp_name << "exp" << i + 1;
		fit_name << "fit" << i + 1;
		file << std::left << std::setw(col_width) << "t(us)"
			<< std::setw(col_width) << exp_name.str()
			<< std::setw(col_width) << fit_name.str();
		exp_name.str(std::string()); exp_name.clear();
		fit_name.str(std::string()); fit_name.clear();
	}
	file << std::endl;

	// Write the calculated PELDOR signals and their fits
	for (size_t t = 0; t < Nmax; ++t) {
		for (size_t i = 0; i < numOfExp; ++i) {
			if (t < experiments[i].timeValues.size()) {
				file << std::left << std::fixed << std::setprecision(5) << std::setw(col_width) << experiments[i].timeValues[t];
				file << std::left << std::fixed << std::setprecision(5) << std::setw(col_width) << experiments[i].signalValues[t];
				file << std::left << std::fixed << std::setprecision(5) << std::setw(col_width) << signals[i][t];
			}
			else {
				file << std::setw(3 * col_width) << " ";
			}
		}
		file << std::endl;
	}
	file.close();
}

void GeneticAlgorithm::record_symmetric_parameters(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system,
	optimization_parameters const& opt_param, genetic_parameters const& genetic_param, output_parameters const& output_param) const
{
	// Create a vector with the symmetry-related sets of angles
	const size_t n_sym = 16;
	std::vector<std::vector<double>> sym_angles; sym_angles.reserve(n_sym);
	sym_angles = peldorCalculator.calculate_symmetric_angles(chromosome.genes, experiments, spin_system, opt_param);

	// Create a new generation
	std::vector<std::vector<double>> genes_sym; genes_sym.reserve(n_sym);
	for (size_t i = 0; i < n_sym; ++i) {
		genes_sym.push_back(chromosome.genes);
		if (opt_param.opt_param_numbers[2] != -1)  genes_sym[i][opt_param.opt_param_numbers[2]] = sym_angles[i][0];
		if (opt_param.opt_param_numbers[4] != -1)  genes_sym[i][opt_param.opt_param_numbers[4]] = sym_angles[i][1];
		if (opt_param.opt_param_numbers[6] != -1)  genes_sym[i][opt_param.opt_param_numbers[6]] = sym_angles[i][2];
		if (opt_param.opt_param_numbers[8] != -1)  genes_sym[i][opt_param.opt_param_numbers[8]] = sym_angles[i][3];
		if (opt_param.opt_param_numbers[10] != -1) genes_sym[i][opt_param.opt_param_numbers[10]] = sym_angles[i][4];
		if (opt_param.opt_param_numbers[14] != -1) genes_sym[i][opt_param.opt_param_numbers[14]] = sym_angles[i][5];
		if (opt_param.opt_param_numbers[16] != -1) genes_sym[i][opt_param.opt_param_numbers[16]] = sym_angles[i][6];
		if (opt_param.opt_param_numbers[18] != -1) genes_sym[i][opt_param.opt_param_numbers[18]] = sym_angles[i][7];
		if (opt_param.opt_param_numbers[20] != -1) genes_sym[i][opt_param.opt_param_numbers[20]] = sym_angles[i][8];
		if (opt_param.opt_param_numbers[22] != -1) genes_sym[i][opt_param.opt_param_numbers[22]] = sym_angles[i][9];
	}
	genetic_parameters genetic_param_sym(genetic_param);
	genetic_param_sym.size_generation = n_sym;
	//std::shared_ptr<Generation> generation_sym(new Generation(genetic_param_sym));
	Generation* generation_sym = new Generation(genetic_param_sym);
	for (size_t i = 0; i < n_sym; ++i) {
		generation_sym->chromosomes.push_back(chromosome);
		generation_sym->chromosomes[i].genes = genes_sym[i];
	}

	// Modify the values of fixed angles to their symmetry-related counterparts
	std::vector<optimization_parameters> opt_param_sym; opt_param_sym.reserve(16);
	for (size_t i = 0; i < n_sym; ++i) {
		opt_param_sym.push_back(opt_param);
		if (opt_param.fixed_param_numbers[2]  != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[2]]  = sym_angles[i][0];
		if (opt_param.fixed_param_numbers[4]  != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[4]]  = sym_angles[i][1];
		if (opt_param.fixed_param_numbers[6]  != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[6]]  = sym_angles[i][2];
		if (opt_param.fixed_param_numbers[8]  != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[8]]  = sym_angles[i][3];
		if (opt_param.fixed_param_numbers[10] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[10]] = sym_angles[i][4];
		if (opt_param.fixed_param_numbers[14] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[14]] = sym_angles[i][5];
		if (opt_param.fixed_param_numbers[16] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[16]] = sym_angles[i][6];
		if (opt_param.fixed_param_numbers[18] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[18]] = sym_angles[i][7];
		if (opt_param.fixed_param_numbers[20] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[20]] = sym_angles[i][8];
		if (opt_param.fixed_param_numbers[22] != -1) opt_param_sym[i].fixed_param_values[opt_param.fixed_param_numbers[22]] = sym_angles[i][9];
	}
	
	// Score the chromosomes
	for (size_t i = 0; i < n_sym; ++i) {
		generation_sym->score_chromosome(i, peldorCalculator, experiments, spin_system, opt_param_sym[i], genetic_param_sym);
	}

	// Create a file
	std::ostringstream filename;
	filename << output_param.directory << "symmetric_parameters.dat";
	std::fstream file;
	file.open(filename.str().c_str(), std::ios::out);
	// Write the column names
	size_t const col_width1 = 20;
	size_t const col_width2 = 12;
	file << std::left << std::setw(col_width1) << "Parameter";
	for (size_t i = 0; i < n_sym; ++i) {
		file << std::left << std::setw(col_width2) << transform_names[i];
	}
	file << std::endl;
	// Fill in the columns with the names of angular parameters and their values
	for (size_t i = 0; i < 10; ++i) {
		file << std::left << std::setw(col_width1) << angle_names[i];
		for (size_t k = 0; k < n_sym; ++k) {
			file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(0) << sym_angles[k][i] * rad2deg;
		}
		file << std::endl;
	}
    // Write the goodness-of-fit values
	if (genetic_param.merit_function == 1)      file << std::left << std::setw(col_width1) << "RMSD";
	else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width1) << "RMSD/PCC";
	else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width1) << "PCC";
	for (size_t k = 0; k < n_sym; ++k) {
		file << std::left << std::setw(col_width2) << std::fixed << std::setprecision(4) << generation_sym->chromosomes[k].fitness;
	}
	file.close();
	delete generation_sym;
}

void GeneticAlgorithm::record_error_plot(PeldorCalculator const& peldorCalculator, Chromosome const& chromosome, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system,
	optimization_parameters const& opt_param, genetic_parameters const& genetic_param, output_parameters const& output_param) const
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);
	
	// Calculate error plots
	const size_t Nplots = output_param.error_plot_variables.size();
	const size_t Nsamples = output_param.error_plot_size;
	size_t Nvar(0);
	for (size_t i = 0; i < Nplots; ++i) {
		// Number of variables in the error plot
		Nvar = output_param.error_plot_variables[i].size();
		// Check if the selected variables of the error plot were optimized
		int idxVar(0);
		bool badVar(false);
		for (size_t j = 0; j < Nvar; ++j) {
			idxVar = output_param.error_plot_variables[i][j] - 1;
			if (opt_param.opt_param_numbers[idxVar] == -1) {
				std::cout << "  Inappropriate variable (was not optimized): " << output_param.error_plot_variables[i][j] << std::endl;
				badVar = true;
				break;
			}
		}
		if (badVar) continue;
		// Create a generation with different values of the error plot variables
		genetic_parameters genetic_param_erp(genetic_param);
		genetic_param_erp.size_generation = output_param.error_plot_size;
		//std::shared_ptr<Generation> generation_erp(new Generation(genetic_param_erp));
		Generation* generation_erp = new Generation(genetic_param_erp);
		// Initialize all chromosomes with the optimized chromosome
		for (size_t k = 0; k < Nsamples; ++k) {
			generation_erp->chromosomes.push_back(chromosome); 
		}
		// Vary the genes corresponding to the error plot variables
		int idxGene(0);
		for (size_t j = 0; j < Nvar; ++j) {
			idxGene = opt_param.opt_param_numbers[output_param.error_plot_variables[i][j] - 1];
			for (size_t k = 0; k < Nsamples; ++k) {
				generation_erp->chromosomes[k].genes[idxGene] = chromosome.create_random_gene(opt_param.opt_param_bounds[idxGene].lower, opt_param.opt_param_bounds[idxGene].upper, urng);
			}
		}	
		// Score the chromosomes of the generation
		generation_erp->score_chromosomes(peldorCalculator, experiments, spin_system, opt_param, genetic_param_erp);
		// Sort chromosomes based on their fitness
		generation_erp->sort_chromosomes();
		// Create a file
		std::ostringstream filename;
		filename << output_param.directory << "error_plot";
		for (size_t j = 0; j < Nvar; ++j) filename << "_" << output_param.error_plot_variables[i][j];
		filename << ".dat";
		std::fstream file;
		file.open(filename.str().c_str(), std::ios::out);
		// Write the column names
		size_t const col_width = 32;
		std::ostringstream col_name;
		for (size_t j = 0; j < Nvar; ++j) {
			col_name << param_names[output_param.error_plot_variables[i][j] - 1];
			file << std::left << std::setw(col_width) << col_name.str();
			col_name.str(std::string()); 
			col_name.clear();
		}
		if      (genetic_param.merit_function == 1) file << std::left << std::setw(col_width) << "RMSD";
		else if (genetic_param.merit_function == 2) file << std::left << std::setw(col_width) << "RMSD/PCC";
		else if (genetic_param.merit_function == 3) file << std::left << std::setw(col_width) << "PCC";
		file << std::endl;
		// Write the error plot
		for (size_t k = 0; k < Nsamples; ++k) {
			for (size_t j = 0; j < Nvar; ++j) {
				idxVar = output_param.error_plot_variables[i][j] - 1;
				idxGene = opt_param.opt_param_numbers[output_param.error_plot_variables[i][j] - 1];
				if (std::count(angular_indices, angular_indices + 20, idxVar)) {
					file << std::left << std::setw(col_width) << std::fixed << std::setprecision(0) << generation_erp->chromosomes[k].genes[idxGene] * rad2deg;
				}
				else {
					file << std::left << std::setw(col_width) << std::fixed << std::setprecision(2) << generation_erp->chromosomes[k].genes[idxGene];
				}
			}
			file << std::left << std::setw(col_width) << std::fixed << std::setprecision(5) << generation_erp->chromosomes[k].fitness;
			file << std::endl;
		}
		file.close();
		delete generation_erp;
	}
}

void GeneticAlgorithm::run_error_plot(std::vector<experiment> const& experiments, std::vector<Spin>& spin_system, optimization_parameters const& opt_param,
	genetic_parameters const& genetic_param, output_parameters& output_param, errorplot_parameters const& errorplot_param) const
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Create a chromosome and set its genes to the previously optimised values
	//std::shared_ptr<Chromosome> chromosome_erp(new Chromosome(opt_param.opt_param_bounds, urng));
	Chromosome* chromosome_erp = new Chromosome(opt_param.opt_param_bounds, urng);
	chromosome_erp->genes = errorplot_param.optimized_genes;
	
	// Initialize the calculator of PELDOR signals
	//std::shared_ptr<PeldorCalculator> peldorCalculator(new PeldorCalculator(genetic_param.num_avg));
	PeldorCalculator* peldorCalculator = new PeldorCalculator(genetic_param.num_avg);
	// Pre-compute the excitation probability of spinA by the detection and pump pulses 
	peldorCalculator->exitation_probabilities_spinA(experiments, spin_system[0]);

	// Set the output parameters
	output_param.directory = errorplot_param.output_directory;
	output_param.error_plot_variables = errorplot_param.error_plot_variables;
	output_param.error_plot_size = errorplot_param.error_plot_size;

	// Record the error plot
	record_error_plot(*peldorCalculator, *chromosome_erp, experiments, spin_system, opt_param, genetic_param, output_param);
	delete peldorCalculator;
	delete chromosome_erp;
}