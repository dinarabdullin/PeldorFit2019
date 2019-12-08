#include "generation.h"
#include "tbb/tbb.h"
#include <array>

Generation::Generation(genetic_parameters const& genetic_param)
{
	// Set the size of a generation
	size = genetic_param.size_generation;
	chromosomes.reserve(size);
	offspring.reserve(size);
}

void Generation::create_initial_generation(std::vector<bound> bounds, genetic_parameters const& genetic_param)
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Create chromosomes
	for (size_t i = 0; i < size; ++i) {
		Chromosome* chromosome = new Chromosome(bounds, urng);
		chromosomes.push_back(*chromosome);
		delete chromosome;
	}
}

void Generation::score_chromosome(int const& idx, PeldorCalculator const& peldorCalculator, std::vector<experiment> const& experiments, 
	std::vector<Spin> const& spin_system, optimization_parameters const& opt_param, genetic_parameters const& genetic_param)
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);
	// Simulate the experimental PELDOR time traces and calculate the goodness-of-fit parameter
	double score(0);
	for (size_t j = 0; j < experiments.size(); ++j) {
		// Calculate the PELDOR signal for each chromosome and for each frequency offset
		std::vector<double> simSignal; simSignal.reserve(experiments[j].signalValues.size());
		simSignal = peldorCalculator.calculate_peldor_signal(chromosomes[idx].genes, experiments[j], spin_system, opt_param, urng, j);
		// Calculate the fitness for each chromosome
		if (genetic_param.merit_function == 1) {
			score += peldorCalculator.rmsd(simSignal, experiments[j].signalValues);
		}
		else if (genetic_param.merit_function == 2) {
			score += peldorCalculator.rmsd_over_pearson(simSignal, experiments[j].signalValues);
		}
		else if (genetic_param.merit_function == 3) {
			score += peldorCalculator.pearson(simSignal, experiments[j].signalValues);
		}
		simSignal.clear();
	}
	chromosomes[idx].fitness = score;
}

void Generation::score_chromosomes(PeldorCalculator const& peldorCalculator, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system, 
	optimization_parameters const& opt_param, genetic_parameters const& genetic_param)
{
	size_t const numOfExp = experiments.size();
	// Calculate the fitness of the chromosomes
	tbb::parallel_for(size_t(0), size, [&](size_t i) {
		// Initialize a random generator
		std::random_device rd;
		std::array<unsigned int, std::mt19937::state_size> seed_data;
		std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
		std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
		std::mt19937 urng(seq);
		// Simulate the experimental PELDOR time traces and calculate the goodness-of-fit parameter
		double score(0);
		for (size_t j = 0; j < numOfExp; ++j) {
			// Calculate the PELDOR signal for each chromosome and for each frequency offset
			std::vector<double> simSignal; simSignal.reserve(experiments[j].signalValues.size());
			simSignal = peldorCalculator.calculate_peldor_signal(chromosomes[i].genes, experiments[j], spin_system, opt_param, urng, j);
			// Calculate the fitness for each chromosome
			if (genetic_param.merit_function == 1) {
				score += peldorCalculator.rmsd(simSignal, experiments[j].signalValues);
			}
			else if (genetic_param.merit_function == 2) {
				score += peldorCalculator.rmsd_over_pearson(simSignal, experiments[j].signalValues);
			}
			else if (genetic_param.merit_function == 3) {
				score += peldorCalculator.pearson(simSignal, experiments[j].signalValues);
			}
			simSignal.clear();
		}
		chromosomes[i].fitness = score;
	});
}

void Generation::sort_chromosomes()
{
	std::sort(chromosomes.begin(), chromosomes.end());
}

void Generation::produce_offspring(std::vector<bound> const& bonds, genetic_parameters const& genetic_param)
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Check if the number of chromosomes is even
	size_t pairs(0);
	if (size % 2 == 0) {
		pairs = size / 2;
	}
	else {
		pairs = (size + 1) / 2;
	}

	// Select the pairs of parents and produce an offspring
	for (size_t i = 0; i < pairs; ++i) {
		// Select parents via tournament selection
		Chromosome parent1 = chromosomes[tournament_selection(urng)];
		Chromosome parent2 = chromosomes[tournament_selection(urng)];
		// Crossover parents
		crossover_chromosomes(parent1, parent2, genetic_param.prob_crossover, urng);
		// Mutate parents
		mutate_chromosome(parent1, genetic_param.prob_mutation, bonds, urng);
		mutate_chromosome(parent2, genetic_param.prob_mutation, bonds, urng);
		// Save new chromosomes
		offspring.push_back(parent1);
		if (2 * (i + 1) <= size) {
			offspring.push_back(parent2);
		}
	}
	// Save new generation
	chromosomes = offspring;
	offspring.clear();
}

int Generation::tournament_selection(std::mt19937& urng) const
{
	std::uniform_int_distribution<int> uniform_distr(0, size - 1);
	int n1 = uniform_distr(urng);
	int n2 = uniform_distr(urng);
	if (chromosomes[n1].fitness < chromosomes[n2].fitness) return n1;
	else                                                   return n2;
}

void Generation::crossover_chromosomes(Chromosome& chromosome1, Chromosome& chromosome2, double const& prob_crossover, std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	if (uniform_distr(urng) <= prob_crossover) {
		std::vector<double> genes_temp; genes_temp.reserve(chromosome1.size); // store initial genes of the 1st chromosome
		for (size_t i = 0; i < chromosome1.size; ++i) genes_temp.push_back(chromosome1.genes[i]);
		std::uniform_int_distribution<int> uniform_distr2(1, chromosome1.size - 2);
		int point = uniform_distr2(urng); // choose a random crossover point
		for (int i = 0; i <= point; ++i) {
			chromosome1.genes[i] = chromosome2.genes[i];
			chromosome2.genes[i] = genes_temp[i];
		}
	}
}

void Generation::mutate_chromosome(Chromosome& chromosome, double const& prob_mutation,	std::vector<bound> const& bonds, std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	for (size_t i = 0; i < chromosome.size; ++i) {
		if (uniform_distr(urng) <= prob_mutation) {
			chromosome.genes[i] = chromosome.create_random_gene(bonds[i].lower, bonds[i].upper, urng);
		}
	}
}