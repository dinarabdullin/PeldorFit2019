#ifndef GENERATION_H
#define GENERATION_H

#include "definitions.h"
#include "chromosome.h"
#include "spin.h"
#include "peldor_calculator.h"
#include <vector>
#include <random>

class Generation
{
public:
	size_t size;
	std::vector<Chromosome> chromosomes;
	std::vector<Chromosome> offspring;

	Generation(genetic_parameters const& genetic_param);

	// Create the first generation
	void create_initial_generation(std::vector<bound> bounds, genetic_parameters const& genetic_param);

	// Calculate the fitness of a single chromosomes
	void score_chromosome(int const& idx, PeldorCalculator const& peldorCalculator, std::vector<experiment> const& experiments, 
		std::vector<Spin> const& spin_system, optimization_parameters const& opt_param, genetic_parameters const& genetic_param);

	// Calculate the fitness of the individual chromosomes
	void score_chromosomes(PeldorCalculator const& peldorCalculator, std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system,
		optimization_parameters const& opt_param, genetic_parameters const& genetic_param);

	// Sort the chromosomes based on their fitness
	void sort_chromosomes();
	
	// Create new chromosomes from the old chromosomes
	void produce_offspring(std::vector<bound> const& bounds, genetic_parameters const& genetic_param);
    
	// Choose two chromosomes based on the tournament_selection
	int tournament_selection(std::mt19937& urng) const;

	// Crossover two chromosomes
	void crossover_chromosomes(Chromosome& chromosome1, Chromosome& chromosome2, double const& prob_crossover, std::mt19937& urng) const;

	// Mutate a chromosome
	void mutate_chromosome(Chromosome& chromosome, double const& prob_mutation, std::vector<bound> const& bounds, std::mt19937& urng) const;
};

#endif
