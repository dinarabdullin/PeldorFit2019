#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "definitions.h"
#include <vector>
#include <random>

class Chromosome
{
public:
	size_t size;
	std::vector<double> genes;
	double fitness;

	Chromosome(std::vector<bound> const& bounds, std::mt19937& urng);

	// Create a random gene
	double create_random_gene(double const& min_bound, double const& max_bound, std::mt19937& urng) const;

	// Overload less operator
	bool operator < (const Chromosome& A) const { return this->fitness < A.fitness; };
};

#endif