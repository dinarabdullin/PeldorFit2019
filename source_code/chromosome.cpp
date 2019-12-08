#include "chromosome.h"

Chromosome::Chromosome(std::vector<bound> const& bounds, std::mt19937& urng)
{
	// Set the size of a chromosome
	size = bounds.size();
	genes.reserve(size);
	// Create initial genes
	double gene(0);
	for (size_t i = 0; i < size; ++i) {
		gene = create_random_gene(bounds[i].lower, bounds[i].upper, urng);
		genes.push_back(gene);
	}
}

double Chromosome::create_random_gene(double const& min_bound, double const& max_bound, std::mt19937& urng) const
{
	std::uniform_real_distribution<> uniform_distr(min_bound, max_bound);
	return uniform_distr(urng);
}