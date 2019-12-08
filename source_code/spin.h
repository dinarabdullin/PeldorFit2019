#ifndef SPIN_H
#define SPIN_H

#include "definitions.h"
#include <vector>
#include <random>

class Spin
{
public:
	double g[3];
	std::vector<int> n;
	std::vector<double> I;
	std::vector<double> A;
	double gStrain[3];
	double AStrain[3];
	double lwpp;
	size_t Nstates; // Number of microstates
	size_t Ncomp; // Number of spectral components
	std::vector<int> Icomp; // Intensities of individual spectral components

	Spin();

	// Calculate the number and intensities of ESR lines
	void initialize();

	// Calculate the g-factor
	double calculate_gfactor(std::vector<double> const& field_direction, std::mt19937& urng) const;

	// Calculate the hyperfine coupling constant
	double calculate_hfc(std::vector<double> const& field_direction, std::mt19937& urng) const;

	double calculate_hfc2(std::vector<double> const& field_direction) const;

	// Calculate the resonance frequency increment due to the inhomogenious broadering
	double calculate_broadering(std::mt19937& urng) const;

	// Calculate the resonance frequencies
	std::vector<double> calculate_resfreq(std::vector<double> const& field_direction, experiment const& exp, std::mt19937& urng) const;
};

#endif