#ifndef PELDOR_CALCULATOR_H
#define PELDOR_CALCULATOR_H

#include "definitions.h"
#include "spin.h"
#include <vector>
#include <tuple>

class PeldorCalculator
{
private:
	int n_fields;

public:
	double threshold;
	std::vector<std::vector<double>> fieldDirA;
	std::vector<double> gValuesA;
	std::vector<std::vector<double>> detProbsA;
	std::vector<std::vector<double>> pumpProbsA;

	PeldorCalculator(int const& num_avg);

	// Simulate the random orientation of a magnetic field
	void initialize_field();

	// Calulate the minimal and maximal resonance frequncies of the spin system for a given field
	std::vector<double> minmax_resfreq(experiment const& exp, std::vector<Spin> const& spin_system, std::mt19937& urng) const;

	// Calculate the spectrum of the spin system
	std::vector<std::vector<double>> calculate_spectrum(experiment const& exp, std::vector<Spin> const& spin_system) const;

	// Calculate the excitation probabilities for the spin A
	void exitation_probabilities_spinA(std::vector<experiment> const& experiments, Spin const& spinA);

	// Calculate the probablity of excitation of a spin by the detection pulses
	double excitation_probability_detection(std::vector<double> const& freqs, experiment const& exp, Spin const& spin) const;

	// Calculate the probablity of excitation of a spin by the pump pulse
	double excitation_probability_pump(std::vector<double> const& freqs, experiment const& exp, Spin const& spin) const;

	// Set the value of the optimization parameter
	double get_parameter_value(std::vector<double> const& opt_param_values, optimization_parameters const& opt_param, int const& nOffset, int const& idx) const;
	
	// Calculate the value of an optimization parameter based on the given type of distribution
	double get_value_from_distribution(std::vector<double> const& opt_param_values, optimization_parameters const& opt_param, std::mt19937& urng,
		int const& idx_mean, int const& idx_width) const;

	// Calculate the direction of a distance vector from the angles
	std::vector<double> direction_from_angles(double const& xi, double const& phi) const;

	// Calculate the length and direction of a distance vector from other two distance vectors
	std::tuple<double, std::vector<double>> direction_from_vectors(double const& dist1, double const& dist2,
		std::vector<double> const& dir1, std::vector<double> const& dir2) const;

	double modulation_depth_correction(double const& w) const;

	// Simulate a PELDOR signal
	std::vector<double> calculate_peldor_signal(std::vector<double> const& opt_param_values, experiment const& exp, std::vector<Spin> const& spin_system,
		optimization_parameters const& opt_param, std::mt19937& urng, int const& nOffset) const;

    // Root-mean-square deviation
	double rmsd(std::vector<double> const& x, std::vector<double> const& y) const;

	// Root-mean-square deviation / Pearson correlation coefficient
	double rmsd_over_pearson(std::vector<double> const& x, std::vector<double> const& y) const;

	// Pearson correlation coefficient
	double pearson(std::vector<double> const& x, std::vector<double> const& y) const;

	// Calculate xi and phi angles from a direction vector
	std::vector<double> angles_from_direction(std::vector<double> const& v) const;

	// Calculate Euler angles alpha,beta, gamma from a rotation matrix
	std::vector<double> angles_from_rotation_matrix(std::vector<std::vector<double>> const& RM) const;

	// Calculate symmetry-related sets of angular parameters
	std::vector<std::vector<double>> symmetric_angles(double const& xi, double const& phi, double const& alpha, double const& beta, double const& gamma) const;

	// Calculate symmetry-related sets of fitting parameters
	std::vector<std::vector<double>> calculate_symmetric_angles(std::vector<double> const& opt_param_values, std::vector<experiment> const& experiments, 
		std::vector<Spin> const& spin_system, optimization_parameters const& opt_param) const;
};

#endif