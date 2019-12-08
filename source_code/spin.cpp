#include "spin.h"
#include <cmath>

Spin::Spin() 
{}

void Spin::initialize()
{
	if (n.empty()) {
		Nstates = 1;
		Ncomp = 1;
		Icomp.reserve(Ncomp);
		Icomp.push_back(1);
	}
	else if (n.size() == 1) {
		Nstates = static_cast<size_t>(pow((2 * I[0] + 1), n[0]));
		Ncomp = static_cast<size_t>(2 * I[0] * n[0] + 1);
		Icomp.reserve(Ncomp);
		int I_index = static_cast<int>(2 * I[0] - 1);
		int n_index = n[0] - 1;
		int intensity;
		for (size_t i = 0; i < Ncomp; i++) {
			intensity = esr_weights[I_index][n_index][i];
			Icomp.push_back(intensity);
		} 
	}
	else if (n.size() == 2) {
		Nstates = static_cast<size_t>( pow((2 * I[0] + 1), n[0]) * pow((2 * I[1] + 1), n[1]) );
		size_t Ncomp1 = static_cast<size_t>(2 * I[0] * n[0] + 1);
		size_t Ncomp2 = static_cast<size_t>(2 *I[1] * n[1] + 1 );		
		Ncomp = Ncomp1 * Ncomp2;
		Icomp.reserve(Ncomp);
		int I_index1 = static_cast<int>(2 * I[0] - 1);
		int n_index1 = n[0] - 1;
		int I_index2 = static_cast<int>(2 * I[1] - 1);
		int n_index2 = n[1] - 1;
		int intensity;
		for (size_t i = 0; i < Ncomp1; i++) {
			for (size_t j = 0; j < Ncomp2; j++) {
				intensity = esr_weights[I_index1][n_index1][i] * esr_weights[I_index2][n_index2][j];
				Icomp.push_back(intensity);
			}
		}
	}
}

double Spin::calculate_gfactor(std::vector<double> const& field_direction, std::mt19937& urng) const
{
	double gValue(0);
	gValue = sqrt(g[0] * g[0] * field_direction[0] * field_direction[0] +
		          g[1] * g[1] * field_direction[1] * field_direction[1] +
				  g[2] * g[2] * field_direction[2] * field_direction[2]);
	if (gStrain) {
		double dg(0);
		dg = g[0] * gStrain[0] * field_direction[0] * field_direction[0] +
			 g[1] * gStrain[1] * field_direction[1] * field_direction[1] +
			 g[2] * gStrain[2] * field_direction[2] * field_direction[2];
		dg /= gValue;
		std::normal_distribution<double> normal_distr(0.0, 1.0);
		gValue += fwhm2sd * dg * normal_distr(urng);
	}
	return gValue;
}

double Spin::calculate_hfc(std::vector<double> const& field_direction, std::mt19937& urng) const
{
	double AValue(0);
	AValue = 1e-3 * sqrt(A[0] * A[0] * field_direction[0] * field_direction[0] +
		                 A[1] * A[1] * field_direction[1] * field_direction[1] +
		                 A[2] * A[2] * field_direction[2] * field_direction[2]);
	if (AStrain) {
		double dA(0);
		dA = 1e-6 * (A[0] * AStrain[0] * field_direction[0] * field_direction[0] +
			         A[1] * AStrain[1] * field_direction[1] * field_direction[1] +
			         A[2] * AStrain[2] * field_direction[2] * field_direction[2]);
		dA /= AValue;
		std::normal_distribution<double> normal_distr(0.0, 1.0);
		AValue += fwhm2sd * dA * normal_distr(urng);
	}
	return AValue;
}

double Spin::calculate_hfc2(std::vector<double> const& field_direction) const
{
	double AValue(0);
	AValue = 1e-3 * sqrt(A[3] * A[3] * field_direction[0] * field_direction[0] +
		                 A[4] * A[4] * field_direction[1] * field_direction[1] +
		                 A[5] * A[5] * field_direction[2] * field_direction[2]);
	return AValue;
}

double Spin::calculate_broadering(std::mt19937& urng) const
{
	std::normal_distribution<double> normal_distr(0.0, 1.0);
	return 1e-3 * 0.5 * lwpp * normal_distr(urng);
}

std::vector<double> Spin::calculate_resfreq(std::vector<double> const& field_direction, experiment const& exp, std::mt19937& urng) const
{
	std::vector<double> resfreq; resfreq.reserve(Ncomp);
	double f(0);
	// Calculate the g-factor
	double gValue = calculate_gfactor(field_direction, urng);
	// For no nuclei
	if (n.empty()) {
		// Calculate the resonance frequency increment due to the inhomogenious broadering
		double df = calculate_broadering(urng);
		// Calculate the resonance frequency 
		f = F * gValue * exp.magnField + df;
		resfreq.push_back(f);
	}
	// For the case of one sort of identical nuclei
	if (n.size() == 1) {
		// Calculate the hyperfine coupling constant
		double AValue = calculate_hfc(field_direction, urng);
		// Nuclear magnetic quantum number
		double m = -I[0] * n[0];
		double df(0);
		for (size_t i = 0; i < Ncomp; ++i) {
			// Calculate the resonance frequency increment due to the inhomogenious broadering
			df = calculate_broadering(urng);
			// Calculate the resonance frequency 
			f = F * gValue * exp.magnField + AValue * m + df;
			resfreq.push_back(f);
			++m;
		}
	}
	// For the case of two sorts of identical nuclei
	if (n.size() == 2) {
		// Calculate the hyperfine coupling constants
		double AValue1 = calculate_hfc(field_direction, urng);
		double AValue2 = calculate_hfc2(field_direction);
		// Number of spectral components for each type of nucleus
		size_t Ncomp1 = static_cast<size_t>(2 * I[0] * n[0] + 1);
		size_t Ncomp2 = static_cast<size_t>(2 * I[1] * n[1] + 1);
		// Nuclear magnetic quantum numbers
		double m1 = -I[0] * n[0];
		double m2 = -I[1] * n[1];
		double df(0);
		for (size_t i = 0; i < Ncomp1; ++i) {
			m2 = -I[1] * n[1];
			for (size_t j = 0; j < Ncomp2; ++j) {
				// Calculate the resonance frequency increment due to the inhomogenious broadering
				df = calculate_broadering(urng);
				f = F * gValue * exp.magnField + AValue1 * m1 + AValue2 * m2 + df;
				resfreq.push_back(f);
				++m2;
			}
			++m1;
		}
	}
	return resfreq;
}