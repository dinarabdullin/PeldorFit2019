#include "peldor_calculator.h"
#include "definitions.h"
#include "rotations.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <array>
#include <algorithm>
#include <memory>
#include <cmath>
#include <algorithm>

PeldorCalculator::PeldorCalculator(int const& num_avg) 
{
	threshold = 1e-3;
	n_fields = num_avg;
	initialize_field();
}

void PeldorCalculator::initialize_field()
{
	// Initialize random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 urng(seq);
	std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
	
	// Calculate random orientations of the magnetic field vector
	fieldDirA.reserve(n_fields);
	double fxi(0), fphi(0);
	std::vector<double> singleDir;
	singleDir.reserve(3);
	for (int f = 0; f < n_fields; ++f) {
		fphi = 0.0 + (2.0*PI - 0.0) * uniform_distr(urng);
		fxi = acos(uniform_distr(urng));
		if (uniform_distr(urng) < 0.5) fxi = PI - fxi;
		singleDir.push_back(sin(fxi)*cos(fphi));
		singleDir.push_back(sin(fxi)*sin(fphi));
		singleDir.push_back(cos(fxi));
		fieldDirA.push_back(singleDir);
		singleDir.clear();
	}
}

std::vector<double> PeldorCalculator::minmax_resfreq(experiment const& exp, std::vector<Spin> const& spin_system, std::mt19937& urng) const
{
	double fmin, fmax;
	double fminF, fmaxF;
	std::vector<double> resfreq;
	for (size_t i = 0; i < spin_system.size(); ++i) {
		for (int f = 0; f < n_fields; ++f) {
			// Compute the resonance frequences for the single orientation of the magnetic field
			resfreq = spin_system[i].calculate_resfreq(fieldDirA[f], exp, urng);
			// Determine the minimal and maximal resonance frequences for the single orientation of the magnetic field
			fminF = *std::min_element(resfreq.begin(), resfreq.end());
			fmaxF = *std::max_element(resfreq.begin(), resfreq.end());
			// Determine the minimal and maximal resonance frequences for all magntic field orientations
			if ((i == 0) && (f == 0)) {
				fmin = fminF;
				fmax = fmaxF;
			}
			else {
				if (fminF < fmin) fmin = fminF;
				if (fmaxF > fmax) fmax = fmaxF;
			}
			resfreq.clear();
		}
	}
	std::vector<double> fminmax; fminmax.reserve(2);
	fminmax.push_back(fmin);
	fminmax.push_back(fmax);
	return fminmax;
}

std::vector<std::vector<double>> PeldorCalculator::calculate_spectrum(experiment const& exp, std::vector<Spin> const& spin_system) const
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Determine the maximal and minimal frequencies
	double fmin, fmax;
	std::vector<double> fminmax; fminmax.reserve(2);
	fminmax = minmax_resfreq(exp, spin_system, urng);
	fmin = fminmax[0];
	fmax = fminmax[1];
	// Increase slightly the frequency range
	fmin -= 0.100;
	fmax += 0.100;

	// Set the frequency increment to 1 MHz
	double df = 0.001;
	// Round the max and min frequencies to MHz
	fmin = ceil(fmin / df) * df;
	fmax = ceil(fmax / df) * df;
	// Calculate the number of increments
	int Nf = static_cast<int>((fmax - fmin) / df);
	// Correlate frequencies with the indices
	int idx0 = static_cast<int>(fmin / df) + 1;
	int idx;
	// The frequency values
	std::vector<double> fValues; fValues.reserve(Nf);
	double f(0);
	for (int i = 0; i < Nf; i++) {
		f = fmin + 0.5 * df * (2 * i + 1);
		fValues.push_back(f);
	}

	// Calculate the probabilities
	std::vector<double> pValues; pValues.reserve(Nf);
	for (int i = 0; i < Nf; i++) pValues.push_back(0.0);
	std::vector<double> resfreq;
	for (size_t i = 0; i < spin_system.size(); ++i) {
		for (int f = 0; f < n_fields; ++f) {
			// Compute the resonance frequences for the single orientation of the magnetic field
			resfreq = spin_system[i].calculate_resfreq(fieldDirA[f], exp, urng);
			// Add the calculated frequencies to the spectrum
			for (size_t k = 0; k < spin_system[i].Ncomp; k++) {
				idx = static_cast<int>(resfreq[k] / df) - idx0;
				pValues[idx] += static_cast<double>(spin_system[i].Icomp[k]) / static_cast<double>(spin_system[i].Nstates);
			}
			resfreq.clear();
		}
	}
	// Normalize the probabilities
	double pMax = *max_element(pValues.begin(), pValues.end());
	for (int i = 0; i < Nf; i++) pValues[i] = pValues[i] / pMax;
	// Create the spectrum
	std::vector<std::vector<double>> spectrum; spectrum.reserve(2);
	spectrum.push_back(fValues);
	spectrum.push_back(pValues);
	return spectrum;
}

void PeldorCalculator::exitation_probabilities_spinA(std::vector<experiment> const& experiments, Spin const& spinA)
{
	// Initialize a random generator
	std::random_device rd;
	std::array<unsigned int, std::mt19937::state_size> seed_data;
	std::generate_n(seed_data.begin(), seed_data.size(), std::ref(rd));
	std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
	std::mt19937 urng(seq);

	// Calculate the excitation probabilities of spin A
	double gValue(0);
	gValuesA.reserve(n_fields);
	std::vector<double> resfreq; resfreq.reserve(spinA.Ncomp);
	double detProb(0), pumpProb(0);
	std::vector<double> detProbs; detProbs.reserve(n_fields);
	std::vector<double> pumpProbs; pumpProbs.reserve(n_fields);
	size_t numOfExp = experiments.size();
	detProbsA.reserve(numOfExp);
	pumpProbsA.reserve(numOfExp);
	for (size_t i = 0; i < numOfExp; ++i) {
		for (int f = 0; f < n_fields; ++f) {
			// Compute an effective g-factor of the spin A
			gValue = spinA.calculate_gfactor(fieldDirA[f], urng);
			// Compute resonance frequencies of the spin A
			resfreq = spinA.calculate_resfreq(fieldDirA[f], experiments[i], urng);
			// Compute the probability of spin A to be excited by the detection pulses
			detProb = excitation_probability_detection(resfreq, experiments[i], spinA);
			// Compute the probability of spin A to be excited by a pump pulse
			pumpProb = excitation_probability_pump(resfreq, experiments[i], spinA);
			// Store the calculated values
			gValuesA.push_back(gValue);
			detProbs.push_back(detProb);
			pumpProbs.push_back(pumpProb);
			resfreq.clear();
		}
		detProbsA.push_back(detProbs);
		pumpProbsA.push_back(pumpProbs);
		detProbs.clear();
		pumpProbs.clear();
	}
}

double PeldorCalculator::excitation_probability_detection(std::vector<double> const& freqs, experiment const& exp, Spin const& spin) const
{
	double detProb(0);
	double detPiHalfBW = 0.25 / exp.detPiHalfLength;
	double detPiBW = 0.50 / exp.detPiLength;
	double freqEff(0), freqEff2(0), prob(0);
	for (size_t k = 0; k < spin.Ncomp; k++) {
		if (detPiHalfBW == detPiBW) {
			freqEff = sqrt(pow(exp.detFreq - freqs[k], 2) + pow(detPiHalfBW, 2));
			prob = pow((detPiHalfBW / freqEff) * sin(2*PI * freqEff * exp.detPiHalfLength), 5);
		}
		else {
			freqEff = sqrt(pow(exp.detFreq - freqs[k], 2) + pow(detPiHalfBW, 2));
			freqEff2 = sqrt(pow(exp.detFreq - freqs[k], 2) + pow(detPiBW, 2));
			prob = (detPiHalfBW / freqEff) * sin(2*PI * freqEff * exp.detPiHalfLength) * 
				   pow((detPiBW / freqEff2) * sin(0.5 * 2*PI * freqEff2 * exp.detPiLength), 4);
		}
		detProb += fabs(static_cast<double>(spin.Icomp[k]) * prob);
	}
	detProb /= static_cast<double>(spin.Nstates);
	return detProb;
}

double PeldorCalculator::excitation_probability_pump(std::vector<double> const& freqs, experiment const& exp, Spin const& spin) const
{
	double pumpProb(0);
	double pumpPiBW = 0.5 / exp.pumpPiLength;
	double freqEff(0), prob(0);
	for (size_t k = 0; k < spin.Ncomp; k++) {
		freqEff = sqrt(pow(exp.pumpFreq - freqs[k], 2) + pow(pumpPiBW, 2));
		prob = pow((pumpPiBW / freqEff) * sin(0.5 * 2*PI * freqEff * exp.pumpPiLength), 2);
		pumpProb += fabs(static_cast<double>(spin.Icomp[k]) * prob);
	}
	pumpProb /= static_cast<double>(spin.Nstates);
	return pumpProb;
}

double PeldorCalculator::get_parameter_value(std::vector<double> const& opt_param_values, optimization_parameters const& opt_param, int const& nOffset, int const& idx) const
{
	double param_value(0);
	// if the parameter value is not optimized 
	if (opt_param.opt_param_numbers[idx] == -1) {
		// ... and if the parameter value is not fixed
		if (opt_param.fixed_param_numbers[idx] == -1) {
			if ((idx <= 23) || (idx == 25) || (idx == 26)) param_value = 0.0;
			else if (idx == 24) param_value = 1.0;
			else if (idx == 27) param_value = 1.0;
		}
		// ... and if the parameter value is fixed
		else {
			param_value = opt_param.fixed_param_values[opt_param.fixed_param_numbers[idx]];
		}
	}
	// if the parameter value is optimized 
	else {
		if ((idx == 27) && (opt_param.param_modes[idx] == 1)) {
			param_value = opt_param_values[opt_param.opt_param_numbers[idx] + nOffset];
		}
		else {
			param_value = opt_param_values[opt_param.opt_param_numbers[idx]];
		}
	}
	return param_value;
}

double PeldorCalculator::get_value_from_distribution(std::vector<double> const& opt_param_values, optimization_parameters const& opt_param,
	std::mt19937& urng, int const& idx_mean, int const& idx_width) const
{
	double mean(0), width(0), value(0);
	mean = get_parameter_value(opt_param_values, opt_param, 0, idx_mean);
	width = get_parameter_value(opt_param_values, opt_param, 0, idx_width);
	if (width == 0) {
		value = mean;
	}
	else {
		// Uniform distribution
		if (opt_param.param_modes[idx_mean] == 0) {
			std::uniform_real_distribution<> uniform_distr(mean - 0.5 * width, mean + 0.5 * width);
			value = uniform_distr(urng);
		}
		// Normal distribution
		else if (opt_param.param_modes[idx_mean] == 1) {
			std::normal_distribution<> normal_distr(mean, width);
			value = normal_distr(urng);
		}
	}
	return value;
}

std::vector<double> PeldorCalculator::direction_from_angles(double const& xi, double const& phi) const
{
	std::vector<double> dirVec; dirVec.reserve(3);
	dirVec.push_back(cos(phi) * sin(xi));
	dirVec.push_back(sin(phi) * sin(xi));
	dirVec.push_back(cos(xi));
	return dirVec;
}

std::tuple<double, std::vector<double>> PeldorCalculator::direction_from_vectors(double const& dist1, double const& dist2, 
	std::vector<double> const& dir1, std::vector<double> const& dir2) const
{
	double dist3(0);
	std::vector<double> dir3; dir3.reserve(3);
	dist3 = sqrt(pow((dist2 * dir2[0] - dist1 * dir1[0]), 2) + pow((dist2 * dir2[1] - dist1 * dir1[1]), 2) + pow((dist2 * dir2[2] - dist1 * dir1[2]), 2));
	dir3.push_back((dist2 * dir2[0] - dist1 * dir1[0]) / dist3);
	dir3.push_back((dist2 * dir2[1] - dist1 * dir1[1]) / dist3);
	dir3.push_back((dist2 * dir2[2] - dist1 * dir1[2]) / dist3);
	return std::make_tuple(dist3, dir3);
}

double PeldorCalculator::modulation_depth_correction(double const& w) const
{
	return exp(-pow((w / 16.0), 2));
}

std::vector<double> PeldorCalculator::calculate_peldor_signal(std::vector<double> const& opt_param_values, experiment const& exp,
	std::vector<Spin> const& spin_system, optimization_parameters const& opt_param, std::mt19937& urng, int const& nOffset) const
{
	// Initialize an array for the simulated signal
	size_t const nPoints = exp.timeValues.size();
	std::vector<double> signalValues; signalValues.reserve(nPoints);
	for (size_t t = 0; t < nPoints; ++t) signalValues.push_back(0.0);
	// Initialize the rotation matrix 
	std::shared_ptr<RotationMatrix> RM(new RotationMatrix(0.0, 0.0, 0.0));
	// Calculate the PELDOR signal for the two-spin sytem
	if (spin_system.size() == 2) {
		// Initialize some variables
		std::uniform_real_distribution<> uniform_distr(0.0, 1.0);
		double dist(0), xi(0), phi(0), alpha(0), beta(0), gamma(0), J(0), ratio(0), pumpEfficiency(0);
		std::vector<double> fieldDirB; fieldDirB.reserve(3);
		std::vector<double> resfreqB; resfreqB.reserve(spin_system[1].Ncomp);
		std::vector<double> distDirAB; distDirAB.reserve(3);
		double gValueA(0), gValueB(0);
		double detProbA(0), detProbB(0), pumpProbA(0), pumpProbB(0);
		double amplitude(0), cosDipolarAngle(0), fdd(0), wmod(0), lambda(0);
		bool excited_AB(false), excited_BA(false);
		// Set the ratio between the first and the second distances
		ratio = get_parameter_value(opt_param_values, opt_param, nOffset, 24);
		// Set the inversion efficiency of the pump pulse
		pumpEfficiency = get_parameter_value(opt_param_values, opt_param, nOffset, 27);
		// Iterate over all directions of the magnetic field
		for (int f = 0; f < n_fields; ++f) {
			// Set the values of all geometric parameters and the J coupling constant
			if (uniform_distr(urng) <= ratio) {
				dist = get_value_from_distribution(opt_param_values, opt_param, urng, 0, 1);
				xi = get_value_from_distribution(opt_param_values, opt_param, urng, 2, 3);
				phi = get_value_from_distribution(opt_param_values, opt_param, urng, 4, 5);
				alpha = get_value_from_distribution(opt_param_values, opt_param, urng, 6, 7);
				beta = get_value_from_distribution(opt_param_values, opt_param, urng, 8, 9);
				gamma = get_value_from_distribution(opt_param_values, opt_param, urng, 10, 11);
			}
			else {
				dist = get_value_from_distribution(opt_param_values, opt_param, urng, 12, 13);
				xi = get_value_from_distribution(opt_param_values, opt_param, urng, 14, 15);
				phi = get_value_from_distribution(opt_param_values, opt_param, urng, 16, 17);
				alpha = get_value_from_distribution(opt_param_values, opt_param, urng, 18, 19);
				beta = get_value_from_distribution(opt_param_values, opt_param, urng, 20, 21);
				gamma = get_value_from_distribution(opt_param_values, opt_param, urng, 22, 23);
			}
			// Set the value of the J coupling constant
			J = get_value_from_distribution(opt_param_values, opt_param, urng, 25, 26);

			// Read out the g value of spin A
			gValueA = gValuesA[f];
			// Read out the probability of spin A to be excited by the detection pulses
			detProbA = detProbsA[nOffset][f];
			// Read out the probability of spin A to be excited by the pump pulse
			pumpProbA = pumpProbsA[nOffset][f];

			// Rotation matrix between the spin A and spin B frames
			RM->reset_angles(alpha, beta, gamma);
			// Calculate the direction of the magnetic field in the spin B frame
			fieldDirB = RM->dot_product(fieldDirA[f], true);
			// Calculate the effective g-factor of  spin B
			gValueB = spin_system[1].calculate_gfactor(fieldDirB, urng);
			// Calculate the resonance frequencies of spin B
			resfreqB = spin_system[1].calculate_resfreq(fieldDirB, exp, urng);
			// Calculate the probability of spin B to be excited by the detection pulses
			detProbB = excitation_probability_detection(resfreqB, exp, spin_system[1]);
			// Calculate the probability of spin B to be excited by the pump pulse
			if (detProbA > threshold) pumpProbB = excitation_probability_pump(resfreqB, exp, spin_system[1]);
			else pumpProbB = 0.0;

			// Calculate the amplitude of the PELDOR signal
			if (detProbA > threshold) amplitude += detProbA;
			if (detProbB > threshold) amplitude += detProbB;

			// Determine whether the spin pair is excited
			excited_AB = ((detProbA > threshold) && (pumpProbB > threshold));
			excited_BA = ((detProbB > threshold) && (pumpProbA > threshold));

			// Calculate the oscillating part of the PELDOR signal
			if (excited_AB || excited_BA) {
				// The direction of the distance vector
				distDirAB = direction_from_angles(xi, phi);
				// The cosine of the dipolar angle
				cosDipolarAngle = fieldDirA[f][0] * distDirAB[0] + fieldDirA[f][1] * distDirAB[1] + fieldDirA[f][2] * distDirAB[2];
				// The dipolar frequency
				fdd = 52.04 * gValueA * gValueB * (1.0 - 3.0 * cosDipolarAngle * cosDipolarAngle) / (2.0023 * 2.0023 * pow(dist, 3));
				// The echo modulation frequency
				wmod = 2 * PI * (fdd + J);
				// The modulation amplutude
				if ((excited_AB) && (excited_BA))       lambda = (detProbA * pumpProbB + detProbB * pumpProbA) * pumpEfficiency;
				else if ((excited_AB) && (!excited_BA)) lambda = (detProbA * pumpProbB) * pumpEfficiency;
				else if ((!excited_AB) && (excited_BA)) lambda = (detProbB * pumpProbA) * pumpEfficiency;
				// The oscillating part of the PELDOR signal
				for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda * (1 - cos(wmod * exp.timeValues[t]));
			}

			// Refresh some variables
			fieldDirB.clear();
			resfreqB.clear();
			distDirAB.clear();
		}
		// Calculate the entire PELDOR signal and normalize it
		double norm = 1.0 / amplitude;
		for (size_t t = 0; t < nPoints; ++t) signalValues[t] = (amplitude - signalValues[t]) * norm;
	}
	// Calculate the PELDOR signal for the three-spin sytem
	else if (spin_system.size() == 3) {
		// Initialize some variables
		double dist_AB(0), xi_AB(0), phi_AB(0), alpha_AB(0), beta_AB(0), gamma_AB(0);
		double dist_AC(0), xi_AC(0), phi_AC(0), alpha_AC(0), beta_AC(0), gamma_AC(0);
        double pumpEfficiency(0), J(0), dist_BC(0);
		double gValueA(0), gValueB(0), gValueC(0);
		double detProbA(0), detProbB(0), detProbC(0), pumpProbA(0), pumpProbB(0), pumpProbC(0);
		double amplitude(0);
		bool excited_AB(false), excited_AC(false), excited_BC(false), excited_BA(false), excited_CA(false), excited_CB(false);
		double cosDipolarAngle_AB(0), cosDipolarAngle_AC(0), cosDipolarAngle_BC(0);
		double fdd_AB(0), fdd_AC(0), fdd_BC(0);
		double wmod_AB(0), wmod_AC(0), wmod_BC(0), wmod_sum(0), wmod_dif(0);
		double lambda_AB(0), lambda_AC(0), lambda_BC(0), lambda_sum(0), lambda_dif(0);
		std::vector<double> fieldDirB; fieldDirB.reserve(3);
		std::vector<double> fieldDirC; fieldDirC.reserve(3);
		std::vector<double> resfreqB; resfreqB.reserve(spin_system[1].Ncomp);
		std::vector<double> resfreqC; resfreqC.reserve(spin_system[2].Ncomp);
		std::vector<double> distDirAB; distDirAB.reserve(3);
		std::vector<double> distDirAC; distDirAC.reserve(3);
		std::vector<double> distDirBC; distDirBC.reserve(3);
		std::tuple<double, std::vector<double>> distVecBC;
		// Set the inversion efficiency of the pump pulse
		pumpEfficiency = get_parameter_value(opt_param_values, opt_param, nOffset, 27);
		// Iterate over all directions of the magnetic field
		for (int f = 0; f < n_fields; ++f) {
			// Set the values of all geometric parameters
			dist_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 0, 1);
			xi_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 2, 3);
			phi_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 4, 5);
			alpha_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 6, 7);
			beta_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 8, 9);
			gamma_AB = get_value_from_distribution(opt_param_values, opt_param, urng, 10, 11);
			dist_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 12, 13);
			xi_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 14, 15);
			phi_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 16, 17);
			alpha_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 18, 19);
			beta_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 20, 21);
			gamma_AC = get_value_from_distribution(opt_param_values, opt_param, urng, 22, 23);
			// Set the value of the J coupling constant
			J = get_value_from_distribution(opt_param_values, opt_param, urng, 25, 26);

			// Read out the g value of spin A
			gValueA = gValuesA[f];
			// Read out the probability of spin A to be excited by the detection pulses
			detProbA = detProbsA[nOffset][f];
			// Read out the probability of spin A to be excited by the pump pulse
			pumpProbA = pumpProbsA[nOffset][f];
			// Spin B:
			// Rotation matrix between the spin A and spin B frames
			RM->reset_angles(alpha_AB, beta_AB, gamma_AB);
			// Calculate the direction of the magnetic field in the spin B frame
			fieldDirB = RM->dot_product(fieldDirA[f], true);
			// Calculate the effective g-factor of  spin B
			gValueB = spin_system[1].calculate_gfactor(fieldDirB, urng);
			// Calculate the resonance frequencies of spin B
			resfreqB = spin_system[1].calculate_resfreq(fieldDirB, exp, urng);
			// Compute the probability of spin B to be excited by the detection pulses
			detProbB = excitation_probability_detection(resfreqB, exp, spin_system[1]);	

			// Rotation matrix between the spin A and spin C frames
			RM->reset_angles(alpha_AC, beta_AC, gamma_AC);
			// Calculate the direction of the magnetic field in the spin C frame
			fieldDirC = RM->dot_product(fieldDirA[f], true);
			// Calculate the effective g-factor of  spin C
			gValueC = spin_system[2].calculate_gfactor(fieldDirC, urng);
			// Calculate the resonance frequencies of spin C
			resfreqC = spin_system[2].calculate_resfreq(fieldDirC, exp, urng);
			// Calculate the probability of spin C to be excited by the detection pulses
			detProbC = excitation_probability_detection(resfreqC, exp, spin_system[2]);
			// Calculate the probability of spin B to be excited by the pump pulse
			if ((detProbA > threshold) || (detProbC > threshold)) pumpProbB = excitation_probability_pump(resfreqB, exp, spin_system[1]);
			else pumpProbB = 0.0;
			// Calculate the probability of spin C to be excited by the pump pulse
			if ((detProbA > threshold) || (detProbB > threshold)) pumpProbC = excitation_probability_pump(resfreqC, exp, spin_system[2]);
			else pumpProbB = 0.0;

			// The direction of the distance vector between spins A and B
			distDirAB = direction_from_angles(xi_AB, phi_AB);
			// The direction of the distance vector between spins A and C
			distDirAC = direction_from_angles(xi_AC, phi_AC);
			// The direction of the distance vector between spins B and C
			distVecBC = direction_from_vectors(dist_AB, dist_AC, distDirAB, distDirAC);
			dist_BC = std::get<0>(distVecBC);
			distDirBC = std::get<1>(distVecBC);

			// Calculate the amplitude of the PELDOR signal
			if (detProbA > threshold) amplitude += detProbA;
			if (detProbB > threshold) amplitude += detProbB;
			if (detProbC > threshold) amplitude += detProbC;

			// Determine which spin pairs are exited
			excited_AB = ((detProbA > threshold) && (pumpProbB > threshold));
			excited_AC = ((detProbA > threshold) && (pumpProbC > threshold));
			excited_BC = ((detProbB > threshold) && (pumpProbC > threshold));
			excited_BA = ((detProbB > threshold) && (pumpProbA > threshold));
			excited_CA = ((detProbC > threshold) && (pumpProbA > threshold));
			excited_CB = ((detProbC > threshold) && (pumpProbB > threshold));

			// Calculate additional terms to the amplitude of the PELDOR signal
			// No multi-spin effects
			if (!opt_param.multiSpin) {
				if (excited_AB || excited_BA) {
					// The cosine of the dipolar angle for spins A and B
					cosDipolarAngle_AB = fieldDirA[f][0] * distDirAB[0] + fieldDirA[f][1] * distDirAB[1] + fieldDirA[f][2] * distDirAB[2];
					// The dipolar frequency for spins A and B
					fdd_AB = 52.04 * gValueA * gValueB * (1.0 - 3.0 * cosDipolarAngle_AB * cosDipolarAngle_AB) / (2.0023 * 2.0023 * pow(dist_AB, 3));
					// The echo modulation frequency for spins A and B
					wmod_AB = 2 * PI * (fdd_AB + J);
					// The modulation amplutude
					if ((excited_AB) && (excited_BA))       lambda_AB = (detProbA * pumpProbB + detProbB * pumpProbA) * pumpEfficiency;
					else if ((excited_AB) && !(excited_BA)) lambda_AB = (detProbA * pumpProbB) * pumpEfficiency;
					else if (!(excited_AB) && (excited_BA)) lambda_AB = (detProbB * pumpProbA) * pumpEfficiency;
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_AB * (1 - cos(wmod_AB * exp.timeValues[t]));
				}
				if (excited_AC || excited_CA) {
					// The cosine of the dipolar angle for spins A and C
					cosDipolarAngle_AC = fieldDirA[f][0] * distDirAC[0] + fieldDirA[f][1] * distDirAC[1] + fieldDirA[f][2] * distDirAC[2];
					// The dipolar frequency for spins A and C
					fdd_AC = 52.04 * gValueA * gValueC * (1.0 - 3.0 * cosDipolarAngle_AC * cosDipolarAngle_AC) / (2.0023 * 2.0023 * pow(dist_AC, 3));
					// The echo modulation frequency for spins A and C
					wmod_AC = 2 * PI * (fdd_AC + J);
					// The oscillating part of the PELDOR signal
					if ((excited_AC) && (excited_CA))       lambda_AC = (detProbA * pumpProbC + detProbC * pumpProbA) * pumpEfficiency;
					else if ((excited_AC) && !(excited_CA)) lambda_AC = (detProbA * pumpProbC) * pumpEfficiency;
					else if (!(excited_AC) && (excited_CA)) lambda_AC = (detProbC * pumpProbA) * pumpEfficiency;
					// The signal values
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_AC * (1 - cos(wmod_AC * exp.timeValues[t]));
				}
				if (excited_BC || excited_CB) {
					// The cosine of the dipolar angle for spins B and C
					cosDipolarAngle_BC = fieldDirA[f][0] * distDirBC[0] + fieldDirA[f][1] * distDirBC[1] + fieldDirA[f][2] * distDirBC[2];
					// The dipolar frequency for spins B and C
					fdd_BC = 52.04 * gValueB * gValueC * (1.0 - 3.0 * cosDipolarAngle_BC * cosDipolarAngle_BC) / (2.0023 * 2.0023 * pow(dist_BC, 3));
					// The echo modulation frequency for spins B and C
					wmod_BC = 2 * PI * (fdd_BC + J);
					// The modulation amplutude
					if ((excited_BC) && (excited_CB))       lambda_BC = (detProbB * pumpProbC + detProbC * pumpProbB) * pumpEfficiency;
					else if ((excited_BC) && !(excited_CB)) lambda_BC = (detProbB * pumpProbC) * pumpEfficiency;
					else if (!(excited_BC) && (excited_CB)) lambda_BC = (detProbC * pumpProbB) * pumpEfficiency;
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_BC * (1 - cos(wmod_BC * exp.timeValues[t]));
				}
			}
			// With multi-spin effects
			else {
				if (excited_AB || excited_BA) {
					// The cosine of the dipolar angle for spins A and B
					cosDipolarAngle_AB = fieldDirA[f][0] * distDirAB[0] + fieldDirA[f][1] * distDirAB[1] + fieldDirA[f][2] * distDirAB[2];
					// The dipolar frequency for spins A and B
					fdd_AB = 52.04 * gValueA * gValueB * (1.0 - 3.0 * cosDipolarAngle_AB * cosDipolarAngle_AB) / (2.0023 * 2.0023 * pow(dist_AB, 3));
					// The echo modulation frequency for spins A and B
					wmod_AB = 2 * PI * (fdd_AB + J);
					// The modulation amplutude
					if ((excited_AB) && (excited_BA))       lambda_AB = (detProbA * pumpProbB + detProbB * pumpProbA) * pumpEfficiency;
					else if ((excited_AB) && !(excited_BA)) lambda_AB = (detProbA * pumpProbB) * pumpEfficiency;
					else if (!(excited_AB) && (excited_BA)) lambda_AB = (detProbB * pumpProbA) * pumpEfficiency;
					if (pumpProbC > threshold) lambda_AB *= (1 - pumpProbC);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_AB * (1 - cos(wmod_AB * exp.timeValues[t]));	
				}
				if (excited_AC || excited_CA) {
					// The cosine of the dipolar angle for spins A and C
					cosDipolarAngle_AC = fieldDirA[f][0] * distDirAC[0] + fieldDirA[f][1] * distDirAC[1] + fieldDirA[f][2] * distDirAC[2];
					// The dipolar frequency for spins A and C
					fdd_AC = 52.04 * gValueA * gValueC * (1.0 - 3.0 * cosDipolarAngle_AC * cosDipolarAngle_AC) / (2.0023 * 2.0023 * pow(dist_AC, 3));
					// The echo modulation frequency for spins A and C
					wmod_AC = 2 * PI * (fdd_AC + J);
					// The modulation amplutude
					if ((excited_AC) && (excited_CA))       lambda_AC = (detProbA * pumpProbC + detProbC * pumpProbA) * pumpEfficiency;
					else if ((excited_AC) && !(excited_CA)) lambda_AC = (detProbA * pumpProbC) * pumpEfficiency;
					else if (!(excited_AC) && (excited_CA)) lambda_AC = (detProbC * pumpProbA) * pumpEfficiency;
					if (pumpProbB > threshold) lambda_AC *= (1 - pumpProbB);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_AC * (1 - cos(wmod_AC * exp.timeValues[t]));
				}
				if (excited_BC || excited_CB) {
					// The cosine of the dipolar angle for spins B and C
					cosDipolarAngle_BC = fieldDirA[f][0] * distDirBC[0] + fieldDirA[f][1] * distDirBC[1] + fieldDirA[f][2] * distDirBC[2];
					// The dipolar frequency for spins B and C
					fdd_BC = 52.04 * gValueB * gValueC * (1.0 - 3.0 * cosDipolarAngle_BC * cosDipolarAngle_BC) / (2.0023 * 2.0023 * pow(dist_BC, 3));
					// The echo modulation frequency for spins B and C
					wmod_BC = 2 * PI * (fdd_BC + J);
					// The modulation amplutude
					if ((excited_BC) && (excited_CB))       lambda_BC = (detProbB * pumpProbC + detProbC * pumpProbB) * pumpEfficiency;
					else if ((excited_BC) && !(excited_CB)) lambda_BC = (detProbB * pumpProbC) * pumpEfficiency;
					else if (!(excited_BC) && (excited_CB)) lambda_BC = (detProbC * pumpProbB) * pumpEfficiency;
					if (pumpProbA > threshold) lambda_BC *= (1 - pumpProbA);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_BC * (1 - cos(wmod_BC * exp.timeValues[t]));
				}
				if (excited_AB && excited_AC) {
					// The sum and difference frequencies
					wmod_sum = wmod_AB + wmod_AC;
					wmod_dif = fabs(wmod_AB - wmod_AC);
					// The modulation amplutude
					lambda_dif = 0.5 * detProbA * pumpProbB * pumpProbC;
					lambda_sum = lambda_dif * modulation_depth_correction(wmod_sum);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_dif * (1 - cos(wmod_dif * exp.timeValues[t])) + lambda_sum * (1 - cos(wmod_sum * exp.timeValues[t]));
				}
				if (excited_BA && excited_BC) {
					// The sum and difference frequencies
					wmod_sum = wmod_AB + wmod_BC;
					wmod_dif = fabs(wmod_AB - wmod_BC);
					// The modulation amplutude
					lambda_dif = 0.5 * detProbB * pumpProbA * pumpProbC;
					lambda_sum = lambda_dif * modulation_depth_correction(wmod_sum);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_dif * (1 - cos(wmod_dif * exp.timeValues[t])) + lambda_sum * (1 - cos(wmod_sum * exp.timeValues[t]));
				}
				if (excited_CA && excited_CB) {
					// The sum and difference frequencies
					wmod_sum = wmod_AC + wmod_BC;
					wmod_dif = fabs(wmod_AC - wmod_BC);
					// The modulation amplutude
					lambda_dif = 0.5 * detProbC * pumpProbA * pumpProbB;
					lambda_sum = lambda_dif * modulation_depth_correction(wmod_sum);
					// The oscillating part of the PELDOR signal
					for (size_t t = 0; t < nPoints; ++t) signalValues[t] += lambda_dif * (1 - cos(wmod_dif * exp.timeValues[t])) + lambda_sum * (1 - cos(wmod_sum * exp.timeValues[t]));
				}

			}

			// Refresh some variables
			fieldDirB.clear();
			fieldDirC.clear();
			resfreqB.clear();
			resfreqC.clear();
			distDirAB.clear();
			distDirAC.clear();
			distDirBC.clear();
		}
		// Calculate the entire PELDOR signal and normalize it
		double norm = 1.0 / amplitude;
		for (size_t t = 0; t < nPoints; ++t) signalValues[t] = (amplitude - signalValues[t]) * norm;
	}
	return signalValues;
}

double PeldorCalculator::rmsd(std::vector<double> const& x, std::vector<double> const& y) const
{
	double rmsd(0);
	size_t const nPoints = y.size();
	for (size_t t = 0; t < nPoints; ++t) rmsd += pow((x[t] - y[t]), 2);
	rmsd = sqrt(rmsd / static_cast<double>(nPoints));
	return rmsd;
}

double PeldorCalculator::rmsd_over_pearson(std::vector<double> const& x, std::vector<double> const& y) const
{
	size_t const nPoints = y.size();
	double y_mean(0), y_variance(0), x_mean(0), x_variance(0);
	double rmsd(0), pcc(0), rmsd_pcc(0);
	for (size_t t = 0; t < nPoints; ++t) {
		x_mean += x[t];
		y_mean += y[t];
	}
	x_mean /= static_cast<double>(nPoints);
	y_mean /= static_cast<double>(nPoints);
	for (size_t t = 0; t < nPoints; ++t) {
		x_variance += (x[t] - x_mean) * (x[t] - x_mean);
		y_variance += (y[t] - y_mean) * (y[t] - y_mean);
		rmsd += (x[t] - y[t]) * (x[t] - y[t]);
		pcc += (x[t] - x_mean) * (y[t] - y_mean);
	}
	rmsd = sqrt(rmsd / static_cast<double>(nPoints));
	pcc /= sqrt(x_variance * y_variance);
	rmsd_pcc = rmsd / pcc;
	return rmsd_pcc;
}

double PeldorCalculator::pearson(std::vector<double> const& x, std::vector<double> const& y) const
{
	size_t const nPoints = y.size();
	double y_mean(0), y_variance(0), x_mean(0), x_variance(0);
	double pcc(0);
	for (size_t t = 0; t < nPoints; ++t) {
		x_mean += x[t];
		y_mean += y[t];
	}
	x_mean /= static_cast<double>(nPoints);
	y_mean /= static_cast<double>(nPoints);
	for (size_t t = 0; t < nPoints; ++t) {
		x_variance += (x[t] - x_mean) * (x[t] - x_mean);
		y_variance += (y[t] - y_mean) * (y[t] - y_mean);
		pcc += (x[t] - x_mean) * (y[t] - y_mean);
	}
	pcc /= sqrt(x_variance * y_variance);
	return pcc;
}

std::vector<double> PeldorCalculator::angles_from_direction(std::vector<double> const& v) const
{
	double xi(0), phi(0);
	std::vector<double> angles; angles.reserve(2);
	xi = acos(v[2]);
	if ((xi == 0) || (xi == PI)) phi = 0;
	else {
		if (v[0] == 0) {
			if (v[1] == sin(xi)) phi = 0.5*PI;
			else                 phi = -0.5*PI;
		}
		else phi = atan2(v[1], v[0]);
	}
	if (xi < 0) xi += PI;
	if (phi < 0) phi += 2 * PI;
	angles.push_back(xi);
	angles.push_back(phi);
	return angles;
}

std::vector<double> PeldorCalculator::angles_from_rotation_matrix(std::vector<std::vector<double>> const& RM) const
{
	double alpha(0), beta(0), gamma(0); 
	std::vector<double> angles; angles.reserve(3);
	beta = acos(RM[2][2]);
	if (beta == 0) {
		gamma = 0;
		alpha = atan2(RM[1][0], RM[1][1]);
	}
	else {
		gamma = atan2(RM[2][0] / sin(beta), RM[2][1] / sin(beta));
		alpha = atan2(RM[0][2] / sin(beta), -RM[1][2] / sin(beta));
	}
	if (alpha < 0) alpha += 2*PI;
	if (beta < 0)  beta += PI;
	if (gamma < 0) gamma += 2*PI;
	angles.push_back(alpha);
	angles.push_back(beta);
	angles.push_back(gamma);
	return angles;
}

std::vector<std::vector<double>> PeldorCalculator::symmetric_angles(double const& xi, double const& phi, double const& alpha, double const& beta, 
	double const& gamma) const
{
	std::vector<std::vector<double>> sym_angles; sym_angles.reserve(16);
	// The direction of the distance vector
	std::vector<double> distDir; distDir.reserve(3);
	distDir = direction_from_angles(xi, phi);
	// Rotation matrices for the rotation between the spin frames
	std::shared_ptr<RotationMatrix> RM(new RotationMatrix(alpha, beta, gamma));
	// Rotation matrices for 180° rotation about the x,y,z axes of the coordinate system
	std::shared_ptr<RotationMatrix> RI(new RotationMatrix(0.0, 0.0, 0.0));
	std::shared_ptr<RotationMatrix> Rx(new RotationMatrix(0.0, PI, 0.0));
	std::shared_ptr<RotationMatrix> Ry(new RotationMatrix(PI, PI, 0.0));
	std::shared_ptr<RotationMatrix> Rz(new RotationMatrix(PI, 0.0, 0.0));
	std::vector<std::shared_ptr<RotationMatrix>> transformations; transformations.reserve(4);
	transformations.push_back(RI);
	transformations.push_back(Rx);
	transformations.push_back(Ry);
	transformations.push_back(Rz);
	// Calculate the symmetry-related sets of fitting parameters
	std::vector<double> distDir_rotated; distDir_rotated.reserve(3);
	std::vector<std::vector<double>> RM_rotated; RM_rotated.reserve(3);
	std::vector<double> new_xi_phi; new_xi_phi.reserve(2);
	std::vector<double> new_alpha_beta_gamma; new_alpha_beta_gamma.reserve(3);
	std::vector<double> new_angles; new_angles.reserve(5);
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			// Rotate the spin A frame
			distDir_rotated = transformations[i]->dot_product(distDir, true);
			// Rotate the spin B frame
			RM_rotated = transformations[i]->matrix_product(RM->matrix_product(transformations[j]->R, false), true);
			// Calculate the new values of the angles
			new_xi_phi = angles_from_direction(distDir_rotated);
			new_alpha_beta_gamma = angles_from_rotation_matrix(RM_rotated);
			// Save the calculated angles
			new_angles.push_back(new_xi_phi[0]);
			new_angles.push_back(new_xi_phi[1]);
			new_angles.push_back(new_alpha_beta_gamma[0]);
			new_angles.push_back(new_alpha_beta_gamma[1]);
			new_angles.push_back(new_alpha_beta_gamma[2]);
			sym_angles.push_back(new_angles);
			new_angles.clear();
		}
	}
	return sym_angles;
}

std::vector<std::vector<double>> PeldorCalculator::calculate_symmetric_angles(std::vector<double> const& opt_param_values, 
	std::vector<experiment> const& experiments, std::vector<Spin> const& spin_system, optimization_parameters const& opt_param) const
{
	// Initialize a vector for the symmetric parameters
	std::vector<std::vector<double>> sym_angles_all; sym_angles_all.reserve(16);
	// Read out the optimized values of the fitting parameters
	double xi_AB = get_parameter_value(opt_param_values, opt_param, 0, 2);
	double phi_AB = get_parameter_value(opt_param_values, opt_param, 0, 4);
	double alpha_AB = get_parameter_value(opt_param_values, opt_param, 0, 6);
	double beta_AB = get_parameter_value(opt_param_values, opt_param, 0, 8);
	double gamma_AB = get_parameter_value(opt_param_values, opt_param, 0, 10);
	double xi_AC = get_parameter_value(opt_param_values, opt_param, 0, 14);
	double phi_AC = get_parameter_value(opt_param_values, opt_param, 0, 16);
	double alpha_AC = get_parameter_value(opt_param_values, opt_param, 0, 18);
	double beta_AC = get_parameter_value(opt_param_values, opt_param, 0, 20);
	double gamma_AC = get_parameter_value(opt_param_values, opt_param, 0, 22);
	// Calculate the symmetry-related sets of angles
	std::vector<std::vector<double>> sym_angles_AB = symmetric_angles(xi_AB, phi_AB, alpha_AB, beta_AB, gamma_AB);
	std::vector<std::vector<double>> sym_angles_AC = symmetric_angles(xi_AC, phi_AC, alpha_AC, beta_AC, gamma_AC);
	std::vector<double> row; row.reserve(10);
	for (size_t i = 0; i < 16; ++i) {
		for (size_t j = 0; j < 5; ++j) row.push_back(sym_angles_AB[i][j]);
		for (size_t j = 0; j < 5; ++j) row.push_back(sym_angles_AC[i][j]);
		//for (size_t j = 0; j < 10; ++j) std::cout << sym_angles_AB[i][j] << std::endl;
		sym_angles_all.push_back(row);
		row.clear();
	}
	return sym_angles_all;
}