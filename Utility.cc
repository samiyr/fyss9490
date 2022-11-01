#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <gsl/gsl_monte_vegas.h>

/// Container for the output data.
struct Row {
	/// Lower bound for the invariant mass bin.
	double sqrt_s_hat_min;
	/// Upper bound for the invariant mass bin.
	double sqrt_s_hat_max;
	/// Value of the cross section, in picobarns.
	double cross_section;
	/// Absolute error estimate of the cross section, in picobarns.
	double error;
	/// The chi^2 value associated with the integral value.
	double chi_squared;
	/// The number of rebinning iterations it took to obtain the integral.
	size_t iterations;

	friend std::ostream& operator<<(std::ostream& os, Row const & r) {
		return os << r.sqrt_s_hat_min << "," << r.sqrt_s_hat_max << "," << r.cross_section << "," << r.error << "," << r.chi_squared << "," << r.iterations;
	}
};

/// Container for parameters given to the VEGAS algorithm.
/// See the type gsl_monte_vegas_params in GSL.
struct VegasParameters {
	/// Stiffness of the rebinning algorithm.
	double alpha;
	/// Number of iterations per rebinning.
	size_t iterations;
	/// Importance sampling mode.
	int mode;
	/// Verbosity level.
	int verbose;

	/// Applies the parameters to an instance of gsl_monte_vegas_params.
	void apply_to(gsl_monte_vegas_params *vegas) {
		vegas->alpha = alpha;
		vegas->iterations = iterations;
		vegas->mode = mode;
		vegas->verbose = verbose;
	}
	
	friend std::ostream& operator<<(std::ostream& os, VegasParameters const &p) {
		os << "==========================" << std::endl;
		os << "====== VEGAS PARAMS ======" << std::endl;
		os << "==========================" << std::endl;
		os << "alpha = " << p.alpha << std::endl;
		os << "iterations = " << p.iterations << std::endl;
		os << "mode = " << p.mode << std::endl;
		os << "verbose = " << p.verbose << std::endl;
		os << "==========================" << std::endl << std::endl;
		return os;
	}
};

/// Container for parameters of the integration routine.
struct IntegrationRoutineParameters {
	/// Number of sampling points for VEGAS.
	size_t points;
	/// The cut-off value for the integral.
	/// When the integral surpasses this value, the 
	/// integration grid is deemed sufficiently adapted.
	double integration_phase_cutoff;
	/// Maximum value of | chi^2 - 1|. When the chi^2-value
	/// goes under this value, and the relative error is small
	/// enough (see max_relative_error), the integration stops.
	double max_chi_sq_deviation;
	/// Maximum value for the relative error of the integral.
	/// When the relative error goes under this value, and the
	/// chi^2 value is sufficiently close to 1 (see max_chi_sq_deviation),
	/// the integration stops.
	double max_relative_error;
	/// Maximum number of rebinning iterations for the fail-safe
	/// mechanism. If the number of rebinning iterations exceeds this value,
	/// the integration is stopped and the integral with a chi^2-value closest
	/// to 1 is returned.
	size_t max_iterations;
	/// Seed for the GSL random number generator.
	int seed;
	/// Maximum number of threads used in the integration routine.
	/// A value of 1 disables parallelization.
	int thread_count;
	/// Output filename.
	std::string filename;

	/// Checks whether parallelization should be enabled.
	/// See thread_count.
	const bool parallelize() const {
		return thread_count > 1;
	}

	friend std::ostream& operator<<(std::ostream& os, IntegrationRoutineParameters const &p) {
		os << "==========================" << std::endl;
		os << "=== INTEGRATION PARAMS ===" << std::endl;
		os << "==========================" << std::endl;
		os << "points = " << p.points << std::endl;
		os << "integration_phase_cutoff = " << p.integration_phase_cutoff << std::endl;
		os << "max_chi_sq_deviation = " << p.max_chi_sq_deviation << std::endl;
		os << "max_relative_error = " << p.max_relative_error << std::endl;
		os << "max_iterations = " << p.max_iterations << std::endl;
		os << "seed = " << p.seed << std::endl;
		os << "thread_count = " << p.thread_count << std::endl;
		os << "filename = " << p.filename << std::endl;
		os << "==========================" << std::endl << std::endl;
		return os;
	}
};

/// Container for parameters common to both the EPA calculation and
/// the full matrix element calculation.
struct CommonParameters {
	/// Square root of the Mandelstam s.
	double sqrt_s;
	/// Lepton mass.
	double m;
	/// Hadron mass.
	double M;
	/// QED coupling constant.
	double alpha;
	/// Dipole form factor parameter.
	double lambda2;
	/// Total magnetic moment of the hadron.
	double mu;

	/// Lower bound for the invariant mass bin.
	double sqrt_s_hat_min;
	/// Upper bound for the invariant mass bin.
	double sqrt_s_hat_max;
	/// Difference between consecutive invariant mass bins.
	double step_size;
	/// Calculates the number of invariant mass bins.
	const int step_count() const {
		return int(ceil((sqrt_s_hat_max - sqrt_s_hat_min) / step_size));
	}

	friend std::ostream& operator<<(std::ostream& os, CommonParameters const &p) {
		os << "==========================" << std::endl;
		os << "====== COMMON PARAMS =====" << std::endl;
		os << "==========================" << std::endl;
		os << "sqrt_s = " << p.sqrt_s << std::endl;
		os << "m = " << p.m << std::endl;
		os << "M = " << p.M << std::endl;
		os << "alpha = " << p.alpha << std::endl;
		os << "lambda2 = " << p.lambda2 << std::endl;
		os << "mu = " << p.mu << std::endl;
		os << "sqrt_s_hat_min = " << p.sqrt_s_hat_min << std::endl;
		os << "sqrt_s_hat_max = " << p.sqrt_s_hat_max << std::endl;
		os << "step_size = " << p.step_size << std::endl;
		os << "step_count = " << p.step_count() << std::endl;
		os << "==========================" << std::endl << std::endl;
		return os;
	}
};

namespace Constants {
	static double muon_mass = 0.105658372;
	static double proton_mass = 0.938272013;
	static double alpha_qed = 0.0072973525693;
	static double proton_magnetic_moment = 2.79;
	static double proton_dipole_lambda2 = 0.71;
	static double GeV2ToPb = 0.389379e9;
}

#endif // UTILITY_H