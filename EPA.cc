#ifndef EPA_H
#define EPA_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_monte_vegas.h>
#include "Utility.cc"
#include <vector>

#define DEBUG false

enum class Flux {
	full, no_magnetic_moment, no_mass_magnetic_moment
};

struct EPAParameters {
	VegasParameters vp;
	IntegrationRoutineParameters ip;
	CommonParameters p;
	Flux flux;

	friend std::ostream& operator<<(std::ostream& os, EPAParameters const &p) {
		return os << p.vp << p.ip << p.p;
	}
};

const double f1(const double y, const double Q2min, const double alpha, const double M, const double lambda2) {
	const double prefactor = alpha / (2.0 * M_PI);
	const double splitting = (1.0 + std::pow(1.0 - y, 2)) / y;
	const double A = 1.0 + lambda2 / Q2min;
	const double A_term = log(A) - 11.0 / 6.0 + 3.0 / A - 3.0 / (2.0 * std::pow(A, 2)) + 1.0 / (3.0 * std::pow(A, 3));

	return prefactor * splitting * A_term;
}

const double f2(const double y, const double Q2min, const double alpha, const double M, const double lambda2) {
	const double prefactor = alpha / (2.0 * M_PI);
	const double prefactor_2 = (2 * y * std::pow(M, 2)) / lambda2;
	const double A = 1.0 + lambda2 / Q2min;
	const double A_term = A - 4.0 * log(A) + 10.0 / 3.0 - 6.0 / A + 2.0 / std::pow(A, 2) - 1.0 / (3.0 * std::pow(A, 3));

	const double first_term = f1(y, Q2min, alpha, M, lambda2);
	return first_term - prefactor * prefactor_2 * A_term;
}

const double f3(const double y, const double Q2min, const double alpha, const double M, const double lambda2, const double mu) {
	const double prefactor = alpha / (2.0 * M_PI);
	const double value = -0.16666666666666666*(lambda2*(-lambda2 + 4*std::pow(M,2))*(256*std::pow(M,8)*(3*std::pow(lambda2,3) + 22*std::pow(lambda2,2)*Q2min + 30*lambda2*std::pow(Q2min,2) + 12*std::pow(Q2min,3))*std::pow(y,2) - std::pow(lambda2,4)*std::pow(mu,2)*Q2min*(11*std::pow(lambda2,2) + 15*lambda2*Q2min + 6*std::pow(Q2min,2))*(2 + (-2 + y)*y) + 32*lambda2*std::pow(M,6)*(-18*std::pow(lambda2,3)*std::pow(y,2) + 11*std::pow(lambda2,2)*Q2min*(4 + y*(-4 + (-11 + std::pow(mu,2))*y)) + 15*lambda2*std::pow(Q2min,2)*(4 + y*(-4 + (-11 + std::pow(mu,2))*y)) + 6*std::pow(Q2min,3)*(4 + y*(-4 + (-11 + std::pow(mu,2))*y))) + 8*std::pow(lambda2,2)*std::pow(M,4)*(18*std::pow(lambda2,3)*std::pow(y,2) - 45*lambda2*std::pow(Q2min,2)*(4 + y*(-4 + (-3 + std::pow(mu,2))*y)) - 18*std::pow(Q2min,3)*(4 + y*(-4 + (-3 + std::pow(mu,2))*y)) + std::pow(lambda2,2)*Q2min*(-124 + 124*y + 101*std::pow(y,2) + std::pow(mu,2)*(-8 + (8 - 35*y)*y))) + 4*std::pow(lambda2,3)*std::pow(M,2)*(-3*std::pow(lambda2,3)*std::pow(y,2) + std::pow(Q2min,3)*(36 + 3*y*(-12 - y + 3*std::pow(mu,2)*y)) + 3*lambda2*std::pow(Q2min,2)*(28 - 28*y - 3*std::pow(y,2) + 2*std::pow(mu,2)*(1 + y*(-1 + 4*y))) + std::pow(lambda2,2)*Q2min*(52 - 52*y - 9*std::pow(y,2) + 2*std::pow(mu,2)*(7 + y*(-7 + 10*y))))) + 3*Q2min*std::pow(lambda2 + Q2min,3)*(-(std::pow(lambda2,5)*std::pow(-2 + y,2)*log(1 + (4*std::pow(M,2))/Q2min)) + 32*std::pow(M,2)*(2*std::pow(lambda2,4)*y + 128*std::pow(M,8)*std::pow(y,2) + 8*std::pow(lambda2,2)*std::pow(M,4)*y*(4 + 5*y) + 3*std::pow(lambda2,3)*std::pow(M,2)*(4 + std::pow(mu,2)*std::pow(y,2)) + 8*lambda2*std::pow(M,6)*(4 + std::pow(mu,2)*std::pow(y,2)))*log(Q2min/(lambda2 + Q2min)) - 2*std::pow(lambda2,5)*std::pow(mu,2)*std::pow(y,2)*log(lambda2 + Q2min) + 16*lambda2*std::pow(M,2)*(2*std::pow(lambda2,2)*std::pow(M,2)*y*(12 + 5*y) + 16*std::pow(M,6)*y*(4 + 15*y) + std::pow(lambda2,3)*(4 + std::pow(mu,2)*std::pow(y,2)) + 16*lambda2*std::pow(M,4)*(4 + std::pow(mu,2)*std::pow(y,2)))*log((lambda2 + Q2min)/Q2min) + 4*std::pow(lambda2,5)*std::pow(mu,2)*y*log((lambda2 + Q2min)/(4*std::pow(M,2) + Q2min)) + std::pow(lambda2,5)*std::pow(mu,2)*std::pow(y,2)*log(Q2min*(4*std::pow(M,2) + Q2min)) + 4*std::pow(lambda2,5)*std::pow(mu,2)*log((4*std::pow(M,2) + Q2min)/(lambda2 + Q2min))))/(lambda2*std::pow(lambda2 - 4*std::pow(M,2),4)*Q2min*std::pow(lambda2 + Q2min,3)*y);
	return prefactor * value;
}

const double f(const double y, const double Q2min, const EPAParameters params) {
	const double alpha = params.p.alpha;
	const double M = params.p.M;
	const double lambda2 = params.p.lambda2;
	const double mu = params.p.mu;
	switch(params.flux) {
		case Flux::no_mass_magnetic_moment: return f1(y, Q2min, alpha, M, lambda2);
		case Flux::no_magnetic_moment: return f2(y, Q2min, alpha, M, lambda2);
		case Flux::full: return f3(y, Q2min, alpha, M, lambda2, mu);
	}
}

const double sigma_hat(const double s, const double alpha, const double m) {
	const double m2 = m * m;
	const double sqrt_term = sqrt(1.0 - 4.0 * m2 / s);

	const double prefactor = 4.0 * M_PI * alpha * alpha / s;
	const double term_1 = 1.0 + 4.0 * m2 / s - 8.0 * m2 * m2 / (s * s);
	const double log_term = log((1.0 + sqrt_term) / (1.0 - sqrt_term));
	const double term_2 = 1.0 + 4.0 * m2 / s;

	return prefactor * (term_1 * log_term - term_2 * sqrt_term);
}

const double integrand(const double y1, const double y2, const double Q2min1, const double Q2min2, const EPAParameters params) {
	const double s = std::pow(params.p.sqrt_s, 2);
	const double f1 = f(y1, Q2min1, params);
	const double f2 = f(y2, Q2min2, params);
	const double cs_hat = sigma_hat(y1 * y2 * s, params.p.alpha, params.p.m);

	const double s_hat = y1 * y2 * s;
	if (s_hat > std::pow(params.p.sqrt_s_hat_max, 2) || s_hat < std::pow(params.p.sqrt_s_hat_min, 2)) { return 0; }

	const double result = f1 * f2 * cs_hat;
	return result;
}

double gsl_integrand_evaluation(double k[], const size_t dim, void *params_in) {
	const double y1 = k[0];
	const double y2 = k[1];

	if (1.0 - y1 < 1e-15 || 1.0 - y2 < 1e-15) { return 0; }

	const struct EPAParameters *params = (struct EPAParameters *)params_in;

	const double M2 = std::pow(params->p.M, 2);
	const double Q2min1 = M2 * std::pow(y1, 2) / (1.0 - y1);
	const double Q2min2 = M2 * std::pow(y2, 2) / (1.0 - y2);

	return integrand(y1, y2, Q2min1, Q2min2, *params);
}

const Row sigma(EPAParameters params, const bool debug = false) {
	const int dim = 2;
	double integral, error;

	double lower[] = {0, 0};
	double upper[] = {1, 1};

	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_monte_function function;

	function.f = &gsl_integrand_evaluation;
	function.dim = dim;
	function.params = &params;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

	gsl_monte_vegas_params vegas_params;
	gsl_monte_vegas_params_get(state, &vegas_params);
	params.vp.apply_to(&vegas_params);
	gsl_monte_vegas_params_set(state, &vegas_params);

	gsl_monte_vegas_integrate(&function, lower, upper, dim, params.ip.points, r, state, &integral, &error);

	const int iter_max = params.ip.max_iterations;
	int iteration = 0;
	bool iteration_limit_reached = false;

	double best_chi_squared = gsl_monte_vegas_chisq(state);
	double best_integral = integral;
	double best_error = error;

	while (true) {
		iteration++;
		if (iteration > iter_max) {
			iteration_limit_reached = true;
			break;
		}
		gsl_monte_vegas_integrate(&function, lower, upper, dim, params.ip.points, r, state, &integral, &error);

		const double chi_squared = gsl_monte_vegas_chisq(state);
		if (abs(chi_squared - 1.0) < abs(best_chi_squared - 1.0)) {
			best_chi_squared = chi_squared;
			best_integral = integral;
			best_error = error;
		}

		if (abs(gsl_monte_vegas_chisq(state) - 1.0) < params.ip.max_chi_sq_deviation && abs(error / integral) < params.ip.max_relative_error) {
			break;
		}
	}

	const double conversion_factor = Constants::GeV2ToPb;

	Row row;
	row.sqrt_s_hat_min = params.p.sqrt_s_hat_min;
	row.sqrt_s_hat_max = params.p.sqrt_s_hat_max;

	double final_integral = integral;
	double final_error = error;
	double final_chi_squared = gsl_monte_vegas_chisq(state);

	if (iteration_limit_reached) {
		final_integral = best_integral;
		final_error = best_error;
		final_chi_squared = best_chi_squared;
	}

	if (debug) {
		std::cout << "Final integral = " << final_integral << " ± " << final_error << std::endl;
	}

	row.cross_section = conversion_factor * final_integral;
	row.error = conversion_factor * final_error;
	row.chi_squared = final_chi_squared;
	row.iterations = iteration;

 	return row;
}

const void start_integration(EPAParameters params) {
	const IntegrationRoutineParameters ip = params.ip;
	const CommonParameters p = params.p;
	const int step_count = p.step_count();

	std::vector<Row> result;
	result.reserve(step_count);

	int calculated_integrals = 0;

	const bool parallelize = ip.parallelize();

	std::cout << params << std::endl;

	std::ofstream file(ip.filename);

	#pragma omp parallel for if(parallelize) num_threads(ip.thread_count)
	for (int i = 0; i < step_count; i++) {
		const double min = fmax(p.sqrt_s_hat_min, p.sqrt_s_hat_min + i * p.step_size);
		const double max = fmin(p.sqrt_s_hat_max, p.sqrt_s_hat_min + (i + 1) * p.step_size);

		EPAParameters bin_params = params;
		bin_params.p.sqrt_s_hat_min = min;
		bin_params.p.sqrt_s_hat_max = max;

		const Row row = sigma(bin_params, DEBUG);
		#pragma omp critical
		{
			result.push_back(row);
			file << row << std::endl;
			file.flush();

			calculated_integrals++;
			std::cout << "Calculated integral " << calculated_integrals << " / " << step_count << " (index " << i << ", bin [" << min << ", " << max << "])" << std::endl;
		}
	}

	file.close();
}

int main() {
	const size_t calls = 10'000'000;
	const Flux flux = Flux::no_mass_magnetic_moment;

	const double sqrt_s = 100'000;
	const double sqrt_s_hat_lower = 12;
	const double sqrt_s_hat_upper = 70;
	const double step_size = 4;

	const double m = Constants::muon_mass;
	const double M = Constants::proton_mass;
	const double mu = Constants::proton_magnetic_moment;
	const double lambda2 = Constants::proton_dipole_lambda2;
	const double alpha = Constants::alpha_qed;

	const std::string filename = "epa_100000_1.csv";
	const double max_chi_sq_deviation = 0.1;
	const double max_relative_error = 0.1;
	const size_t max_iterations = 200;
	const int seed = 1234567;

	const int thread_count = 8;

	VegasParameters vp = {1.5, 5, GSL_VEGAS_MODE_STRATIFIED, DEBUG ? 0 : -1};
	IntegrationRoutineParameters ip = {calls, 0, max_chi_sq_deviation, max_relative_error, max_iterations, seed, thread_count, filename};
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, sqrt_s_hat_lower, sqrt_s_hat_upper, step_size};
	EPAParameters params = {vp, ip, p, flux};

	start_integration(params);
}

#endif // EPA_H