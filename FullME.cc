#ifndef FULL_ME_H
#define FULL_ME_H

#include "Contraction.cc"
#include "Utility.cc"
#include <vector>
#include "LorentzVector.cc"
#include <gsl/gsl_monte_vegas.h>
#include <cassert>
#include <cmath>
#include <optional>
#include <iomanip>
#include <fstream>
#include <limits>

#define DEBUG false

/// A quick hack for including experimental cuts for ATLAS.
/// 0: no ATLAS cuts
/// 1: invariant mass bin [12, 17], pT > 6 GeV and |eta| < 2.4
/// 2: invariant mass bin [17, 22], pT > 6 GeV and |eta| < 2.4
/// 3: invariant mass bin [22, 30], pT > 6 GeV and |eta| < 2.4
/// 4: invariant mass bin [30, 70], pT > 10 GeV and |eta| < 2.4
#define EXPERIMENTAL_CUT_SET 0

#if DEBUG
double Q1_squared_min = std::numeric_limits<double>::infinity();
double Q1_squared_min_corresponding_Q2 = 0;
double Q2_squared_min = std::numeric_limits<double>::infinity();
double Q2_squared_min_corresponding_Q1 = 0;
bool sqrt_s_hat_range_found = false;
double sqrt_s_hat_min = std::numeric_limits<double>::infinity();
double sqrt_s_hat_max = 0;

double r_p1_iter_min = std::numeric_limits<double>::infinity();
double r_p1_iter_max = 0;

double cos_p1_iter_min = std::numeric_limits<double>::infinity();
double cos_p1_iter_max = -std::numeric_limits<double>::infinity();

double phi_p1_iter_min = std::numeric_limits<double>::infinity();
double phi_p1_iter_max = -std::numeric_limits<double>::infinity();

double Q1_squared_iter_min = std::numeric_limits<double>::infinity();
double Q1_squared_iter_max = 0;

double r_p2_iter_min = std::numeric_limits<double>::infinity();
double r_p2_iter_max = 0;

double cos_p2_iter_min = std::numeric_limits<double>::infinity();
double cos_p2_iter_max = -std::numeric_limits<double>::infinity();

double Q2_squared_iter_min = std::numeric_limits<double>::infinity();
double Q2_squared_iter_max = 0;

double phi_p2_iter_min = std::numeric_limits<double>::infinity();
double phi_p2_iter_max = -std::numeric_limits<double>::infinity();

double cos_k1_iter_min = std::numeric_limits<double>::infinity();
double cos_k1_iter_max = -std::numeric_limits<double>::infinity();

double phi_k1_iter_min = std::numeric_limits<double>::infinity();
double phi_k1_iter_max = -std::numeric_limits<double>::infinity();

double sqrt_s_hat_iter_min = std::numeric_limits<double>::infinity();
double sqrt_s_hat_iter_max = 0;
#endif

struct PhaseSpacePoint {
	LorentzVector p1_out;
	LorentzVector p2_out;
	LorentzVector k1_out;
	LorentzVector k2_out;

	double r_p1;
	double cos_p1;
	double sin_p1;
	double phi_p1;

	double r_p2;
	double cos_p2;
	double sin_p2;
	double phi_p2;

	double r_k1;
	double cos_k1;
	double sin_k1;
	double phi_k1;

	double epsilon;
	double u;
	double rho2;
};

struct PhaseSpaceParameters {
	double s;
	double pz;
	double pz_squared;

	double m2;
	double M2;

	LorentzVector p1;
	LorentzVector p2;
	LorentzVector in_sum;
	double p1p2;
	
	PhaseSpaceParameters(CommonParameters p) {
		m2 = std::pow(p.m, 2);
		M2 = std::pow(p.M, 2);
		s = std::pow(p.sqrt_s, 2);
		pz_squared = s / 4.0 - M2;
		pz = sqrt(pz_squared);

		p1 = LorentzVector(p.sqrt_s / 2.0, 0, 0, pz);
		p2 = LorentzVector(p.sqrt_s / 2.0, 0, 0, -pz);
		in_sum = LorentzVector(p.sqrt_s, 0, 0, 0);
		p1p2 = 0.5 * (s - 2 * M2);
	}

	friend std::ostream& operator<<(std::ostream& os, PhaseSpaceParameters const &p) {
		os << "==========================" << std::endl;
		os << "=== PHASE SPACE PARAMS ===" << std::endl;
		os << "==========================" << std::endl;
		os << "s = " << p.s << std::endl;
		os << "pz = " << p.pz << std::endl;
		os << "pz_squared = " << p.pz_squared << std::endl;
		os << "m2 = " << p.m2 << std::endl;
		os << "M2 = " << p.M2 << std::endl;
		os << "p1 = " << p.p1 << std::endl;
		os << "p2 = " << p.p2 << std::endl;
		os << "in_sum = " << p.in_sum << std::endl;
		os << "p1p2 = " << p.p1p2 << std::endl;
		os << "==========================" << std::endl << std::endl;
		return os;
	}
};

struct FullMEParameters {
	VegasParameters phase1;
	VegasParameters phase2;
	IntegrationRoutineParameters ip;
	CommonParameters p;
	PhaseSpaceParameters php;

	friend std::ostream& operator<<(std::ostream& os, FullMEParameters const &p) {
		return os << p.phase1 << p.phase2 << p.ip << p.p << p.php;
	}
};

const bool r1_check(const double r1, const FullMEParameters p, const double epsilon, const double u, const double rho2, const double tolerance = 1e-6) {
	const double value = epsilon - sqrt(p.php.m2 + std::pow(r1, 2)) - sqrt(p.php.m2 + rho2 + 2 * u * r1 + std::pow(r1, 2));
	return abs(value) < tolerance;
}

const std::vector<std::optional<PhaseSpacePoint>> constructPhaseSpacePoints(const double x[], const FullMEParameters p) {
	std::optional<PhaseSpacePoint> point1 = std::nullopt;
	std::optional<PhaseSpacePoint> point2 = std::nullopt;

	const double r_p1 = x[0];
	const double cos_p1 = x[1];
	const double phi_p1 = x[2];

	const double r_p2 = x[3];
	const double cos_p2 = x[4];
	const double phi_p2 = x[5];
	
	const double cos_k1 = x[6];
	const double phi_k1 = x[7];

	const double sin_p1 = sqrt(1 - std::pow(cos_p1, 2));
	const double sin_p2 = sqrt(1 - std::pow(cos_p2, 2));
	const double sin_k1 = sqrt(1 - std::pow(cos_k1, 2));

	const ThreeVector p1_vector = ThreeVector::constructSphericalVector(r_p1, cos_p1, phi_p1);
	const LorentzVector p1(p.p.M, p1_vector);

	const ThreeVector p2_vector = ThreeVector::constructSphericalVector(r_p2, cos_p2, phi_p2);
	const LorentzVector p2(p.p.M, p2_vector);

	const double u(cos_k1 * (r_p1 * cos_p1 + r_p2 * cos_p2) + sin_k1 * (r_p1 * cos(phi_k1 - phi_p1) * sin_p1 + r_p2 * cos(phi_k1 - phi_p2) * sin_p2));
	const double epsilon = p.p.sqrt_s - p1.E - p2.E;
	const double rho2 = (p1_vector + p2_vector).mag_squared();

	assert(epsilon > 0);
	assert(rho2 > 0);

	const double epsilon2 = std::pow(epsilon, 2);
	const double term_1 = epsilon2 - rho2;
	const double term_2 = epsilon2 - std::pow(u, 2);
	const double sqrt_term = epsilon * sqrt(std::pow(term_1, 2) - 4.0 * p.php.m2 * term_2);

	const double r_k1_plus = (- u * term_1 + sqrt_term) / (2.0 * term_2);
	const double r_k1_minus = (- u * term_1 - sqrt_term) / (2.0 * term_2);

	const bool r_k1_plus_is_solution = r1_check(r_k1_plus, p, epsilon, u, rho2);
	const bool r_k1_minus_is_solution = r1_check(r_k1_minus, p, epsilon, u, rho2);

	const bool r_k1_plus_check = r_k1_plus > 0 && r_k1_plus_is_solution;
	const bool r_k1_minus_check = r_k1_minus > 0 && r_k1_minus_is_solution;

	if (r_k1_plus_check) {
		const ThreeVector k1_vector_plus = ThreeVector::constructSphericalVector(r_k1_plus, cos_k1, phi_k1);
		const LorentzVector k1_plus(p.p.m, k1_vector_plus);

		const LorentzVector k2_plus(p.p.m, p.php.in_sum.spatial() - p1_vector - p2_vector - k1_vector_plus);

		point1 = {p1, p2, k1_plus, k2_plus, r_p1, cos_p1, sin_p1, phi_p1, r_p2, cos_p2, sin_p2, phi_p2, r_k1_plus, cos_k1, sin_k1, phi_k1, epsilon, u, rho2};
	}
	
	if (r_k1_minus_check) {
		const ThreeVector k1_vector_minus = ThreeVector::constructSphericalVector(r_k1_minus, cos_k1, phi_k1);
		const LorentzVector k1_minus(p.p.m, k1_vector_minus);

		const LorentzVector k2_minus(p.p.m, p.php.in_sum.spatial() - p1_vector - p2_vector - k1_vector_minus);

		point2 = {p1, p2, k1_minus, k2_minus, r_p1, cos_p1, sin_p1, phi_p1, r_p2, cos_p2, sin_p2, phi_p2, r_k1_minus, cos_k1, sin_k1, phi_k1, epsilon, u, rho2};
	}

	return {point1, point2};
}

const double integrand_product_term(const PhaseSpacePoint point, const FullMEParameters p) {
	const LorentzVector p1_out = point.p1_out;
	const LorentzVector p2_out = point.p2_out;
	const LorentzVector k1_out = point.k1_out;
	const LorentzVector k2_out = point.k2_out;

	const LorentzVector q1 = p.php.p1 - p1_out;
	const LorentzVector q2 = p.php.p2 - p2_out;

	const double Q1_squared = -dot(q1, q1);
	const double Q2_squared = -dot(q2, q2);

	assert(Q1_squared > 0);
	assert(Q2_squared > 0);

	assert(!p1_out.contains_nan());
	assert(!p2_out.contains_nan());
	assert(!k1_out.contains_nan());
	assert(!k2_out.contains_nan());
	
	assert(p1_out.E > 0);
	assert(p2_out.E > 0);
	assert(k1_out.E > 0);
	assert(k2_out.E > 0);

	const double s_hat = dot(k1_out + k2_out);
	const double sqrt_s_hat = sqrt(s_hat);
	const double r_k1 = point.r_k1;

	#if DEBUG
	if (sqrt_s_hat < sqrt_s_hat_min) {
		sqrt_s_hat_min = sqrt_s_hat;
	}

	if (sqrt_s_hat > sqrt_s_hat_max) {
		sqrt_s_hat_max = sqrt_s_hat;
	}
	#endif

	if (sqrt_s_hat > p.p.sqrt_s_hat_max || sqrt_s_hat < p.p.sqrt_s_hat_min) { return 0; }

	#if EXPERIMENTAL_CUT_SET > 0
	const double eta_1 = k1_out.eta();
	const double eta_2 = k2_out.eta();

	if (abs(eta_1) >= 2.4 || abs(eta_2) >= 2.4) { return 0; }

	const double pT_1 = k1_out.pT();
	const double pT_2 = k2_out.pT();	
	#endif

	#if EXPERIMENTAL_CUT_SET == 1 || EXPERIMENTAL_CUT_SET == 2 || EXPERIMENTAL_CUT_SET == 3
	if (pT_1 <= 6 || pT_2 <= 6) { return 0; }
	#elif EXPERIMENTAL_CUT_SET == 4
	if (pT_1 <= 10 || pT_2 <= 10) { return 0; }
	#endif

	#if DEBUG
	if (sqrt_s_hat_min < p.p.sqrt_s_hat_min && sqrt_s_hat_max > p.p.sqrt_s_hat_max) {
		sqrt_s_hat_range_found = true;
	}
	#endif

	assert(s_hat >= 4.0 * p.php.m2);
	assert(p.php.p1.on_shell(p.p.M));
	assert(p.php.p2.on_shell(p.p.M));
	assert(p1_out.on_shell(p.p.M));
	assert(p2_out.on_shell(p.p.M));
	assert(k1_out.on_shell(p.p.m));
	assert(k2_out.on_shell(p.p.m));
	
	#if DEBUG
	const LorentzVector total = p.php.p1 + p.php.p2 - p1_out - p2_out - k1_out - k2_out;
	assert(total.is_zero());
	#endif
	
	#if DEBUG
	const double r_p1 = point.r_p1;
	const double cos_p1 = point.cos_p1;
	const double phi_p1 = point.phi_p1;

	const double r_p2 = point.r_p2;
	const double cos_p2 = point.cos_p2;
	const double phi_p2 = point.phi_p2;

	const double cos_k1 = point.cos_k1;
	const double phi_k1 = point.phi_k1;

	if (Q1_squared < Q1_squared_min) {
		Q1_squared_min = Q1_squared;
		Q1_squared_min_corresponding_Q2 = Q2_squared;
	}
	if (Q2_squared < Q2_squared_min) {
		Q2_squared_min = Q2_squared;
		Q2_squared_min_corresponding_Q1 = Q1_squared;
	}

	if (r_p1 < r_p1_iter_min) {
		r_p1_iter_min = r_p1;
	}
	if (r_p1 > r_p1_iter_max) {
		r_p1_iter_max = r_p1;
	}

	if (cos_p1 < cos_p1_iter_min) {
		cos_p1_iter_min = cos_p1;
	}
	if (cos_p1 > cos_p1_iter_max) {
		cos_p1_iter_max = cos_p1;
	}

	if (r_p2 < r_p2_iter_min) {
		r_p2_iter_min = r_p2;
	}
	if (r_p2 > r_p2_iter_max) {
		r_p2_iter_max = r_p2;
	}

	if (cos_p2 < cos_p2_iter_min) {
		cos_p2_iter_min = cos_p2;
	}
	if (cos_p2 > cos_p2_iter_max) {
		cos_p2_iter_max = cos_p2;
	}

	if (Q1_squared < Q1_squared_iter_min) {
		Q1_squared_iter_min = Q1_squared;
	}
	if (Q1_squared > Q1_squared_iter_max) {
		Q1_squared_iter_max = Q1_squared;
	}
	
	if (Q2_squared < Q2_squared_iter_min) {
		Q2_squared_iter_min = Q2_squared;
	}
	if (Q2_squared > Q2_squared_iter_max) {
		Q2_squared_iter_max = Q2_squared;
	}

	if (phi_p1 < phi_p1_iter_min) {
		phi_p1_iter_min = phi_p1;
	}
	if (phi_p1 > phi_p1_iter_max) {
		phi_p1_iter_max = phi_p1;
	}

	if (phi_p2 < phi_p2_iter_min) {
		phi_p2_iter_min = phi_p2;
	}
	if (phi_p2 > phi_p2_iter_max) {
		phi_p2_iter_max = phi_p2;
	}

	if (cos_k1 < cos_k1_iter_min) {
		cos_k1_iter_min = cos_k1;
	}
	if (cos_k1 > cos_k1_iter_max) {
		cos_k1_iter_max = cos_k1;
	}

	if (phi_k1 < phi_k1_iter_min) {
		phi_k1_iter_min = phi_k1;
	}
	if (phi_k1 > phi_k1_iter_max) {
		phi_k1_iter_max = phi_k1;
	}

	if (sqrt_s_hat < sqrt_s_hat_iter_min) {
		sqrt_s_hat_iter_min = sqrt_s_hat;
	}
	if (phi_k1 > sqrt_s_hat_iter_max) {
		sqrt_s_hat_iter_max = sqrt_s_hat;
	}
	#endif

	const double term_1 = 1.0 + Q1_squared / p.p.lambda2;
	const double GE21 = 1.0 / std::pow(term_1, 4);

	const double term_2 = 1.0 + Q2_squared / p.p.lambda2;
	const double GE22 = 1.0 / std::pow(term_2, 4);

	const double y1_inv = p.php.p1p2 / dot(q1, p.php.p2);
	const double y2_inv = p.php.p1p2 / dot(q2, p.php.p1);

	const LorentzVector P1 = p.php.p1 - y1_inv * q1;
	const LorentzVector P2 = p.php.p2 - y2_inv * q2;

	const double k1P1 = dot(k1_out, P1);
	const double k1P2 = dot(k1_out, P2);
	const double k1q2 = dot(k1_out, q2);
	const double P1q2 = dot(P1, q2);
	const double q2q2 = dot(q2, q2);
	const double P1P2 = dot(P1, P2);
	const double k1q1 = dot(k1_out, q1);
	const double P2q1 = dot(P2, q1);
	const double q1q2 = dot(q1, q2);
	const double q1q1 = dot(q1, q1);
	const double P1P1 = dot(P1, P1);
	const double P2P2 = dot(P2, P2);
	const double P1q1 = dot(P1, q1);
	const double P2q2 = dot(P2, q2);

	const double c = contraction(k1P1, k1P2, k1q2, P1q2, q2q2, P1P2, k1q1, P2q1, q1q2, q1q1, P1P1, P2P2, P1q1, P2q2, p.p.m, p.p.M, p.p.mu, GE21, GE22);

	const double r = r_k1;
	const double r2 = std::pow(r_k1, 2);
	const double epsilon = point.epsilon;
	const double u = point.u;
	const double rho2 = point.rho2;
	const double f = abs((r * epsilon) / k1_out.E + u) / sqrt(p.php.m2 + rho2 + 2 * u * r + r2);

	const double result = (c * r2) / (k1_out.E * k2_out.E * f * std::pow(Q1_squared, 2) * std::pow(Q2_squared, 2));

	return result;
}

const double integrand(const std::vector<PhaseSpacePoint> points, const FullMEParameters p) {
	if (points.empty()) { return 0; }

	const LorentzVector p1_out = points[0].p1_out;
	const LorentzVector p2_out = points[0].p2_out;

	double sum = 0.0;
	for (auto point : points) {
		const double value = integrand_product_term(point, p);
		sum += value;
	}

	const double result = (p1_out.momentum_mag_squared() * p2_out.momentum_mag_squared() * sum) / (p1_out.E * p2_out.E);
	return result;
}

double gsl_integrand_evaluation(double k[], const size_t dim, void *params) {
	FullMEParameters *p = (FullMEParameters *)params;
	const std::vector<std::optional<PhaseSpacePoint>> optional_points = constructPhaseSpacePoints(k, *p);
	
	std::vector<PhaseSpacePoint> points;
	for (auto optional : optional_points) {
		if (optional) { points.push_back(*optional); }
	}

	const double result = integrand(points, *p);
	return result;
}

const Row integrated_amplitude(FullMEParameters params, const int i, const bool debug = false) {
	const int dim = 8;

	double integral, error;

	const double l = 0.0;
	const double u = params.p.sqrt_s / 2.0;

	const double tau = 2.0 * M_PI;
	double lower[] = {l, -1, 0, l, -1, 0, -1, 0};
	double upper[] = {u, 1, tau, u, 1, tau, 1, tau};

	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_monte_function function;

	function.f = &gsl_integrand_evaluation;
	function.dim = dim;
	function.params = &params;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, params.ip.seed);

	gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

	gsl_monte_vegas_params vegas_params;
	gsl_monte_vegas_params_get(state, &vegas_params);
	params.phase1.apply_to(&vegas_params);
	gsl_monte_vegas_params_set(state, &vegas_params);

	gsl_monte_vegas_integrate(&function, lower, upper, dim, params.ip.points, r, state, &integral, &error);

	size_t iteration = 0;
	bool iteration_limit_reached = false;

	double best_chi_squared = gsl_monte_vegas_chisq(state);
	double best_integral = integral;
 	double best_error = error;

 	#if DEBUG
	std::cout << std::endl << "Q1_squared_min = " << Q1_squared_min << " with Q2_squared = " << Q1_squared_min_corresponding_Q2 << std::endl;
	std::cout << "Q2_squared_min = " << Q2_squared_min << " with Q1_squared = " << Q2_squared_min_corresponding_Q1 << std::endl << std::endl;
	#endif

	bool allowed_region_found = false;

	while (true) {
		iteration++;
		if (iteration > params.ip.max_iterations) {
			iteration_limit_reached = true;
			break;
		}

		#if DEBUG
		std::cout << "Iteration " << iteration << " (index " << i << "):" << std::endl;
		#else
		std::cout << "Iteration " << iteration << " (index " << i << ")" << std::endl;
		#endif

		gsl_monte_vegas_integrate(&function, lower, upper, dim, params.ip.points, r, state, &integral, &error);

		const double chi_squared = gsl_monte_vegas_chisq(state);
		if (abs(chi_squared - 1.0) < abs(best_chi_squared - 1.0)) {
			best_chi_squared = chi_squared;
			best_integral = integral;
			best_error = error;
		}

		if (!allowed_region_found && integral > params.ip.integration_phase_cutoff) {
			allowed_region_found = true;

			gsl_monte_vegas_params vegas_params;
			gsl_monte_vegas_params_get(state, &vegas_params);
			params.phase2.apply_to(&vegas_params);
			gsl_monte_vegas_params_set(state, &vegas_params);
		}

		#if DEBUG_INTEGRATION_REGION
		std::cout << std::endl << "global Q1_squared_min = " << Q1_squared_min << " with Q2_squared = " << Q1_squared_min_corresponding_Q2 << std::endl;
		std::cout << "global Q2_squared_min = " << Q2_squared_min << " with Q1_squared = " << Q2_squared_min_corresponding_Q1 << std::endl;
		std::cout << "global sqrt_s_hat = " << sqrt_s_hat_min << " ... " << sqrt_s_hat_max << std::endl << std::endl;
	
		std::cout << "local r_p1 = " << r_p1_iter_min << " ... " << r_p1_iter_max << std::endl;
		std::cout << "local cos_p1 = " << cos_p1_iter_min << " ... " << cos_p1_iter_max << std::endl;
		std::cout << "local phi_p1 = " << phi_p1_iter_min << " ... " << phi_p1_iter_max << std::endl;
		std::cout << "local Q1_squared = " << Q1_squared_iter_min << " ... " << Q1_squared_iter_max << std::endl;
		std::cout << "local r_p2 = " << r_p2_iter_min << " ... " << r_p2_iter_max << std::endl;
		std::cout << "local cos_p2 = " << cos_p2_iter_min << " ... " << cos_p2_iter_max << std::endl;
		std::cout << "local phi_p2 = " << phi_p2_iter_min << " ... " << phi_p2_iter_max << std::endl;
		std::cout << "local Q2_squared = " << Q2_squared_iter_min << " ... " << Q2_squared_iter_max << std::endl;
		std::cout << "local cos_k1 = " << cos_k1_iter_min << " ... " << cos_k1_iter_max << std::endl;
		std::cout << "local phi_k1 = " << phi_k1_iter_min << " ... " << phi_k1_iter_max << std::endl << std::endl;

		r_p1_iter_min = std::numeric_limits<double>::infinity();
		r_p1_iter_max = 0;

		cos_p1_iter_min = std::numeric_limits<double>::infinity();
		cos_p1_iter_max = -std::numeric_limits<double>::infinity();

		r_p2_iter_min = std::numeric_limits<double>::infinity();
		r_p2_iter_max = 0;

		cos_p2_iter_min = std::numeric_limits<double>::infinity();
		cos_p2_iter_max = -std::numeric_limits<double>::infinity();

		Q1_squared_iter_min = std::numeric_limits<double>::infinity();
		Q1_squared_iter_max = 0;

		Q2_squared_iter_min = std::numeric_limits<double>::infinity();
		Q2_squared_iter_max = 0;

		phi_p1_iter_min = std::numeric_limits<double>::infinity();
		phi_p1_iter_max = -std::numeric_limits<double>::infinity();

		phi_p2_iter_min = std::numeric_limits<double>::infinity();
		phi_p2_iter_max = -std::numeric_limits<double>::infinity();

		cos_k1_iter_min = std::numeric_limits<double>::infinity();
		cos_k1_iter_max = -std::numeric_limits<double>::infinity();

		phi_k1_iter_min = std::numeric_limits<double>::infinity();
		phi_k1_iter_max = -std::numeric_limits<double>::infinity();

		sqrt_s_hat_iter_min = std::numeric_limits<double>::infinity();
		sqrt_s_hat_iter_max = 0;
		#endif

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

	const double prefactor = std::pow(params.p.alpha, 4) / (64.0 * std::pow(M_PI, 4));
	const double flux_factor = sqrt(std::pow(params.php.p1p2, 2) - std::pow(params.p.M, 4));
	const double total_factor = prefactor * conversion_factor / flux_factor;

	row.cross_section = total_factor * final_integral;
	row.error = total_factor * final_error;
	row.chi_squared = final_chi_squared;
	row.iterations = iteration;

	return row;
}

void start_integration(FullMEParameters params) {
	const int step_count = params.p.step_count();

	std::vector<Row> result;
	result.reserve(step_count);

	int calculated_integrals = 0;

	const bool parallelize = params.ip.parallelize();

	std::cout << params << std::endl;

	std::ofstream file(params.ip.filename);

	#pragma omp parallel for if(parallelize) num_threads(params.ip.thread_count)
	for (int i = 0; i < step_count; i++) {
		const double min = fmax(params.p.sqrt_s_hat_min, params.p.sqrt_s_hat_min + i * params.p.step_size);
		const double max = fmin(params.p.sqrt_s_hat_max, params.p.sqrt_s_hat_min + (i + 1) * params.p.step_size);
		
		FullMEParameters bin_params = params;
		bin_params.p.sqrt_s_hat_min = min;
		bin_params.p.sqrt_s_hat_max = max;

		const Row row = integrated_amplitude(bin_params, i, DEBUG);
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
	const size_t calls = 2'000'000'000;

	const double sqrt_s = 13'000;
	const double sqrt_s_hat_lower = 12;
	const double sqrt_s_hat_upper = 70;
	const double step_size = 4;

	const double m = Constants::muon_mass;
	const double M = Constants::proton_mass;
	const double mu = Constants::proton_magnetic_moment;
	const double lambda2 = Constants::proton_dipole_lambda2;
	const double alpha = Constants::alpha_qed;

	const std::string filename = "full_me_13000.csv";
	const double max_chi_sq_deviation = 0.1;
	const double max_relative_error = 0.1;
	const double integration_phase_cutoff = 1;
	const size_t max_iterations = 200;
	const int seed = 1234567;

	const int thread_count = 16;

	VegasParameters phase1 = {2.0, 2, GSL_VEGAS_MODE_STRATIFIED, DEBUG ? 0 : -1};
	VegasParameters phase2 = {0.1, 20, GSL_VEGAS_MODE_STRATIFIED, DEBUG ? 0 : -1};
	IntegrationRoutineParameters ip = {calls, integration_phase_cutoff, max_chi_sq_deviation, max_relative_error, max_iterations, seed, thread_count, filename};
	#if EXPERIMENTAL_CUT_SET == 0
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, sqrt_s_hat_lower, sqrt_s_hat_upper, step_size};
	#elif EXPERIMENTAL_CUT_SET == 1
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, 12, 17, 5};
	#elif EXPERIMENTAL_CUT_SET == 2
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, 17, 22, 5};
	#elif EXPERIMENTAL_CUT_SET == 3
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, 22, 30, 8};
	#elif EXPERIMENTAL_CUT_SET == 4
	CommonParameters p = {sqrt_s, m, M, alpha, lambda2, mu, 30, 70, 40};
	#endif
	PhaseSpaceParameters php(p);

	FullMEParameters params = {phase1, phase2, ip, p, php};

	start_integration(params);
}

#endif // FULL_ME_H