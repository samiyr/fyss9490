#ifndef LORENTZ_VECTOR_H
#define LORENTZ_VECTOR_H

#include <iostream>
#include <cmath>
#include <vector>
#include "ThreeVector.cc"

/// A Lorentz 4-vector
struct LorentzVector {
	double E;
	double px;
	double py;
	double pz;

	LorentzVector() : E(0), px(0), py(0), pz(0) {}
	LorentzVector(const double E_in, const double px_in, const double py_in, const double pz_in) : E(E_in), px(px_in), py(py_in), pz(pz_in) {}
	LorentzVector(const double m_in, const double p_in[3]) : E(sqrt(m_in * m_in + p_in[0] * p_in[0] + p_in[1] * p_in[1] + p_in[2] * p_in[2])), px(p_in[0]), py(p_in[1]), pz(p_in[2]) {}
	LorentzVector(const double m_in, const std::vector<double> p_in) : E(sqrt(m_in * m_in + p_in[0] * p_in[0] + p_in[1] * p_in[1] + p_in[2] * p_in[2])), px(p_in[0]), py(p_in[1]), pz(p_in[2]) {}
	LorentzVector(const double m_in, const ThreeVector spatial) : E(sqrt(m_in * m_in + spatial.mag_squared())), px(spatial.x), py(spatial.y), pz(spatial.z) {}

	friend std::ostream& operator<<(std::ostream& os, LorentzVector const & v) {
		return os << v.E << "\t" << v.px << "\t" << v.py << "\t" << v.pz;
	}

	/// Returns the square of the norm, given by the metric.
	const double mag_squared() const;
	/// Returns the norm, given by the metric.
	const double mag() const {
		return sqrt(mag_squared());
	}
	/// Returns the square of the norm of the spatial components.
	const double momentum_mag_squared() const {
		return spatial().mag_squared();
	}
	/// Returns the norm of the spatial components.
	const double momentum_mag() const {
		return sqrt(momentum_mag_squared());
	}
	/// Checks whether the vector is the zero vector, up to some tolerance.
	const bool is_zero(const double tolerance = 1e-6) const {
		return (abs(E) < tolerance && abs(px) < tolerance && abs(py) < tolerance && abs(pz) < tolerance);
	}
	/// Checks whether any vector component is nan.
	const bool contains_nan() const {
		return std::isnan(E) || std::isnan(px) || std::isnan(py) || std::isnan(pz);
	}
	/// Returns a 3-vector containing the spatial components.
	const ThreeVector spatial() const {
		return ThreeVector(px, py, pz);
	}
	/// Checks whether the vector is on-shell with mass m (p^2 = m^2), up to some tolerance.
	const bool on_shell(const double m, const double tolerance = 1e-6) const {
		return abs(mag_squared() - std::pow(m, 2)) < tolerance;
	}
	/// Returns the squared transverse momentum.
	const double pT2() const {
		return std::pow(px, 2) + std::pow(py, 2);
	}
	/// Returns the transverse momentum.
	const double pT() const {
		return sqrt(pT2());
	}
	/// Returns the pseudorapidity.
	const double eta() const {
		const double pL = pz;
		const double p = momentum_mag();
		return 0.5 * log((p + pL) / (p - pL));
	}
};

const double dot(const LorentzVector lhs, const LorentzVector rhs) {
	return lhs.E * rhs.E - (lhs.px * rhs.px + lhs.py * rhs.py + lhs.pz * rhs.pz);
}
const double dot(const LorentzVector lhs) {
	return dot(lhs, lhs);
}

const double LorentzVector::mag_squared() const {
	return dot(*this, *this);
}

const LorentzVector operator+(const LorentzVector& lhs, const LorentzVector& rhs) {
	return LorentzVector(lhs.E + rhs.E, lhs.px + rhs.px, lhs.py + rhs.py, lhs.pz + rhs.pz);
}
const LorentzVector operator-(const LorentzVector& lhs, const LorentzVector& rhs) {
	return LorentzVector(lhs.E - rhs.E, lhs.px - rhs.px, lhs.py - rhs.py, lhs.pz - rhs.pz);
}
const LorentzVector operator*(const double lhs, const LorentzVector& rhs) {
	return LorentzVector(lhs * rhs.E, lhs * rhs.px, lhs * rhs.py, lhs * rhs.pz);
}

#endif // LORENTZ_VECTOR_H