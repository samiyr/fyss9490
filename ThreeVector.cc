#ifndef THREE_VECTOR_H
#define THREE_VECTOR_H

/// A Cartesian 3-vector
struct ThreeVector {
	double x;
	double y;
	double z;

	ThreeVector() : x(0), y(0), z(0) {}
	ThreeVector(const double x_in, const double y_in, const double z_in) : x(x_in), y(y_in), z(z_in) {}
	ThreeVector(const double in[3]) : x(in[0]), y(in[1]), z(in[2]) {}
	ThreeVector(const std::vector<double> in) : x(in[0]), y(in[1]), z(in[2]) {}

	/// Constructs a Cartesian 3-vector from the given spherical components.
	/// r >= 0: radial component
	/// -1 <= cos <= 1: cosine of the polar angle 0 <= theta <= pi
	/// 0 <= phi <= 2pi: azimuthal angle
	static const ThreeVector constructSphericalVector(const double r, const double cosine, const double phi);

	/// Returns the square of the norm,
	const double mag_squared() const;
	/// Returns the norm.
	const double mag() const {
		return sqrt(mag_squared());
	}
	friend std::ostream& operator<<(std::ostream& os, ThreeVector const & v) {
		return os << v.x << "\t" << v.y << "\t" << v.z;
	}
};

const double dot(const ThreeVector lhs, const ThreeVector rhs) {
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
const double dot(const ThreeVector lhs) {
	return dot(lhs, lhs);
}
const double ThreeVector::mag_squared() const {
	return dot(*this, *this);
}
const ThreeVector operator+(const ThreeVector& lhs, const ThreeVector& rhs) {
	return ThreeVector(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}
const ThreeVector operator-(const ThreeVector& lhs, const ThreeVector& rhs) {
	return ThreeVector(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

const ThreeVector ThreeVector::constructSphericalVector(const double r, const double cosine, const double phi) {
	// Since 0 <= theta <= pi, the sine is uniquely determined by sqrt(1 - cos^2)
	double sine = sqrt(1 - std::pow(cosine, 2));
	double x = r * cos(phi) * sine;
	double y = r * sin(phi) * sine;
	double z = r * cosine;
	return ThreeVector(x, y, z);
}

#endif // THREE_VECTOR_H