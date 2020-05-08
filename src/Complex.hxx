#pragma once

#include <cmath>
#include <sstream>
#include <exception>
#include <string>

class Complex
{
private:
	double m_real;
	double m_imag;

public:
	Complex(double x = 0., double y = 0., std::string mode = "cartesian");
	~Complex();

	double getReal() const;
	double getImag() const;
	double getMag() const;
	double getMagSquared() const;
	double getArg() const;

	void setReal(double x);
	void setImag(double y);

	std::string toString() const;
	std::string toFileFormat() const;

	bool operator==(Complex const &z) const;
	bool operator!=(Complex const &z) const;

	Complex inverse() const;
	Complex conj() const;

	void operator+=(Complex const &z);
	void operator-=(Complex const &z);
	void operator*=(Complex const &z);
	void operator/=(Complex const &z);

	Complex operator+(Complex const &z) const;
	Complex operator-(Complex const &z) const;
	Complex operator*(Complex const &z) const;
	Complex operator/(Complex const &z) const;

	Complex operator-() const;
	Complex operator*(double x) const;

	static Complex exp(Complex const &z);
	static Complex ln(Complex const &z);

	static double dotProduct(Complex const &z1, Complex const &z2);
	static double crossProduct(Complex const &z1, Complex const &z2);
};