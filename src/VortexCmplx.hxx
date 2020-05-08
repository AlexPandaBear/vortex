#pragma once

#include "Complex.hxx"

class VortexCmplx
{
private:
	Complex m_pos;

	double m_circ;
	double m_rad;

public:
	VortexCmplx(Complex const &position, double circulation, double regRadius);
	VortexCmplx(double x, double y, double circulation, double regRadius, std::string mode = "cartesian");
	~VortexCmplx();

	Complex getPosition() const;
	double getCirculation() const;
	double getRegRadius() const;

	void setPosition(Complex newPosition);
	void setCirculation(double newCirculation);
	void setRegRadius(double newRegRadius);

	std::string toString() const;

	Complex computeInducedVelocityAt(Complex position) const;
	static Complex computeInducedVelocityAt(Complex position, Complex vtx_pos, double vtx_circ, double vtx_rad);

	void move(Complex deltaPosition);
};