#pragma once

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

class Vortex
{
private:
	double m_x;
	double m_y;

	double m_circ;
	double m_rad;

	size_t m_fluid_id;

	bool m_X_periodic;
	double m_X_period;

	double m_global_factor;
	double m_common_factor;

public:
	Vortex(double x, double y, double circulation, double regRadius, size_t fluid_id = 0, bool X_periodic = false, double X_period = 0.);
	~Vortex();

	double getX() const;
	double getY() const;
	double getCirculation() const;
	double getRegRadius() const;
	size_t getFluidId() const;
	bool getXPeriodicity() const;
	double getXPeriod() const;

	void setX(double newX);
	void setY(double newY);
	void setCirculation(double newCirculation);
	void setRegRadius(double newRegRadius);
	void setFluidId(size_t newId);
	void setXPeriodicity(bool newXPeriodicity);
	void setXPeriod(double newXPeriod);

	std::string toString() const;

	void prepareFactors();
	double computeXInducedVelocityAt(double x, double y) const;
	double computeYInducedVelocityAt(double x, double y) const;

	void move(double deltaX, double deltaY);
};