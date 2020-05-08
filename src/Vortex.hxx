#pragma once

#include <cmath>
#include <string>
#include <sstream>

class Vortex
{
private:
	double m_x;
	double m_y;

	double m_circ;
	double m_rad;

public:
	Vortex(double x, double y, double circulation, double regRadius);
	~Vortex();

	double getX() const;
	double getY() const;
	double getCirculation() const;
	double getRegRadius() const;

	void setX(double newX);
	void setY(double newY);
	void setCirculation(double newCirculation);
	void setRegRadius(double newRegRadius);

	std::string toString() const;

	double computeXInducedVelocityAt(double x, double y) const;
	double computeYInducedVelocityAt(double x, double y) const;
	static double computeXInducedVelocityAt(double x, double y, double vtx_x, double vtx_y, double vtx_circ, double vtx_rad);
	static double computeYInducedVelocityAt(double x, double y, double vtx_x, double vtx_y, double vtx_circ, double vtx_rad);

	void move(double deltaX, double deltaY);
};