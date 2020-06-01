#pragma once

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

/**
 * A class representing a point vortex for the Vortex Method
 */
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
	/**
	 * The constructor of the class
	 * 
	 * @param x The x coordinate of the vortex
	 * 
	 * @param y The y coordinate of the vortex
	 * 
	 * @param circulation The velocity circulation of the vortex
	 * 
	 * @param regRadius The regularization radius inside which the induced velocity will be restrained in order to avoid singularities
	 * 
	 * @param fluid_id Identification number for the fluid the vortex id part of (Optional - 0 by default)
	 * 
	 * @param X_periodic Boolean indicating if the problem is x-periodic, ie if the vortex represents an infinite row of vortices along the x-axis (Optional - false by default)
	 * 
	 * @param X_period The value of the spatial period along the x-axis if X_periodic is true (Optional - 0. by default - must be positive)
	 */
	Vortex(double x, double y, double circulation, double regRadius, size_t fluid_id = 0, bool X_periodic = false, double X_period = 0.);

	/**
	 * The destructor of the class
	 */
	~Vortex();

	/**
	 * Getter to the x-coordinate of the instance
	 *
	 * @returns The x-coordinate of the instance as a double
	 */
	double getX() const;

	/**
	 * Getter to the y-coordinate of the instance
	 *
	 * @returns The y-coordinate of the instance as a double
	 */
	double getY() const;

	/**
	 * Getter to the circulation of the instance
	 *
	 * @returns The velocity circulation of the instance as a double
	 */
	double getCirculation() const;

	/**
	 * Getter to the regularization radius of the instance
	 *
	 * @returns The regularization radius of the instance as a double
	 */
	double getRegRadius() const;

	/**
	 * Getter to the fluid of the instance
	 *
	 * @returns The identification number of the fluid the instance belongs to as a size-type
	 */
	size_t getFluidId() const;

	/**
	 * Getter to the x-periodicity of the instance
	 *
	 * @returns The x-periodicity of the instance as a boolean
	 */
	bool getXPeriodicity() const;

	/**
	 * Getter to the x-period of the instance
	 *
	 * @returns The period along the x-axis of the instance if x-periodic as a double
	 */
	double getXPeriod() const;

	/**
	 * Setter for the x-coordinate of the instance
	 *
	 * @param newX The new x-coordinate of the instance
	 */
	void setX(double newX);
	
	/**
	 * Setter for the y-coordinate of the instance
	 *
	 * @param newY The new y-coordinate of the instance
	 */
	void setY(double newY);

	/**
	 * Setter for the circulation of the instance
	 *
	 * @param newCirculation The new velocity circulation of the instance
	 */
	void setCirculation(double newCirculation);

	/**
	 * Setter for the regularization radius of the instance
	 *
	 * @param newRegRadius The new regularization radius of the instance
	 */
	void setRegRadius(double newRegRadius);

	/**
	 * Setter for the fluid of the instance
	 *
	 * @param newId The new fluid identification number of the instance
	 */
	void setFluidId(size_t newId);

	/**
	 * Setter for the x-periodicity of the instance
	 *
	 * @param newXPeriodicity The new x-periodicity boolean of the instance
	 */
	void setXPeriodicity(bool newXPeriodicity);

	/**
	 * Setter for the x-period of the instance
	 *
	 * @param newXPeriod The new spatial period along the x-axis of the instance
	 */
	void setXPeriod(double newXPeriod);

	/**
	 * Method returning a printable image of the current state of the instance
	 *
	 * @returns A std::string summing up the current state of the instance
	 */
	std::string toString() const;

	/**
	 * Method that pre-computes useful quantities for the computation of x-periodic induced velocities
	 * 
	 * @warning Calling this method before computing periodic induced velocities is mandatory !
	 */
	void prepareFactors();

	/**
	 * Method computing the velocity induced by this instance along the x-axis at a specific point in space
	 *
	 * @param x The x-coordinate of the point where the velocity will be computed
	 *
	 * @param y The y-coordinate of the point where the velocity will be computed
	 *
	 * @returns The projection on the x-axis of the velocity induced by this instance
	 *
	 * @warning If this instance is x-periodic, calling the prepareFactors method before is required
	 */
	double computeXInducedVelocityAt(double x, double y) const;
	
	/**
	 * Method computing the velocity induced by this instance along the y-axis at a specific point in space
	 *
	 * @param x The x-coordinate of the point where the velocity will be computed
	 *
	 * @param y The y-coordinate of the point where the velocity will be computed
	 *
	 * @returns The projection on the y-axis of the velocity induced by this instance
	 *
	 * @warning If this instance is x-periodic, calling the prepareFactors method before is required
	 */
	double computeYInducedVelocityAt(double x, double y) const;

	/**
	 * Method moving the coordinates of the instance when called
	 *
	 * @param deltaX The amount to add to the x-coordinate
	 *
	 * @param deltaY The amount to add to the y-coordinate
	 */
	void move(double deltaX, double deltaY);
};