#pragma once

#include <vector>
#include <thread>
#include <sstream>
#include <iostream>
#include "DataManager.hxx"
#include "DataAnalyst.hxx"
#include "Vortex.hxx"

/**
 * A class implementing the kernel of the Vortex Method
 */
class SimKernel
{
private:
	bool m_X_periodic;
	bool m_Y_periodic;

	double m_X_period;
	double m_Y_period;

	std::string m_method;

	std::vector<Vortex> v_vtx;
	std::vector<double> v_time;

	size_t m_nb_vtx;
	size_t m_nb_steps;

	bool readyToSim();
	
	void integrate(DataManager &dm, size_t step);
	
	void computeEEStep(DataManager &dm, size_t step, size_t firstVtx, size_t lastVtx) const;
	void computeEEStep_multithread(DataManager &dm, size_t step, size_t nb_threads);
	
	void computeRK4Substep(std::vector<Vortex> &workingCopy, std::vector<double> &k_u, std::vector<double> &k_v, size_t step, size_t firstVtx, size_t lastVtx, size_t substep) const;
	void computeRK4Step_multithread(DataManager &dm, size_t step, size_t nb_threads);
	
	void computeEAStep(DataManager &dm, size_t step);
	
	void computeEBStep(DataManager &dm, size_t step);
	
	void computeSVStep(DataManager &dm, size_t step);
	
	void computeSVIStep(DataManager &dm, size_t step);

	void printSimProgression(size_t step) const;

public:
	/**
	 * The constructor of the class
	 *
	 * @param x_periodic A boolean telling if the simulation is x-periodic
	 *
	 * @param y_periodic A boolean telling if the simulation is y-periodic \todo
	 *
	 * @param x_period The value of the eventual spatial period along the x-axis
	 *
	 * @param y_period The value of the eventual spatial period along the y-axis
	 *
	 * @param method The name of the time stepping integration method (Possible values : 'euler', 'rk4', 'eulerA', 'eulerB', 'sv', 'svi')
	 */
	SimKernel(bool x_periodic, bool y_periodic, double x_period, double y_period, std::string method);
	
	/**
	 * The destructor of the class	 
	 */
	~SimKernel();

	/**
	 * Setter for the x-periodicity of the instance
	 * 
	 * @param periodic A boolean telling if computation should be periodic along the x-axis
	 *
	 * @param period The value of the eventual spatial period
	 */
	void setXPeriodicityTo(bool periodic, double period);

	/**
	 * Setter for the y-periodicity of the instance
	 * 
	 * @param periodic A boolean telling if computation should be periodic along the y-axis
	 *
	 * @param period The value of the eventual spatial period
	 */
	void setYPeriodicityTo(bool periodic, double period);

	/**
	 * Setter for the time stepping integration method to use
	 *
	 * @param method The name of the method (Possible values : 'euler', 'rk4', 'eulerA', 'eulerB', 'sv', 'svi')
	 */
	void setMethodTo(std::string method);
	
	/**
	 * Method adding a vortex to the list of vortex that will be used for computation
	 *
	 * @param v The vortex to add to the list
	 */
	void addVortex(Vortex &v);

	/**
	 * Method building the list of the time-values used for computation
	 *
	 * @param t0 The time at which the simulation will start
	 *
	 * @param tEnd The time at which the simulation will end (must be greater than t0)
	 *
	 * @param steps The number of time steps to simulate from t0 to tEnd (must be greater than 1)
	 */
	void buildTimeSample(double t0, double tEnd, size_t steps);

	/**
	 * Method simulating the flow with all the parameters previously defined
	 *
	 * @param dm The DataManager object in which the simulation will be stored
	 *
	 * @param nb_threads The number of threads to use to speed-up computation (must be greater than 1)
	 *
	 * @warning All the vortices must be defined and added to the simulation (with the addVortex method), and the time sample must be defined (with the buildTimeSample method) before calling this method
	 */
	void sim(DataManager &dm, size_t nb_threads);
};