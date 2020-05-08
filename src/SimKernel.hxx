#pragma once

#include <vector>
#include <thread>
#include <sstream>
#include <iostream>
#include "DataManager.hxx"
#include "DataAnalyst.hxx"
#include "Vortex.hxx"

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
	void computeSVStep(DataManager &dm, size_t step);

public:
	SimKernel(bool x_periodic, bool y_periodic, double x_period, double y_period, std::string method);
	~SimKernel();

	void setXPeriodicityTo(bool periodic, double period);
	void setYPeriodicityTo(bool periodic, double period);
	void setMethodTo(std::string method);
	
	void addVortex(Vortex &v);
	void buildTimeSample(double t0, double tEnd, size_t steps);

	void sim(DataManager &dm, size_t nb_threads);
};