#pragma once

#include "SimKernel.hxx"
#include "DataManager.hxx"
#include "DataAnalyst.hxx"

class SimManager
{
private:
	std::string name;
	SimKernel kernel;
	DataManager data;
	DataAnalyst afterprocessor;

public:
	SimManager();
	~SimManager();

	void setName(std::string new_name);
	void addVtx(double x, double y, double circulation, double radius, size_t fluidId);
	void buildTimeSample(double t0, double tEnd, size_t nb_steps);
	void setXPeriodicityTo(bool periodic, double period);
	void setYPeriodicityTo(bool periodic, double period);
	void setMethodTo(std::string method);

	void sim(size_t nb_threads);

	std::string getName() const;
	size_t getNbVtx() const;
	size_t getNbSteps() const;
	std::vector<double> getTimeVector() const;

	std::vector<double> getXsAt(size_t step) const;
	std::vector<double> getYsAt(size_t step) const;
	std::vector<double> getUsAt(size_t step) const;
	std::vector<double> getVsAt(size_t step) const;
	std::vector<double> getCirculationsAt(size_t step) const;

	std::vector<double> getVtxXs(size_t vtx) const;
	std::vector<double> getVtxYs(size_t vtx) const;
	std::vector<double> getVtxUs(size_t vtx) const;
	std::vector<double> getVtxVs(size_t vtx) const;
	std::vector<double> getVtxCirculations(size_t vtx) const;

	std::map<size_t, double> computeCompositionAt(double x, double y, size_t step, double radius) const;
	std::vector<double> computeVelocityAt(double x, double y, size_t step, bool x_periodic, double x_period) const;
	double computeVorticityAt(double x, double y, size_t step, double h, bool x_periodic, double x_period) const;

	std::vector<double> computeHamiltonianEvolution(size_t nb_threads) const;

	void saveSim(std::string fileNameWithPath) const;
	void loadSim(std::string fileNameWithPath);
};