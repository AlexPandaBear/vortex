#pragma once

#include <memory>
#include <fstream>
#include "Complex.hxx"

class DataManager
{
private:
	size_t m_nb_timeSamples;
	size_t m_nb_vtx;

	std::unique_ptr<double[]> ptr_time;
	std::unique_ptr<double[]> ptr_x;
	std::unique_ptr<double[]> ptr_y;
	std::unique_ptr<double[]> ptr_u;
	std::unique_ptr<double[]> ptr_v;
	std::unique_ptr<double[]> ptr_circulations;
	std::unique_ptr<double[]> ptr_radiuses;

public:
	DataManager(size_t nb_steps, size_t nb_vtx);
	~DataManager();

	size_t getNbVtx() const;
	size_t getNbSteps() const;
	
	double getTimeAt(size_t timeStep) const;
	double getXAt(size_t timeStep, size_t vortexID) const;
	double getYAt(size_t timeStep, size_t vortexID) const;
	double getUAt(size_t timeStep, size_t vortexID) const;
	double getVAt(size_t timeStep, size_t vortexID) const;
	double getCirculationAt(size_t timeStep, size_t vortexID) const;
	double getRegRadiusAt(size_t timeStep, size_t vortexID) const;

	void storeTimeAt(size_t timeStep, double time);
	void storeXAt(size_t timeStep, size_t vortexID, double x);
	void storeYAt(size_t timeStep, size_t vortexID, double y);
	void storeUAt(size_t timeStep, size_t vortexID, double u);
	void storeVAt(size_t timeStep, size_t vortexID, double v);
	void storeCirculationAt(size_t timeStep, size_t vortexID, double circulation);
	void storeRegRadiusAt(size_t timeStep, size_t vortexID, double radius);

	void saveState(std::string fileNameWithPath) const;
	void loadFile(std::string fileNameWithPath);

	void reset(size_t nb_timeSamples, size_t nb_vtx);
};