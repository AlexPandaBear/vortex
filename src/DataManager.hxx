#pragma once

#include <memory>
#include <fstream>
#include "Complex.hxx"

/**
 * A class managing simulation data storage
 */
class DataManager
{
private:
	size_t m_nb_timeSamples;
	size_t m_nb_vtx;

	std::unique_ptr<double[]> ptr_time;
	std::unique_ptr<size_t[]> ptr_fluid_id;
	std::unique_ptr<double[]> ptr_x;
	std::unique_ptr<double[]> ptr_y;
	std::unique_ptr<double[]> ptr_u;
	std::unique_ptr<double[]> ptr_v;
	std::unique_ptr<double[]> ptr_circulations;
	std::unique_ptr<double[]> ptr_radiuses;

public:
	/**
	 * The constructor of the class
	 *
	 * @param nb_steps The number of time steps of the simulation to save
	 *
	 * @param nb_vtx The number of vortices used in the simulation to save
	 */
	DataManager(size_t nb_steps, size_t nb_vtx);

	/**
	 * The destructor of the class
	 */
	~DataManager();

	/**
	 * Getter to the number of vortices data this instance can store
	 *
	 * @returns The number of vortices data this instance is currently built to store
	 */
	size_t getNbVtx() const;

	/**
	 * Getter to the number of steps this instance can store
	 *
	 * @returns The number of time steps this instance is currently built to store
	 */
	size_t getNbSteps() const;
	
	/**
	 * Getter to the time corresponding to a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 */
	double getTimeAt(size_t timeStep) const;
	
	/**
	 * Getter to the fluid identification number of a specific vortex in the simulation stored in the instance
	 *
	 * @param vortexID The number of the vortex
	 */
	size_t getFluidId(size_t vortexID) const;

	/**
	 * Getter to the x-coordinate of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getXAt(size_t timeStep, size_t vortexID) const;
	
	/**
	 * Getter to the y-coordinate of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getYAt(size_t timeStep, size_t vortexID) const;

	/**
	 * Getter to the projection of the velocity along the x-axis of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getUAt(size_t timeStep, size_t vortexID) const;

	/**
	 * Getter to the projection of the velocity along the y-axis of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getVAt(size_t timeStep, size_t vortexID) const;

	/**
	 * Getter to the velocity circulation of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getCirculationAt(size_t timeStep, size_t vortexID) const;

	/**
	 * Getter to the regularization radius of a specific vortex at a specific time step in the simulation stored in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 */
	double getRegRadiusAt(size_t timeStep, size_t vortexID) const;

	/**
	 * Setter for the time corresponding to a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param time The time value to store
	 */
	void storeTimeAt(size_t timeStep, double time);

	/**
	 * Setter for the fluid identification number of a specific vortex to store in the instance
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param fluidId The value to store
	 */
	void storeFluidId(size_t vortexID, size_t fluidId);

	/**
	 * Setter for the x-coordinate of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param x The value to store
	 */
	void storeXAt(size_t timeStep, size_t vortexID, double x);
	
	/**
	 * Setter for the y-coordinate of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param y The value to store
	 */
	void storeYAt(size_t timeStep, size_t vortexID, double y);

	/**
	 * Setter for the projection along the x-axis of the velocity of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param u The value to store
	 */
	void storeUAt(size_t timeStep, size_t vortexID, double u);

	/**
	 * Setter for the projection along the y-axis of the velocity of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param v The value to store
	 */
	void storeVAt(size_t timeStep, size_t vortexID, double v);
	
	/**
	 * Setter for the velocity circulation of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param circulation The value to store
	 */
	void storeCirculationAt(size_t timeStep, size_t vortexID, double circulation);

	/**
	 * Setter for the regularization radius of a specific vortex at a specific time step to store in the instance
	 *
	 * @param timeStep The number of the time step
	 *
	 * @param vortexID The number of the vortex
	 *
	 * @param radius The value to store
	 */
	void storeRegRadiusAt(size_t timeStep, size_t vortexID, double radius);

	/**
	 * Method writing the current state of the instance in a file in order to be able to reconstruct the same instance later with the loadFile method
	 *
	 * @param fileNameWithPath The name of the file with its path
	 */
	void saveState(std::string fileNameWithPath) const;

	/**
	 * Method loading the data previously stored in a file with the method saveState in order to restore a previous state or another instance
	 *
	 * @param fileNameWithPath The name of the file with its path
	 */
	void loadFile(std::string fileNameWithPath);

	/**
	 * Method resizing the different arrays in order to store a simulation with a new number of vortices or a new number of time steps
	 *
	 * @param nb_timeSamples The new number of time samples
	 *
	 * @param nb_vtx The new number of vortices
	 *
	 * @warning nb_timeSamples = nb_steps + 1
	 */
	void reset(size_t nb_timeSamples, size_t nb_vtx);
};