#pragma once

#include "SimKernel.hxx"
#include "DataManager.hxx"
#include "DataAnalyst.hxx"

/**
 * A class managing the simulation and coordinating the SimKernel, the DataManager and the DataAnalyst objects
 */
class SimManager
{
private:
	std::string m_name;
	bool m_x_periodic;
	double m_x_period;

	SimKernel m_kernel;
	DataManager m_data;
	DataAnalyst m_afterprocessor;

public:
	/**
	 * The constructor of the class
	 */
	SimManager();

	/**
	 * The destructor of the class
	 */
	~SimManager();

	/**
	 * Setter for the name of the simulation
	 *
	 * @param new_name The name to give to the simulation
	 */
	void setName(std::string new_name);

	/**
	 * Method to add a vortex to the set of vortices that will be used for the simulation, which is by default empty when the SimManager object is created
	 *
	 * @param x The x-coordinate of the vortex to add
	 *
	 * @param y The y-coordinate of the vortex to add
	 *
	 * @param circulation The velocity circulation around the vortex to add
	 *
	 * @param radius The regularization radius around the vortex where the induced velocities will be restrained to avoid singularities (must be greater than 0)
	 *
	 * @param fluidId The number of identification of the fluid the vortex to add belongs to
	 */
	void addVtx(double x, double y, double circulation, double radius, size_t fluidId);

	/**
	 * Method defining all the time steps which will be used for the simulation
	 *
	 * @param t0 The time when the simulation will start
	 *
	 * @param tEnd The time when the simulation will end (must be greater than t0)
	 *
	 * @param nb_steps The number of time steps that will be computed to go from t0 to tEnd (must be greater than 1)
	 */
	void buildTimeSample(double t0, double tEnd, size_t nb_steps);

	/**
	 * Method defining the eventual periodicity along the x-axis of the flow to simulate
	 *
	 * @param periodic A boolean telling if the flow is periodic
	 *
	 * @param period The value of the eventual spatial period (must be greater than 0)
	 */
	void setXPeriodicityTo(bool periodic, double period);
	
	/**
	 * Method defining the eventual periodicity along the y-axis of the flow to simulate
	 *
	 * @param periodic A boolean telling if the flow is periodic
	 *
	 * @param period The value of the eventual spatial period (must be greater than 0)
	 */
	void setYPeriodicityTo(bool periodic, double period);

	/**
	 * Method defining the temporal integration method to use for the simulation
	 *
	 * @param method The name of the method to use (Possible values : 'euler', 'rk4', 'eulerA', 'eulerB', 'sv', 'svi')
	 */
	void setMethodTo(std::string method);

	/**
	 * Method performing the simulation using all the parameters previously defined
	 *
	 * @param nb_threads The number of threads to use to speed-up the computation (must be greater than 1)
	 *
	 * @warning The method addVtx must be called for each vortex to use in the simulation before calling this method
	 *
	 * @warning The method buildTimeSample must be called once before calling this method
	 *
	 * @warning If the methods setXPeriodicity and setYPeriodicity are not used before calling this method, the flow will by default be simulated as non-periodic
	 *
	 * @warning This method can be long to execute, as the computation time is proportional to the number of steps and the square of the number of vortices
	 */
	void sim(size_t nb_threads);

	/**
	 * Getter to the name given to the simulation
	 *
	 * @returns The name of the simulation as a std::string
	 */
	std::string getName() const;

	/**
	 * Getter to the number of vortices currently added to the simulation
	 *
	 * @returns The number of vortices of the simulation
	 */
	size_t getNbVtx() const;

	/**
	 * Getter to the number of steps of the simulation
	 *
	 * @returns The number of steps of the simulation
	 */
	size_t getNbSteps() const;

	/**
	 * Getter to the time values of each time step of the simulation
	 *
	 * @returns A std::vector containing all the time values of the simulation
	 */
	std::vector<double> getTimeVector() const;

	/**
	 * Getter to the x-coordinates of each vortex at a specific time-step
	 *
	 * @param step The number of the time step
	 *
	 * @returns A std::vector containing all the x-coordinates of the vortices in the same order as the one of addition to the simulation
	 */
	std::vector<double> getXsAt(size_t step) const;

	/**
	 * Getter to the y-coordinates of each vortex at a specific time-step
	 *
	 * @param step The number of the time step
	 *
	 * @returns A std::vector containing all the y-coordinates of the vortices in the same order as the one of addition to the simulation
	 */
	std::vector<double> getYsAt(size_t step) const;

	/**
	 * Getter to the velocity of each vortex projected on the x-axis at a specific time-step
	 *
	 * @param step The number of the time step
	 *
	 * @returns A std::vector containing all the velocity projections of the vortices in the same order as the one of addition to the simulation
	 */
	std::vector<double> getUsAt(size_t step) const;

	/**
	 * Getter to the velocity of each vortex projected on the y-axis at a specific time-step
	 *
	 * @param step The number of the time step
	 *
	 * @returns A std::vector containing all the velocity projections of the vortices in the same order as the one of addition to the simulation
	 */
	std::vector<double> getVsAt(size_t step) const;

	/**
	 * Getter to the velocity circulation of each vortex at a specific time-step
	 *
	 * @param step The number of the time step
	 *
	 * @returns A std::vector containing all the velocity circulations of the vortices in the same order as the one of addition to the simulation
	 */
	std::vector<double> getCirculationsAt(size_t step) const;

	/**
	 * Getter to all the succesive x-coordinates of a specific vortex during the simulation
	 *
	 * @param vtx The number of the vortex (ie the number of addition to the simulation of the vortex)
	 *
	 * @returns A std::vector containing the x-coordinate of the vortex at each time step
	 */
	std::vector<double> getVtxXs(size_t vtx) const;

	/**
	 * Getter to all the succesive y-coordinates of a specific vortex during the simulation
	 *
	 * @param vtx The number of the vortex (ie the number of addition to the simulation of the vortex)
	 *
	 * @returns A std::vector containing the y-coordinate of the vortex at each time step
	 */
	std::vector<double> getVtxYs(size_t vtx) const;

	/**
	 * Getter to all the succesive velocity projections on the x-axis of a specific vortex during the simulation
	 *
	 * @param vtx The number of the vortex (ie the number of addition to the simulation of the vortex)
	 *
	 * @returns A std::vector containing the velocity projection on the x-axis of the vortex at each time step
	 */
	std::vector<double> getVtxUs(size_t vtx) const;

	/**
	 * Getter to all the succesive velocity projections on the y-axis of a specific vortex during the simulation
	 *
	 * @param vtx The number of the vortex (ie the number of addition to the simulation of the vortex)
	 *
	 * @returns A std::vector containing the velocity projection on the y-axis of the vortex at each time step
	 */
	std::vector<double> getVtxVs(size_t vtx) const;

	/**
	 * Getter to all the succesive values of velocity circulation of a specific vortex during the simulation
	 *
	 * @param vtx The number of the vortex (ie the number of addition to the simulation of the vortex)
	 *
	 * @returns A std::vector containing the velocity circulations of the vortex at each time step
	 */
	std::vector<double> getVtxCirculations(size_t vtx) const;

	/**
	 * Method computing the composition at a specific point in space and a specific time step
	 * 
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param step The time step at which the computation will be done
	 *
	 * @param radius The radius around the point (x,y) to consider to perform the computation
	 */
	std::map<size_t, double> computeCompositionAt(double x, double y, size_t step, double radius) const;

	/**
	 * Method computing the velocity at a specific point in space and a specific time step
	 * 
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param step The time step at which the computation will be done
	 *
	 * @param x_periodic A boolean telling if the velocity should be computed considering an x-periodic flow
	 *
	 * @param x_period The eventual period along the x-axis of the flow (must be positive)
	 *
	 * @returns A std::vector containing 4 elements : the projection on the x-axis, the projection on the y-axis, the module and the argument of the computed velocity
	 *
	 * @throws std::invalid_argument Thrown if one of the components of the velocity computed is NaN
	 */
	std::vector<double> computeVelocityAt(double x, double y, size_t step, bool x_periodic, double x_period) const;

	/**
	 * Method computing the vorticity at a specific point in space and a specific time step
	 * 
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param step The time step at which the computation will be done
	 *
	 * @param h The spatial step to use for the numerical derivation of the velocity field
	 *
	 * @param x_periodic A boolean telling if the velocity should be computed considering an x-periodic flow
	 *
	 * @param x_period The eventual period along the x-axis of the flow (must be positive)
	 *
	 * @returns The vorticity computed (only one component, the one on the z-axis, is returned because the vorticity vector if always along the z-axis)
	 */
	double computeVorticityAt(double x, double y, size_t step, double h, bool x_periodic, double x_period) const;

	/**
	 * Method computing the evolution of the Hamiltonian over the simulation
	 *
	 * @param nb_threads The number of threads to use to perform the computation (must be greater than 1)
	 *
	 * @returns A std::vector containing the value of the Hamiltonian computed for each time step
	 *
	 * @warning This method can be long to exectute for big simulations since the computation time is proportional to the number of time steps and the square of the number of vortices
	 */
	std::vector<double> computeHamiltonianEvolution(size_t nb_threads) const;

	/**
	 * Method saving the current state of the instance to a file in order to restore it later with the method loadSim
	 *
	 * @param fileNameWithPath The name of the file to save in with the full path leading to it
	 */
	void saveSim(std::string fileNameWithPath) const;

	/**
	 * Method restoring a previous state or another instance saved in a file with the method saveSim
	 *
	 * @param fileNameWithPath The nameof the file to read with the full path leading to it
	 */
	void loadSim(std::string fileNameWithPath);
};