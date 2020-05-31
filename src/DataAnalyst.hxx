#pragma once

#include <vector>
#include <thread>
#include <map>
#include "DataManager.hxx"
#include "Vortex.hxx"

/**
 * A class providing different methods to post-process simulation data
 *
 * @todo Put a DataManager reference as member => no more static methods
 */
class DataAnalyst
{
private:
	static void computeHamiltonianEvolutionBetween(DataManager const &dm, size_t start_step, size_t end_step, std::vector<double> &v_H, bool x_periodic, double x_period);

public:
	/**
	 * The constructor of the class
	 */
	DataAnalyst();

	/**
	 * The destructor of the class
	 */
	~DataAnalyst();

	/**
	 * Method computing the Hamiltonian of the point vortex system at a given time step
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param time_index The targetted time step to perform the computation
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The enventual spatial period along the x-axis (must be positive)
	 *
	 * @returns The value of the Hamiltonian
	 */
	static double computeHamiltonianAt(DataManager const &dm, size_t time_index, bool x_periodic, double x_period);
	
	/**
	 * Method computing the evolution of the Hamiltonian of the point vortex system over the simulation
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param nb_threads The number of threads to use for computation (must be greater than 1)
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The enventual spatial period along the x-axis (must be positive)
	 *
	 * @returns A std::vector containing the value of the Hamiltonian at each time step
	 */
	static std::vector<double> computeHamiltonianEvolution(DataManager const &dm, size_t nb_threads, bool x_periodic, double x_period);

	/**
	 * Method computing the composition of the flow at a given point in space ans at a given time step
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param step The time step to use for the computation
	 * 
	 * @param radius The radius around the point (x,y) to consider in order to perform the computation
	 *
	 * @returns A std::map<size_t,double> where the keys are the fluid identification numbers and the values are the proportion of this fluid at the point of computation
	 */
	static std::map<size_t, double> computeCompositionAt(DataManager const &dm, double x, double y, size_t step, double radius);

	/**
	 * Method computing the evolution of the composition at a given point in space over the simulation
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param x The x-coordinate of the point where the computations will be done
	 *
	 * @param y The y-coordinate of the point where the computations will be done
	 *
	 * @param radius The radius aroud the point (x,y) to consider in order to perform the computations
	 *
	 * @returns A std::vector containing all the std::map<size_t,double> computed for each time step by the computeCompositionAt method
	 */
	static std::vector<std::map<size_t, double>> computeCompositionEvolutionAt(DataManager const &dm, double x, double y, double radius);
	
	/**
	 * Method computing the projection on the x-axis of the velocity induced by the whole system of point vortices at a specific point in space and at a specific time step
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param time_index The time step where the computation will be done
	 *
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @returns The computed velocity projection on the x-axis
	 */
	static double computeUAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period);
	
	/**
	 * Method computing the evolution of the velocity projected on the x-axis at a given point in space over the simulation
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param x The x-coordinate of the point where the computations will be done
	 *
	 * @param y The y-coordinate of the point where the computations will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @returns A std::vector of the computed velocity projections for each time step
	 */
	static std::vector<double> computeUEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period);
	
	/**
	 * Method computing the projection on the y-axis of the velocity induced by the whole system of point vortices at a specific point in space and at a specific time step
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param time_index The time step where the computation will be done
	 *
	 * @param x The x-coordinate of the point where the computation will be done
	 *
	 * @param y The y-coordinate of the point where the computation will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @returns The computed velocity projection on the y-axis
	 */
	static double computeVAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period);
	
	/**
	 * Method computing the evolution of the velocity projected on the y-axis at a given point in space over the simulation
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param x The x-coordinate of the point where the computations will be done
	 *
	 * @param y The y-coordinate of the point where the computations will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @returns A std::vector of the computed velocity projections for each time step
	 */
	static std::vector<double> computeVEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period);

	/**
	 * Method computing the vorticity at a given point in space and a given time step
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param time_index The time step where the computation will be done
	 *
	 * @param x The x-coordinate of the point where the computations will be done
	 *
	 * @param y The y-coordinate of the point where the computations will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @param h The spatial step to use for the numerical derivation of the velocity field
	 *
	 * @returns The computed vorticity at the given point in time and space
	 */
	static double computeVorticityAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period, double h);
	
	/**
	 * Method computing the evolution of the vorticity at a given point in space over the simulation
	 *
	 * @param dm The DataManager object where the simulation is stored
	 *
	 * @param x The x-coordinate of the point where the computations will be done
	 *
	 * @param y The y-coordinate of the point where the computations will be done
	 *
	 * @param x_periodic A boolean telling if the flow is periodic along the x-axis or not
	 *
	 * @param x_period The eventual spatial period along the x-axis (must be positive)
	 *
	 * @param h The spatial step to use for the numerical derivation of the velocity field
	 *
	 * @returns A std::vector of the computed vorticities for each time step
	 */
	static std::vector<double> computeVorticityEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period, double h);
};