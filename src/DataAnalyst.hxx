#pragma once

#include <vector>
#include <thread>
#include <map>
#include "DataManager.hxx"
#include "Vortex.hxx"

class DataAnalyst
{
public:
	DataAnalyst();
	~DataAnalyst();

	static double computeHamiltonianAt(DataManager const &dm, size_t time_index);
	static void computeHamiltonianEvolutionBetween(DataManager const &dm, size_t start_step, size_t end_step, std::vector<double> &v_H);
	static std::vector<double> computeHamiltonianEvolution(DataManager const &dm, size_t nb_threads);

	static std::map<size_t, double> computeCompositionAt(DataManager const &dm, double x, double y, size_t step, double radius);
	static std::vector<std::map<size_t, double>> computeCompositionEvolutionAt(DataManager const &dm, double x, double y, double radius);
	
	static double computeUAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period);
	static std::vector<double> computeUEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period);
	
	static double computeVAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period);
	static std::vector<double> computeVEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period);

	static double computeVorticityAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period, double h);
	static std::vector<double> computeVorticityEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period, double h);
};