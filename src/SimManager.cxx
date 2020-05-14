#include "SimManager.hxx"

SimManager::SimManager() :
	name(""),
	kernel(false, false, 0., 0., "eulerExp"),
	data(0, 0),
	afterprocessor() {}

SimManager::~SimManager() {}

void SimManager::setName(std::string new_name)
{
	name = new_name;
}

void SimManager::addVtx(double x, double y, double circulation, double radius, size_t fluidId)
{
	Vortex vtx(x, y, circulation, radius, fluidId);
	kernel.addVortex(vtx);
}

void SimManager::buildTimeSample(double t0, double tEnd, size_t nb_steps)
{
	kernel.buildTimeSample(t0, tEnd, nb_steps);
}

void SimManager::setXPeriodicityTo(bool periodic, double period)
{
	kernel.setXPeriodicityTo(periodic, period);
}

void SimManager::setYPeriodicityTo(bool periodic, double period)
{
	kernel.setYPeriodicityTo(periodic, period);
}

void SimManager::setMethodTo(std::string method)
{
	kernel.setMethodTo(method);
}

void SimManager::sim(size_t nb_threads)
{
	kernel.sim(data, nb_threads);
}

std::string SimManager::getName() const
{
	return name;
}

size_t SimManager::getNbVtx() const
{
	return data.getNbVtx();
}

size_t SimManager::getNbSteps() const
{
	return data.getNbSteps();
}


std::vector<double> SimManager::getTimeVector() const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps()+1; t++)
	{
		res.push_back(data.getTimeAt(t));
	}

	return res;
}

std::vector<double> SimManager::getXsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < data.getNbVtx(); i++)
	{
		res.push_back(data.getXAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getYsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < data.getNbVtx(); i++)
	{
		res.push_back(data.getYAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getUsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < data.getNbVtx(); i++)
	{
		res.push_back(data.getUAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getVsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < data.getNbVtx(); i++)
	{
		res.push_back(data.getVAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getCirculationsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < data.getNbVtx(); i++)
	{
		res.push_back(data.getCirculationAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getVtxXs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps()+1; t++)
	{
		res.push_back(data.getXAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxYs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps()+1; t++)
	{
		res.push_back(data.getYAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxUs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps()+1; t++)
	{
		res.push_back(data.getUAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxVs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps()+1; t++)
	{
		res.push_back(data.getVAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxCirculations(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < data.getNbSteps(); t++)
	{
		res.push_back(data.getCirculationAt(t, vtx));
	}

	return res;
}

std::map<size_t, double> SimManager::computeCompositionAt(double x, double y, size_t step, double radius) const
{
	return DataAnalyst::computeCompositionAt(data, x, y, step, radius);
}

std::vector<double> SimManager::computeVelocityAt(double x, double y, size_t step, bool x_periodic, double x_period) const
{
	double u = DataAnalyst::computeUAt(data, step, x, y, x_periodic, x_period);
	double v = DataAnalyst::computeVAt(data, step, x, y, x_periodic, x_period);
	std::vector<double> res;
	res.push_back(u);
	res.push_back(v);
	res.push_back(sqrt(u*u + v*v));
	res.push_back(Complex(u,v).getArg());
	return res;
}

double SimManager::computeVorticityAt(double x, double y, size_t step, double h, bool x_periodic, double x_period) const
{
	return DataAnalyst::computeVorticityAt(data, step, x, y, x_periodic, x_period, h);
}

std::vector<double> SimManager::computeHamiltonianEvolution(size_t nb_threads) const
{
	return DataAnalyst::computeHamiltonianEvolution(data, nb_threads);
}

void SimManager::saveSim(std::string fileNameWithPath) const
{
	data.saveState(fileNameWithPath);
}

void SimManager::loadSim(std::string fileNameWithPath)
{
	data.loadFile(fileNameWithPath);
}