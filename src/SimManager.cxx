#include "SimManager.hxx"

SimManager::SimManager() :
	m_name(""),
	m_x_periodic(false),
	m_x_period(0.),
	m_kernel(false, false, 0., 0., "euler"),
	m_data(0, 0),
	m_afterprocessor() {}

SimManager::~SimManager() {}

void SimManager::setName(std::string new_name)
{
	m_name = new_name;
}

void SimManager::addVtx(double x, double y, double circulation, double radius, size_t fluidId)
{
	Vortex vtx(x, y, circulation, radius, fluidId);
	m_kernel.addVortex(vtx);
}

void SimManager::buildTimeSample(double t0, double tEnd, size_t nb_steps)
{
	m_kernel.buildTimeSample(t0, tEnd, nb_steps);
}

void SimManager::setXPeriodicityTo(bool periodic, double period)
{
	m_kernel.setXPeriodicityTo(periodic, period);
}

void SimManager::setYPeriodicityTo(bool periodic, double period)
{
	m_x_periodic = periodic;
	m_x_period = period;
	m_kernel.setYPeriodicityTo(periodic, period);
}

void SimManager::setMethodTo(std::string method)
{
	m_kernel.setMethodTo(method);
}

void SimManager::sim(size_t nb_threads)
{
	m_kernel.sim(m_data, nb_threads);
}

std::string SimManager::getName() const
{
	return m_name;
}

size_t SimManager::getNbVtx() const
{
	return m_data.getNbVtx();
}

size_t SimManager::getNbSteps() const
{
	return m_data.getNbSteps();
}


std::vector<double> SimManager::getTimeVector() const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps()+1; t++)
	{
		res.push_back(m_data.getTimeAt(t));
	}

	return res;
}

std::vector<double> SimManager::getXsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < m_data.getNbVtx(); i++)
	{
		res.push_back(m_data.getXAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getYsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < m_data.getNbVtx(); i++)
	{
		res.push_back(m_data.getYAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getUsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < m_data.getNbVtx(); i++)
	{
		res.push_back(m_data.getUAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getVsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < m_data.getNbVtx(); i++)
	{
		res.push_back(m_data.getVAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getCirculationsAt(size_t step) const
{
	std::vector<double> res;

	for (size_t i = 0; i < m_data.getNbVtx(); i++)
	{
		res.push_back(m_data.getCirculationAt(step, i));
	}

	return res;
}

std::vector<double> SimManager::getVtxXs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps()+1; t++)
	{
		res.push_back(m_data.getXAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxYs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps()+1; t++)
	{
		res.push_back(m_data.getYAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxUs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps()+1; t++)
	{
		res.push_back(m_data.getUAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxVs(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps()+1; t++)
	{
		res.push_back(m_data.getVAt(t, vtx));
	}

	return res;
}

std::vector<double> SimManager::getVtxCirculations(size_t vtx) const
{
	std::vector<double> res;

	for (size_t t = 0; t < m_data.getNbSteps(); t++)
	{
		res.push_back(m_data.getCirculationAt(t, vtx));
	}

	return res;
}

std::map<size_t, double> SimManager::computeCompositionAt(double x, double y, size_t step, double radius) const
{
	return DataAnalyst::computeCompositionAt(m_data, x, y, step, radius);
}

std::vector<double> SimManager::computeVelocityAt(double x, double y, size_t step, bool x_periodic, double x_period) const
{
	double u = DataAnalyst::computeUAt(m_data, step, x, y, x_periodic, x_period);
	double v = DataAnalyst::computeVAt(m_data, step, x, y, x_periodic, x_period);
	std::vector<double> res;
	res.push_back(u);
	res.push_back(v);
	res.push_back(sqrt(u*u + v*v));
	res.push_back(Complex(u,v).getArg());
	return res;
}

double SimManager::computeVorticityAt(double x, double y, size_t step, double h, bool x_periodic, double x_period) const
{
	return DataAnalyst::computeVorticityAt(m_data, step, x, y, x_periodic, x_period, h);
}

std::vector<double> SimManager::computeHamiltonianEvolution(size_t nb_threads) const
{
	return DataAnalyst::computeHamiltonianEvolution(m_data, nb_threads, m_x_periodic, m_x_period);
}

void SimManager::saveSim(std::string fileNameWithPath) const
{
	m_data.saveState(fileNameWithPath);
}

void SimManager::loadSim(std::string fileNameWithPath)
{
	m_data.loadFile(fileNameWithPath);
}