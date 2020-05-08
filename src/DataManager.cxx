#include "DataManager.hxx"

DataManager::DataManager(size_t nb_steps, size_t nb_vtx) :
	m_nb_timeSamples(nb_steps+1),
	m_nb_vtx(nb_vtx),
	ptr_time(new double[m_nb_timeSamples]),
	ptr_x(new double[m_nb_timeSamples*m_nb_vtx]),
	ptr_y(new double[m_nb_timeSamples*m_nb_vtx]),
	ptr_u(new double[m_nb_timeSamples*m_nb_vtx]),
	ptr_v(new double[m_nb_timeSamples*m_nb_vtx]),
	ptr_circulations(new double[m_nb_timeSamples*m_nb_vtx]),
	ptr_radiuses(new double[m_nb_timeSamples*m_nb_vtx]) {}

DataManager::~DataManager() {}

size_t DataManager::getNbVtx() const
{
	return m_nb_vtx;
}

size_t DataManager::getNbSteps() const
{
	return m_nb_timeSamples-1;
}

double DataManager::getTimeAt(size_t timeStep) const
{
	return ptr_time[timeStep];
}

double DataManager::getXAt(size_t timeStep, size_t vortexID) const
{
	return ptr_x[vortexID + m_nb_vtx*timeStep];
}

double DataManager::getYAt(size_t timeStep, size_t vortexID) const
{
	return ptr_y[vortexID + m_nb_vtx*timeStep];
}

double DataManager::getUAt(size_t timeStep, size_t vortexID) const
{
	return ptr_u[vortexID + m_nb_vtx*timeStep];
}

double DataManager::getVAt(size_t timeStep, size_t vortexID) const
{
	return ptr_v[vortexID + m_nb_vtx*timeStep];
}

double DataManager::getCirculationAt(size_t timeStep, size_t vortexID) const
{
	return ptr_circulations[vortexID + m_nb_vtx*timeStep];
}

double DataManager::getRegRadiusAt(size_t timeStep, size_t vortexID) const
{
	return ptr_radiuses[vortexID + m_nb_vtx*timeStep];
}

void DataManager::storeTimeAt(size_t timeStep, double time)
{
	ptr_time[timeStep] = time;
}

void DataManager::storeXAt(size_t timeStep, size_t vortexID, double x)
{
	ptr_x[vortexID + m_nb_vtx*timeStep] = x;
}

void DataManager::storeYAt(size_t timeStep, size_t vortexID, double y)
{
	ptr_y[vortexID + m_nb_vtx*timeStep] = y;
}

void DataManager::storeUAt(size_t timeStep, size_t vortexID, double u)
{
	ptr_u[vortexID + m_nb_vtx*timeStep] = u;
}

void DataManager::storeVAt(size_t timeStep, size_t vortexID, double v)
{
	ptr_v[vortexID + m_nb_vtx*timeStep] = v;
}

void DataManager::storeCirculationAt(size_t timeStep, size_t vortexID, double circulation)
{
	ptr_circulations[vortexID + m_nb_vtx*timeStep] = circulation;
}

void DataManager::storeRegRadiusAt(size_t timeStep, size_t vortexID, double radius)
{
	ptr_radiuses[vortexID + m_nb_vtx*timeStep] = radius;
}

void DataManager::saveState(std::string fileNameWithPath) const
{
	std::ofstream outfile;
	outfile.open(fileNameWithPath);

	if (!outfile.is_open())
	{
		throw std::invalid_argument("Unable to save data in file " + fileNameWithPath);
	}

	outfile << m_nb_timeSamples << "\n";
	outfile << m_nb_vtx << "\n";

	for (size_t t = 0; t < m_nb_timeSamples; t++)
	{
		outfile << ptr_time[t] << "\n";
	}

	size_t const size(m_nb_vtx*m_nb_timeSamples);
	double x, y;

	for (size_t index = 0; index < size; index++)
	{
		x = ptr_x[index];
		y = ptr_y[index];
		outfile << x << "\n";
		outfile << y << "\n";
	}

	for (size_t index = 0; index < size; index++)
	{
		x = ptr_u[index];
		y = ptr_v[index];
		outfile << x << "\n";
		outfile << y << "\n";
	}

	for (size_t index = 0; index < size; index++)
	{
		outfile << ptr_circulations[index] << "\n";
	}

	for (int index = 0; index < size; index++)
	{
		outfile << ptr_radiuses[index] << "\n";
	}

	outfile.close();
}

void DataManager::loadFile(std::string fileNameWithPath)
{
	std::ifstream dataFile;
	dataFile.open(fileNameWithPath);

	if (!dataFile.is_open())
	{
		throw std::invalid_argument("File " + fileNameWithPath + " can't be opened");
	}

	std::string dataStr;

	getline(dataFile, dataStr);
	m_nb_timeSamples = atoi(dataStr.c_str());

	getline(dataFile, dataStr);
	m_nb_vtx = atoi(dataStr.c_str());

	reset(m_nb_timeSamples-1, m_nb_vtx);

	for (size_t t = 0; t < m_nb_timeSamples; t++)
	{
		getline(dataFile, dataStr);
		ptr_time[t] = atof(dataStr.c_str());
	}

	size_t const size(m_nb_vtx*m_nb_timeSamples);
	double x, y;

	for (size_t index = 0; index < size; index++)
	{
		getline(dataFile, dataStr);
		ptr_x[index] = atof(dataStr.c_str());

		getline(dataFile, dataStr);
		ptr_y[index] = atof(dataStr.c_str());
	}

	for (size_t index = 0; index < size; index++)
	{
		getline(dataFile, dataStr);
		ptr_u[index] = atof(dataStr.c_str());

		getline(dataFile, dataStr);
		ptr_v[index] = atof(dataStr.c_str());
	}

	for (size_t index = 0; index < size; index++)
	{
		getline(dataFile, dataStr);
		ptr_circulations[index] = atof(dataStr.c_str());
	}

	for (size_t index = 0; index < size; index++)
	{
		getline(dataFile, dataStr);
		ptr_radiuses[index] = atof(dataStr.c_str());
	}

	dataFile.close();
}

void DataManager::reset(size_t nb_steps, size_t nb_vtx)
{
	m_nb_timeSamples = nb_steps+1;
	m_nb_vtx = nb_vtx;
	ptr_time.reset(new double[m_nb_timeSamples]);
	ptr_x.reset(new double[m_nb_timeSamples*m_nb_vtx]);
	ptr_y.reset(new double[m_nb_timeSamples*m_nb_vtx]);
	ptr_u.reset(new double[m_nb_timeSamples*m_nb_vtx]);
	ptr_v.reset(new double[m_nb_timeSamples*m_nb_vtx]);
	ptr_circulations.reset(new double[m_nb_timeSamples*m_nb_vtx]);
	ptr_radiuses.reset(new double[m_nb_timeSamples*m_nb_vtx]);
}