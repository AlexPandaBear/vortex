#include "SimKernel.hxx"

SimKernel::SimKernel(bool x_periodic, bool y_periodic, double x_period, double y_period, std::string method) :
	m_X_periodic(x_periodic),
	m_Y_periodic(y_periodic),
	m_X_period(x_period),
	m_Y_period(y_period),
	m_method(method) {}

SimKernel::~SimKernel() {}

void SimKernel::addVortex(Vortex &v)
{
	v_vtx.push_back(v);
}

void SimKernel::buildTimeSample(double t0, double tEnd, size_t steps)
{
	double deltaT(tEnd - t0);
	if (!(deltaT > 0.))
	{
		throw std::invalid_argument("t0 must be smaller than tEnd");
	}

	for (size_t i = 0; i <= steps; i++)
	{
		v_time.push_back(t0 + ((double) i/steps) * deltaT);
	}
}

void SimKernel::setXPeriodicityTo(bool periodic, double period)
{
	m_X_periodic = periodic;
	m_X_period = period;
}

void SimKernel::setYPeriodicityTo(bool periodic, double period)
{
	m_Y_periodic = periodic;
	m_Y_period = period;
}

void SimKernel::setMethodTo(std::string method)
{
	m_method = method;
}

bool SimKernel::readyToSim()
{
	if (m_X_period)
	{
		double x;

		for (auto &v : v_vtx)
		{
			x = v.getX();
			if (x > m_X_period || x < 0.)
			{
				return false;
			}
		}
	}

	if (m_Y_period)
	{
		double y;

		for (auto &v : v_vtx)
		{
			y = v.getY();
			if (y > m_Y_period || y < 0.)
			{
				return false;
			}
		}
	}

	return true;
}

void SimKernel::integrate(DataManager &dm, size_t step)
{
	double dt(v_time[step+1] - v_time[step]);

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		v_vtx[i].move(dm.getUAt(step, i) * dt, dm.getVAt(step, i) * dt);
		dm.storeXAt(step+1, i, v_vtx[i].getX());
		dm.storeYAt(step+1, i, v_vtx[i].getY());
		dm.storeCirculationAt(step+1, i, v_vtx[i].getCirculation());
		dm.storeRegRadiusAt(step+1, i, v_vtx[i].getRegRadius());
	}
}

void SimKernel::computeEEStep(DataManager &dm, size_t step, size_t firstVtx, size_t lastVtx) const
{
	double x, y, u, v;
	Vortex vtx(0,0,0,0);

	for (size_t i = firstVtx; i < lastVtx; i++)
	{
		x = v_vtx[i].getX();
		y = v_vtx[i].getY();
		u = 0.;
		v = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = v_vtx[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u);
		dm.storeVAt(step, i, v);
	}
}

void SimKernel::computeEEStep_multithread(DataManager &dm, size_t step, size_t nb_threads)
{
	if (nb_threads < 1)
	{
		std::ostringstream ss;
		ss << "Unable to perform computation on " << nb_threads << " threads";
		throw std::invalid_argument(ss.str());
	}

	std::vector<std::thread> v_threads;
	size_t firstVtx(0), lastVtx(m_nb_vtx/nb_threads);

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads.push_back(std::thread(&SimKernel::computeEEStep, this, std::ref(dm), step, firstVtx, lastVtx));
		firstVtx = lastVtx;
		lastVtx += m_nb_vtx/nb_threads;
	}

	lastVtx = m_nb_vtx;
	computeEEStep(dm, step, firstVtx, lastVtx);

	for (auto &th : v_threads) {th.join();}
	integrate(dm, step);
}

void SimKernel::computeRK4Substep(std::vector<Vortex> &workingCopy, std::vector<double> &k_u, std::vector<double> &k_v, size_t step, size_t firstVtx, size_t lastVtx, size_t substep) const
{
	double x, y;
	Vortex vtx(0,0,0,0);

	for (size_t i = firstVtx; i < lastVtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		k_u[i] = 0.;
		k_v[i] = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				k_u[i] += vtx.computeXInducedVelocityAt(x, y);
				k_v[i] += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				k_u[i] += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				k_u[i] += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				k_u[i] += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				k_u[i] += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
				k_v[i] += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				k_v[i] += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				k_v[i] += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				k_v[i] += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}
	}
}

void SimKernel::computeRK4Step_multithread(DataManager &dm, size_t step, size_t nb_threads)
{
	if (nb_threads < 1)
	{
		std::ostringstream ss;
		ss << "Unable to perform computation on " << nb_threads << " threads";
		throw std::invalid_argument(ss.str());
	}

	//association of vortices to threads

	std::vector<size_t> v_separations;

	for (size_t i = 0; i < nb_threads; i++)
	{
		v_separations.push_back(i*m_nb_vtx/nb_threads);
	}

	//initialisation

	std::vector<double> k1_u(m_nb_vtx, 0.), k2_u(m_nb_vtx, 0.), k3_u(m_nb_vtx, 0.), k4_u(m_nb_vtx, 0.);
	std::vector<double> k1_v(m_nb_vtx, 0.), k2_v(m_nb_vtx, 0.), k3_v(m_nb_vtx, 0.), k4_v(m_nb_vtx, 0.);

	double dt(v_time[step+1] - v_time[step]);

	std::vector<Vortex> workingCopy;
	for (auto v : v_vtx) {workingCopy.push_back(v);}

	//computation of k1

	std::vector<std::thread> v_threads_k1;

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads_k1.push_back(std::thread(&SimKernel::computeRK4Substep, this, std::ref(workingCopy), std::ref(k1_u), std::ref(k1_v), step, v_separations[i], v_separations[i+1], 1));
	}
	computeRK4Substep(workingCopy, k1_u, k1_v, step, v_separations[nb_threads-1], m_nb_vtx, 1);

	for (auto &th : v_threads_k1) {th.join();}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		workingCopy[i].move(k1_u[i]*0.5*dt, k1_v[i]*0.5*dt);
	}

	//computation of k2

	std::vector<std::thread> v_threads_k2;

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads_k2.push_back(std::thread(&SimKernel::computeRK4Substep, this, std::ref(workingCopy), std::ref(k2_u), std::ref(k2_v), step, v_separations[i], v_separations[i+1], 2));
	}
	computeRK4Substep(workingCopy, k2_u, k2_v, step, v_separations[nb_threads-1], m_nb_vtx, 2);

	for (auto &th : v_threads_k2) {th.join();}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		workingCopy[i].setX(v_vtx[i].getX() + k2_u[i]*0.5*dt);
		workingCopy[i].setY(v_vtx[i].getY() + k2_v[i]*0.5*dt);
	}

	//computation of k3

	std::vector<std::thread> v_threads_k3;

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads_k3.push_back(std::thread(&SimKernel::computeRK4Substep, this, std::ref(workingCopy), std::ref(k3_u), std::ref(k3_v), step, v_separations[i], v_separations[i+1], 3));
	}
	computeRK4Substep(workingCopy, k3_u, k3_v, step, v_separations[nb_threads-1], m_nb_vtx, 3);

	for (auto &th : v_threads_k3) {th.join();}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		workingCopy[i].setX(v_vtx[i].getX() + k3_u[i]*dt);
		workingCopy[i].setY(v_vtx[i].getY() + k3_v[i]*dt);
	}

	//computation of k4

	std::vector<std::thread> v_threads_k4;

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads_k4.push_back(std::thread(&SimKernel::computeRK4Substep, this, std::ref(workingCopy), std::ref(k4_u), std::ref(k4_v), step, v_separations[i], v_separations[i+1], 4));
	}
	computeRK4Substep(workingCopy, k4_u, k4_v, step, v_separations[nb_threads-1], m_nb_vtx, 4);

	for (auto &th : v_threads_k4) {th.join();}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		dm.storeUAt(step, i, (k1_u[i] + k2_u[i]*2. + k3_u[i]*2. + k4_u[i])/6.);
		dm.storeVAt(step, i, (k1_v[i] + k2_v[i]*2. + k3_v[i]*2. + k4_v[i])/6.);
	}

	integrate(dm, step);
}

void SimKernel::computeEAStep(DataManager &dm, size_t step)
{
	double x, y, u, v;
	Vortex vtx(0,0,0,0);

	double dt(v_time[step+1] - v_time[step]);

	std::vector<Vortex> workingCopy;
	for (auto v : v_vtx) {workingCopy.push_back(v);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		u = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(dt*dm.getUAt(step, i), 0.);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		v = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeVAt(step, i, v);
	}

	integrate(dm, step);
}

void SimKernel::computeEBStep(DataManager &dm, size_t step)
{
	double x, y, u, v;
	Vortex vtx(0,0,0,0);

	double dt(v_time[step+1] - v_time[step]);

	std::vector<Vortex> workingCopy;
	for (auto v : v_vtx) {workingCopy.push_back(v);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		v = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeVAt(step, i, v);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(0., dt*dm.getVAt(step, i));}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		u = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u);
	}

	integrate(dm, step);
}

void SimKernel::computeSVStep(DataManager &dm, size_t step)
{
	double x, y, u, v;
	Vortex vtx(0,0,0,0);

	double dt(v_time[step+1] - v_time[step]);

	std::vector<Vortex> workingCopy;
	for (auto v : v_vtx) {workingCopy.push_back(v);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		u = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(0.5*dt*dm.getUAt(step, i), 0.);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		v = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeVAt(step, i, v);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(0., dt*dm.getVAt(step, i));}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		u = dm.getUAt(step, i);

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u/2.);
	}

	integrate(dm, step);
}

void SimKernel::computeSVIStep(DataManager &dm, size_t step)
{
	double x, y, u, v;
	Vortex vtx(0,0,0,0);

	double dt(v_time[step+1] - v_time[step]);

	std::vector<Vortex> workingCopy;
	for (auto v : v_vtx) {workingCopy.push_back(v);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		v = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeVAt(step, i, v);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(0., 0.5*dt*dm.getVAt(step, i));}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		u = 0.;

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				u += vtx.computeXInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				u += vtx.computeXInducedVelocityAt(x + m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x + 2*m_X_period, y);
				u += vtx.computeXInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeUAt(step, i, u);
	}

	for (size_t i = 0; i < m_nb_vtx; i++) {workingCopy[i].move(dt*dm.getUAt(step, i), 0.);}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		x = workingCopy[i].getX();
		y = workingCopy[i].getY();
		v = dm.getVAt(step, i);

		for (size_t j = 0; j < m_nb_vtx; j++)
		{
			vtx = workingCopy[j];
			if (j!=i)
			{
				v += vtx.computeYInducedVelocityAt(x, y);
			}

			if ((i!=j) && m_X_periodic)
			{
				v += vtx.computeYInducedVelocityAt(x + m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x + 2*m_X_period, y);
				v += vtx.computeYInducedVelocityAt(x - 2*m_X_period, y);
			}
		}

		dm.storeVAt(step, i, v/2.);
	}

	integrate(dm, step);
}

void SimKernel::printSimProgression(size_t step) const
{
	std::cout << "\rComputing step " << step+1 << " out of " << m_nb_steps << " -- " << 100.*(step+1)/m_nb_steps << "% completed     " << std::flush;
}

void SimKernel::sim(DataManager &dm, size_t nb_threads)
{
	m_nb_vtx = v_vtx.size();
	m_nb_steps = v_time.size()-1;

	if (!readyToSim())
	{
		throw std::invalid_argument("Bad choice of vortices");
	}

	dm.reset(m_nb_steps, m_nb_vtx);

	for (size_t t = 0; t < m_nb_steps+1; t++)
	{
		dm.storeTimeAt(t, v_time[t]);
	}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		dm.storeFluidId(i, v_vtx[i].getFluidId());
	}

	for (size_t i = 0; i < m_nb_vtx; i++)
	{
		dm.storeXAt(0, i, v_vtx[i].getX());
		dm.storeYAt(0, i, v_vtx[i].getY());
		dm.storeCirculationAt(0, i, v_vtx[i].getCirculation());
		dm.storeRegRadiusAt(0, i, v_vtx[i].getRegRadius());
	}


	if (m_method == "euler")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeEEStep_multithread(dm, t, nb_threads);
		}
	}

	else if (m_method == "rk4")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeRK4Step_multithread(dm, t, nb_threads);
		}
	}

	else if (m_method == "eulerA")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeEAStep(dm, t);
		}
	}

	else if (m_method == "eulerB")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeEBStep(dm, t);
		}
	}

	else if (m_method == "sv")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeSVStep(dm, t);
		}
	}

	else if (m_method == "svi")
	{
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimProgression(t);
			computeSVIStep(dm, t);
		}
	}

	else
	{
		throw std::invalid_argument("Unknown numerical method " + m_method);
	}



/*
	switch (m_method)
	{
		case "euler":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeEEStep_multithread(dm, t, nb_threads);
			}

		case "rk4":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeRK4Step_multithread(dm, t, nb_threads);
			}

		case "eulerA":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeEAStep(dm, t);
			}

		case "eulerB":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeEBStep(dm, t);
			}

		case "sv":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeSVStep(dm, t);
			}

		case "svi":
			for (size_t t = 0; t < m_nb_steps; t++)
			{
				printSimProgression(t);
				computeSVIStep(dm, t);
			}

		default:
			throw std::invalid_argument("Unknown numerical method " + m_method);
	}
*/
	
	std::cout << std::endl;
}