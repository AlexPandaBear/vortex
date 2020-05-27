#include "DataAnalyst.hxx"

DataAnalyst::DataAnalyst() {}

DataAnalyst::~DataAnalyst() {}

double DataAnalyst::computeHamiltonianAt(DataManager const &dm, size_t time_index, bool x_periodic, double x_period)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	double H(0.);
	double deltaX, deltaY;

	if (x_periodic)
	{
		if (!(x_period > 0.))
		{
			std::ostringstream ss;
			ss << "Unable to perform computation with this value of x-period : " << x_period;
			throw std::invalid_argument(ss.str());
		}

		double factor(2*M_PI/x_period);

		for (size_t i = 0; i < nb_vtx; i++)
		{
			for (size_t j = i+1; j < nb_vtx; j++)
			{
				deltaX = dm.getXAt(time_index, i) - dm.getXAt(time_index, j);
				deltaY = dm.getYAt(time_index, i) - dm.getYAt(time_index, j);
				H += dm.getCirculationAt(time_index, i)*dm.getCirculationAt(time_index, j) * log(cosh(factor*deltaY) - cos(factor*deltaX));
			}
		}
	}

	else
	{
		for (size_t i = 0; i < nb_vtx; i++)
		{
			for (size_t j = i+1; j < nb_vtx; j++)
			{
				deltaX = dm.getXAt(time_index, i) - dm.getXAt(time_index, j);
				deltaY = dm.getYAt(time_index, i) - dm.getYAt(time_index, j);
				H += dm.getCirculationAt(time_index, i)*dm.getCirculationAt(time_index, j)*0.5*log(deltaX*deltaX + deltaY*deltaY);
			}
		}
	}

	return H/(4*M_PI);
}

void DataAnalyst::computeHamiltonianEvolutionBetween(DataManager const &dm, size_t start_step, size_t end_step, std::vector<double> &v_H, bool x_periodic, double x_period)
{
	for (size_t t = start_step; t < end_step; t++)
	{
		v_H[t] = computeHamiltonianAt(dm, t, x_periodic, x_period);
	}
}

std::vector<double> DataAnalyst::computeHamiltonianEvolution(DataManager const &dm, size_t nb_threads, bool x_periodic, double x_period)
{
	if (nb_threads < 1)
	{
		std::ostringstream ss;
		ss << "Unable to perform computation on " << nb_threads << " threads";
		throw std::invalid_argument(ss.str());
	}

	size_t nb_steps = dm.getNbSteps();
	std::vector<double> v_H(nb_steps+1, 0.);
	
	std::vector<std::thread> v_threads;
	size_t start_step(0), end_step(nb_steps/nb_threads);

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads.push_back(std::thread(&DataAnalyst::computeHamiltonianEvolutionBetween, std::ref(dm), start_step, end_step, std::ref(v_H), x_periodic, x_period));
		start_step = end_step;
		end_step += nb_steps/nb_threads;
	}

	end_step = nb_steps+1;
	computeHamiltonianEvolutionBetween(dm, start_step, end_step, v_H, x_periodic, x_period);

	for (auto &th : v_threads) {th.join();}

	return v_H;
}

std::map<size_t, double> DataAnalyst::computeCompositionAt(DataManager const &dm, double x, double y, size_t step, double radius)
{
	size_t nb_vtx = dm.getNbVtx();
	size_t fluid_id;
	std::map<size_t, double> map_compo;

	for (size_t i = 0; i < nb_vtx; i++)
	{
		fluid_id = dm.getFluidId(i);
		if (map_compo.find(fluid_id) == map_compo.end())
		{
			map_compo[fluid_id] = 0.;
		}
	}

	double deltaX, deltaY;
	double dist, dist2, rad2(radius*radius);

	for (size_t i = 0; i < nb_vtx; i++)
	{
		deltaX = dm.getXAt(step, i) - x;
		deltaY = dm.getYAt(step, i) - y;
		dist2 = deltaX*deltaX + deltaY*deltaY;
		
		if (dist2 < rad2)
		{
			dist = sqrt(dist2);
			fluid_id = dm.getFluidId(i);
			map_compo[fluid_id] += 1. - dist/radius;
		}
	}

	double sum(0.);
	for (std::map<size_t, double>::iterator it = map_compo.begin(); it != map_compo.end(); it++) {sum += it->second;}
	for (std::map<size_t, double>::iterator it = map_compo.begin(); it != map_compo.end(); it++) {it->second /= sum;}

	return map_compo;
}

std::vector<std::map<size_t, double>> DataAnalyst::computeCompositionEvolutionAt(DataManager const &dm, double x, double y, double radius)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	std::vector<std::map<size_t, double>> v_C;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_C.push_back(computeCompositionAt(dm, x, y, t, radius));
	}
	return v_C;
}

double DataAnalyst::computeUAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	double u(0.);
	double vtx_x, vtx_y, vtx_circ, vtx_rad;

	Vortex vtx(0., 0., 0., 0.);
	
	for (size_t i = 0; i < nb_vtx; i++)
	{
		vtx_x = dm.getXAt(time_index, i);
		vtx_y = dm.getYAt(time_index, i);
		vtx_circ = dm.getCirculationAt(time_index, i);
		vtx_rad = dm.getRegRadiusAt(time_index, i);

		vtx = Vortex(vtx_x, vtx_y, vtx_circ, vtx_rad, 0, x_periodic, x_period);

		if (x_periodic)
		{
			vtx.prepareFactors();
		}
		
		u += vtx.computeXInducedVelocityAt(x, y);
	}

	return u;
}

std::vector<double> DataAnalyst::computeUEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	std::vector<double> v_U;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_U.push_back(computeUAt(dm, t, x, y, x_periodic, x_period));
	}
	return v_U;
}

double DataAnalyst::computeVAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	double v(0.);
	double vtx_x, vtx_y, vtx_circ, vtx_rad;

	Vortex vtx(0., 0., 0., 0.);
	
	for (size_t i = 0; i < nb_vtx; i++)
	{
		vtx_x = dm.getXAt(time_index, i);
		vtx_y = dm.getYAt(time_index, i);
		vtx_circ = dm.getCirculationAt(time_index, i);
		vtx_rad = dm.getRegRadiusAt(time_index, i);

		vtx = Vortex(vtx_x, vtx_y, vtx_circ, vtx_rad, 0, x_periodic, x_period);

		if (x_periodic)
		{
			vtx.prepareFactors();
		}
		
		v += vtx.computeYInducedVelocityAt(x, y);
	}

	return v;
}

std::vector<double> DataAnalyst::computeVEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	std::vector<double> v_V;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_V.push_back(computeVAt(dm, t, x, y, x_periodic, x_period));
	}
	return v_V;
}

double DataAnalyst::computeVorticityAt(DataManager const &dm, size_t time_index, double x, double y, bool x_periodic, double x_period, double h)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	double u_top(computeUAt(dm, time_index, x, y + h, x_periodic, x_period));
	double u_bottom(computeUAt(dm, time_index, x, y - h, x_periodic, x_period));
	
	double v_left(computeVAt(dm, time_index, x - h, y, x_periodic, x_period));
	double v_right(computeVAt(dm, time_index, x + h, y, x_periodic, x_period));

	return 0.5*(v_right - v_left - u_top + u_bottom)/h;

}

std::vector<double> DataAnalyst::computeVorticityEvolutionAt(DataManager const &dm, double x, double y, bool x_periodic, double x_period, double h)
{
	size_t nb_steps = dm.getNbSteps();
	size_t nb_vtx = dm.getNbVtx();

	std::vector<double> v_O;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_O.push_back(computeVorticityAt(dm, t, x, y, x_periodic, x_period, h));
	}
	return v_O;
}