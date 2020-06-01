#include "DataAnalyst.hxx"

DataAnalyst::DataAnalyst(DataManager const &dm) :
	ptr_data(&dm) {}

DataAnalyst::~DataAnalyst() {}

void DataAnalyst::resetDataPtr(DataManager const &dm)
{
	ptr_data = &dm;
}

double DataAnalyst::computeHamiltonianAt(size_t time_index, bool x_periodic, double x_period) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

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
				deltaX = ptr_data->getXAt(time_index, i) - ptr_data->getXAt(time_index, j);
				deltaY = ptr_data->getYAt(time_index, i) - ptr_data->getYAt(time_index, j);
				H += ptr_data->getCirculationAt(time_index, i)*ptr_data->getCirculationAt(time_index, j) * log(cosh(factor*deltaY) - cos(factor*deltaX));
			}
		}
	}

	else
	{
		for (size_t i = 0; i < nb_vtx; i++)
		{
			for (size_t j = i+1; j < nb_vtx; j++)
			{
				deltaX = ptr_data->getXAt(time_index, i) - ptr_data->getXAt(time_index, j);
				deltaY = ptr_data->getYAt(time_index, i) - ptr_data->getYAt(time_index, j);
				H += ptr_data->getCirculationAt(time_index, i)*ptr_data->getCirculationAt(time_index, j)*0.5*log(deltaX*deltaX + deltaY*deltaY);
			}
		}
	}

	return H/(4*M_PI);
}

void DataAnalyst::computeHamiltonianEvolutionBetween(size_t start_step, size_t end_step, std::vector<double> &v_H, bool x_periodic, double x_period) const
{
	for (size_t t = start_step; t < end_step; t++)
	{
		v_H[t] = computeHamiltonianAt(t, x_periodic, x_period);
	}
}

std::vector<double> DataAnalyst::computeHamiltonianEvolution(size_t nb_threads, bool x_periodic, double x_period) const
{
	if (nb_threads < 1)
	{
		std::ostringstream ss;
		ss << "Unable to perform computation on " << nb_threads << " threads";
		throw std::invalid_argument(ss.str());
	}

	size_t nb_steps = ptr_data->getNbSteps();
	std::vector<double> v_H(nb_steps+1, 0.);
	
	std::vector<std::thread> v_threads;
	size_t start_step(0), end_step(nb_steps/nb_threads);

	for (size_t i = 0; i < nb_threads-1; i++)
	{
		v_threads.push_back(std::thread(&DataAnalyst::computeHamiltonianEvolutionBetween, this, start_step, end_step, std::ref(v_H), x_periodic, x_period));
		start_step = end_step;
		end_step += nb_steps/nb_threads;
	}

	end_step = nb_steps+1;
	computeHamiltonianEvolutionBetween(start_step, end_step, v_H, x_periodic, x_period);

	for (auto &th : v_threads) {th.join();}

	return v_H;
}

std::map<size_t, double> DataAnalyst::computeCompositionAt(double x, double y, size_t step, double radius) const
{
	size_t nb_vtx = ptr_data->getNbVtx();
	size_t fluid_id;
	std::map<size_t, double> map_compo;

	for (size_t i = 0; i < nb_vtx; i++)
	{
		fluid_id = ptr_data->getFluidId(i);
		if (map_compo.find(fluid_id) == map_compo.end())
		{
			map_compo[fluid_id] = 0.;
		}
	}

	double deltaX, deltaY;
	double dist, dist2, rad2(radius*radius);

	for (size_t i = 0; i < nb_vtx; i++)
	{
		deltaX = ptr_data->getXAt(step, i) - x;
		deltaY = ptr_data->getYAt(step, i) - y;
		dist2 = deltaX*deltaX + deltaY*deltaY;
		
		if (dist2 < rad2)
		{
			dist = sqrt(dist2);
			fluid_id = ptr_data->getFluidId(i);
			map_compo[fluid_id] += 1. - dist/radius;
		}
	}

	double sum(0.);
	for (std::map<size_t, double>::iterator it = map_compo.begin(); it != map_compo.end(); it++) {sum += it->second;}
	for (std::map<size_t, double>::iterator it = map_compo.begin(); it != map_compo.end(); it++) {it->second /= sum;}

	return map_compo;
}

std::vector<std::map<size_t, double>> DataAnalyst::computeCompositionEvolutionAt(double x, double y, double radius) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	std::vector<std::map<size_t, double>> v_C;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_C.push_back(computeCompositionAt(x, y, t, radius));
	}
	return v_C;
}

double DataAnalyst::computeUAt(size_t time_index, double x, double y, bool x_periodic, double x_period) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	double u(0.);
	double vtx_x, vtx_y, vtx_circ, vtx_rad;

	Vortex vtx(0., 0., 0., 0.);
	
	for (size_t i = 0; i < nb_vtx; i++)
	{
		vtx_x = ptr_data->getXAt(time_index, i);
		vtx_y = ptr_data->getYAt(time_index, i);
		vtx_circ = ptr_data->getCirculationAt(time_index, i);
		vtx_rad = ptr_data->getRegRadiusAt(time_index, i);

		vtx = Vortex(vtx_x, vtx_y, vtx_circ, vtx_rad, 0, x_periodic, x_period);

		if (x_periodic)
		{
			vtx.prepareFactors();
		}
		
		u += vtx.computeXInducedVelocityAt(x, y);
	}

	return u;
}

std::vector<double> DataAnalyst::computeUEvolutionAt(double x, double y, bool x_periodic, double x_period) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	std::vector<double> v_U;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_U.push_back(computeUAt(t, x, y, x_periodic, x_period));
	}
	return v_U;
}

double DataAnalyst::computeVAt(size_t time_index, double x, double y, bool x_periodic, double x_period) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	double v(0.);
	double vtx_x, vtx_y, vtx_circ, vtx_rad;

	Vortex vtx(0., 0., 0., 0.);
	
	for (size_t i = 0; i < nb_vtx; i++)
	{
		vtx_x = ptr_data->getXAt(time_index, i);
		vtx_y = ptr_data->getYAt(time_index, i);
		vtx_circ = ptr_data->getCirculationAt(time_index, i);
		vtx_rad = ptr_data->getRegRadiusAt(time_index, i);

		vtx = Vortex(vtx_x, vtx_y, vtx_circ, vtx_rad, 0, x_periodic, x_period);

		if (x_periodic)
		{
			vtx.prepareFactors();
		}
		
		v += vtx.computeYInducedVelocityAt(x, y);
	}

	return v;
}

std::vector<double> DataAnalyst::computeVEvolutionAt(double x, double y, bool x_periodic, double x_period) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	std::vector<double> v_V;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_V.push_back(computeVAt(t, x, y, x_periodic, x_period));
	}
	return v_V;
}

double DataAnalyst::computeVorticityAt(size_t time_index, double x, double y, bool x_periodic, double x_period, double h) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	double u_top(computeUAt(time_index, x, y + h, x_periodic, x_period));
	double u_bottom(computeUAt(time_index, x, y - h, x_periodic, x_period));
	
	double v_left(computeVAt(time_index, x - h, y, x_periodic, x_period));
	double v_right(computeVAt(time_index, x + h, y, x_periodic, x_period));

	return 0.5*(v_right - v_left - u_top + u_bottom)/h;

}

std::vector<double> DataAnalyst::computeVorticityEvolutionAt(double x, double y, bool x_periodic, double x_period, double h) const
{
	size_t nb_steps = ptr_data->getNbSteps();
	size_t nb_vtx = ptr_data->getNbVtx();

	std::vector<double> v_O;
	for (size_t t = 0; t < nb_steps+1; t++)
	{
		v_O.push_back(computeVorticityAt(t, x, y, x_periodic, x_period, h));
	}
	return v_O;
}