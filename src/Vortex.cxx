#include "Vortex.hxx"

Vortex::Vortex(double x, double y, double circulation, double regRadius, size_t fluid_id, bool X_periodic, double X_period) :
	m_x(x),
	m_y(y),
	m_circ(circulation),
	m_rad(regRadius),
	m_fluid_id(fluid_id),
	m_X_periodic(X_periodic),
	m_X_period(X_period),
	m_global_factor(0.),
	m_common_factor(0.) {}

Vortex::~Vortex() {}

double Vortex::getX() const
{
	return m_x;
}

double Vortex::getY() const
{
	return m_y;
}

double Vortex::getCirculation() const
{
	return m_circ;
}

double Vortex::getRegRadius() const
{
	return m_rad;
}

size_t Vortex::getFluidId() const
{
	return m_fluid_id;
}

bool Vortex::getXPeriodicity() const
{
	return m_X_periodic;
}

double Vortex::getXPeriod() const
{
	return m_X_period;
}


void Vortex::setX(double newX)
{
	m_x = newX;
}

void Vortex::setY(double newY)
{
	m_y = newY;
}

void Vortex::setCirculation(double newCirculation)
{
	m_circ = newCirculation;
}

void Vortex::setRegRadius(double newRegRadius)
{
	m_rad = newRegRadius;
}

void Vortex::setFluidId(size_t newId)
{
	m_fluid_id = newId;
}

void Vortex::setXPeriodicity(bool newXPeriodicity)
{
	m_X_periodic = newXPeriodicity;
}

void Vortex::setXPeriod(double newXPeriod)
{
	m_X_period = newXPeriod;
}


std::string Vortex::toString() const
{
	std::ostringstream ss;
	ss << "Vtx: fluid = " << m_fluid_id << " ; x = " << m_x << " ; y = " << m_y << " ; circulation = " << m_circ << " ; regularization radius = " << m_rad;
	return  ss.str();
}


void Vortex::prepareFactors()
{
	m_global_factor = 0.5 * m_circ / m_X_period;
	m_common_factor = 2 * M_PI / m_X_period;
}

double Vortex::computeXInducedVelocityAt(double x, double y) const
{
	double deltaX(x-m_x), deltaY(y-m_y);
	double r2(deltaX*deltaX + deltaY*deltaY);

	if (r2 < m_rad*m_rad)
	{
		return - (m_circ * deltaY) / (2 * M_PI * m_rad);
	}

	if (m_X_periodic)
	{
		double den(cosh(m_common_factor * deltaY) - cos(m_common_factor * deltaX));
		
		if (den == 0.) //Should not happen but in some edge cases...
		{
			return 0.;
		}
		
		return - m_global_factor * sinh(m_common_factor * deltaY) / den;
	}
	
	return - (m_circ * deltaY) / (2 * M_PI * r2);
}

double Vortex::computeYInducedVelocityAt(double x, double y) const
{
	double deltaX(x-m_x), deltaY(y-m_y);
	double r2(deltaX*deltaX + deltaY*deltaY);

	if (r2 < m_rad*m_rad)
	{
		return (m_circ * deltaX) / (2 * M_PI * m_rad);
	}

	if (m_X_periodic)
	{
		double den(cosh(m_common_factor * deltaY) - cos(m_common_factor * deltaX));
		
		if (den == 0.) //Should not happen but in some edge cases...
		{
			return 0.;
		}
		
		return m_global_factor * sin(m_common_factor * deltaX) / den;
	}
	
	return (m_circ * deltaX) / (2 * M_PI * r2);
}


void Vortex::move(double deltaX, double deltaY)
{
	m_x += deltaX;
	m_y += deltaY;

	if (m_X_periodic)
	{
		if (m_x < 0.)
		{
			m_x += m_X_period;
		}
		else if (m_x > m_X_period)
		{
			m_x -= m_X_period;
		}
	}
}
