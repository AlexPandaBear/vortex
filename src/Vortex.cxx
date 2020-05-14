#include "Vortex.hxx"

Vortex::Vortex(double x, double y, double circulation, double regRadius, size_t fluid_id) :
	m_x(x),
	m_y(y),
	m_circ(circulation),
	m_rad(regRadius),
	m_fluid_id(fluid_id) {}

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


std::string Vortex::toString() const
{
	std::ostringstream ss;
	ss << "Vtx: fluid = " << m_fluid_id << " ; x = " << m_x << " ; y = " << m_y << " ; circulation = " << m_circ << " ; regularization radius = " << m_rad;
	return  ss.str();
}


double Vortex::computeXInducedVelocityAt(double x, double y) const
{
	return computeXInducedVelocityAt(x, y, m_x, m_y, m_circ, m_rad);
}

double Vortex::computeYInducedVelocityAt(double x, double y) const
{
	return computeYInducedVelocityAt(x, y, m_x, m_y, m_circ, m_rad);
}

double Vortex::computeXInducedVelocityAt(double x, double y, double vtx_x, double vtx_y, double vtx_circ, double vtx_rad)
{
	double r2((x-vtx_x)*(x-vtx_x) + (y-vtx_y)*(y-vtx_y));
	if (r2 < vtx_rad*vtx_rad)
	{
		return -(vtx_circ*(y-vtx_y))/(2*M_PI*vtx_rad);
	}
	return -(vtx_circ*(y-vtx_y))/(2*M_PI*r2);
}

double Vortex::computeYInducedVelocityAt(double x, double y, double vtx_x, double vtx_y, double vtx_circ, double vtx_rad)
{
	double r2((x-vtx_x)*(x-vtx_x) + (y-vtx_y)*(y-vtx_y));
	if (r2 < vtx_rad*vtx_rad)
	{
		return (vtx_circ*(x-vtx_x))/(2*M_PI*vtx_rad);
	}
	return (vtx_circ*(x-vtx_x))/(2*M_PI*r2);
}


void Vortex::move(double deltaX, double deltaY)
{
	m_x += deltaX;
	m_y += deltaY;
}
