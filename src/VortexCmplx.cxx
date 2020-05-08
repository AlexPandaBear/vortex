#include "VortexCmplx.hxx"

VortexCmplx::VortexCmplx(Complex const &position, double circulation, double regRadius) :
	m_pos(position),
	m_circ(circulation),
	m_rad(regRadius) {}

VortexCmplx::VortexCmplx(double x, double y, double circulation, double regRadius, std::string mode) :
	m_pos(x, y, mode),
	m_circ(circulation),
	m_rad(regRadius) {}

VortexCmplx::~VortexCmplx() {}

Complex VortexCmplx::getPosition() const
{
	return m_pos;
}

double VortexCmplx::getCirculation() const
{
	return m_circ;
}

double VortexCmplx::getRegRadius() const
{
	return m_rad;
}

void VortexCmplx::setPosition(Complex newPosition)
{
	m_pos = newPosition;
}

void VortexCmplx::setCirculation(double newCirculation)
{
	m_circ = newCirculation;
}

void VortexCmplx::setRegRadius(double newRegRadius)
{
	m_rad = newRegRadius;
}

std::string VortexCmplx::toString() const
{
	std::ostringstream ss;
	ss << "Vtx: position = " << m_pos.toString() << " ; circulation = " << m_circ << " ; regularization radius = " << m_rad;
	return  ss.str();
}

Complex VortexCmplx::computeInducedVelocityAt(Complex position) const
{
	return computeInducedVelocityAt(position, m_pos, m_circ, m_rad);
}

Complex VortexCmplx::computeInducedVelocityAt(Complex position, Complex vtx_pos, double vtx_circ, double vtx_rad)
{
	Complex i(0,1);
	Complex dist(position - vtx_pos);

	if (dist.getMag() < vtx_rad)
	{
		return  i*dist*vtx_circ / (2*M_PI*vtx_rad*vtx_rad);
	}

	return (-(i/(dist)) * (vtx_circ/(2*M_PI))).conj();
}

void VortexCmplx::move(Complex deltaPosition)
{
	m_pos += deltaPosition;
}