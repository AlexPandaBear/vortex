#include "Complex.hxx"

Complex::Complex(double x, double y, std::string mode) :
	m_real(x),
	m_imag(y)
{
	if (mode == "cartesian")
	{}
	else if (mode == "polar")
	{
		m_real = x*cos(y);
		m_imag = x*sin(y);
	}
	else
	{
		throw std::invalid_argument("Can't build complex number in mode " + mode);
	}
}

Complex::~Complex() {}

double Complex::getReal() const
{
	return m_real;
}

double Complex::getImag() const
{
	return m_imag;
}

double Complex::getMag() const
{
	return sqrt(this->getMagSquared());
}

double Complex::getMagSquared() const
{
	return m_real*m_real + m_imag*m_imag;
}

double Complex::getArg() const
{
	if (m_real == 0.)
	{
		if (m_imag == 0.)
		{
			return 0.; 
		}
		else if (m_imag > 0.)
		{
			return M_PI/2.;
		}
		else if (m_imag < 0.)
		{
			return -M_PI/2.;
		}
		throw std::invalid_argument("Imaginary part is NaN");
	}

	else if (m_real > 0.)
	{
		return atan(m_imag / m_real);
	}

	else if (m_real < 0.)
	{
		if (m_imag >= 0.)
		{
			return M_PI + atan(m_imag / m_real);
		}
		else if (m_imag < 0.)
		{
			return - M_PI + atan(m_imag / m_real);
		}
		throw std::invalid_argument("Imaginary part is NaN");
	}
	
	throw std::invalid_argument("Real part is NaN");
}

void Complex::setReal(double x)
{
	m_real = x;
}

void Complex::setImag(double y)
{
	m_imag = y;
}

std::string Complex::toString() const
{
	std::ostringstream ss;
	ss << m_real << " + " << m_imag << "i";
	return ss.str();
}

std::string Complex::toFileFormat() const
{
	std::ostringstream ss;
	ss << m_real << " " << m_imag;
	return ss.str();
}

bool Complex::operator==(Complex const &z) const
{
	if (m_real == z.getReal() && m_imag == z.getImag())
	{
		return true;
	}
	return false;
}

bool Complex::operator!=(Complex const &z) const
{
	return !(*this == z);
}

Complex Complex::inverse() const
{
	double norm2(this->getMagSquared());
	if (norm2 == 0.)
	{
		throw std::invalid_argument("Can't divide by zero");
	}
	return Complex(m_real/norm2, -m_imag/norm2);
}

Complex Complex::conj() const
{
	return Complex(m_real, -m_imag);
}

void Complex::operator+=(Complex const &z)
{
	m_real += z.getReal();
	m_imag += z.getImag();
}

void Complex::operator-=(Complex const &z)
{
	m_real -= z.getReal();
	m_imag -= z.getImag();
}

void Complex::operator*=(Complex const &z)
{
	m_real = m_real*z.getReal() - m_imag*z.getImag();
	m_imag = m_real*z.getImag() + m_imag*z.getReal();
}

void Complex::operator/=(Complex const &z)
{
	Complex inv(z.inverse());
	m_real = m_real*inv.getReal() - m_imag*inv.getImag();
	m_imag = m_real*inv.getImag() + m_imag*inv.getReal();
}

Complex Complex::operator+(Complex const &z) const
{
	return Complex(m_real + z.getReal(), m_imag + z.getImag());
}

Complex Complex::operator-(Complex const &z) const
{
	return Complex(m_real - z.getReal(), m_imag - z.getImag());
}

Complex Complex::operator*(Complex const &z) const
{
	return Complex(m_real*z.getReal() - m_imag*z.getImag(), m_real*z.getImag() + m_imag*z.getReal());
}

Complex Complex::operator/(Complex const &z) const
{
	Complex res(*this);
	return res*z.inverse();
}

Complex Complex::operator-() const
{
	return Complex(-m_real, -m_imag);
}

Complex Complex::operator*(double x) const
{
	return Complex(x*m_real, x*m_imag);
}

Complex Complex::exp(Complex const &z)
{
	double mod(std::exp(z.getReal()));
	return Complex(mod*cos(z.getImag()), mod*sin(z.getImag()));
}

Complex Complex::ln(Complex const &z)
{
	double mod(z.getMag());
	if (mod == 0.)
	{
		throw std::invalid_argument("Can't take log of zero");
	}
	return Complex(log(mod), z.getArg());
}

double Complex::dotProduct(Complex const &z1, Complex const &z2)
{
	return z1.getReal()*z2.getReal() + z1.getImag()*z2.getImag();
}

double Complex::crossProduct(Complex const &z1, Complex const &z2)
{
	return z1.getReal()*z2.getImag() - z1.getImag()*z2.getReal();
}
