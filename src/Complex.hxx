#pragma once

#include <cmath>
#include <sstream>
#include <exception>
#include <string>

/**
 *A C++ class representing complex numbers
 */
class Complex
{
private:
	double m_real;
	double m_imag;

public:
	/**
	 * The constructor of the class
	 *
	 * @param x the real part of the complex number
	 *	 
	 * @param y the imaginary part of the complex number
	 * 
	 * @param mode optional parameter allowing to create complex numbers using x as the module and y as the argument if equal to "polar"
	 */
	Complex(double x = 0., double y = 0., std::string mode = "cartesian");
	
	/**
	 *The destructor of the class
	 */
	~Complex();

	/**
	 * Getter to the real part of the instance
	 * 
	 * @returns The real part as a double
	 */
	double getReal() const;
	
	/**
	 * Getter to the imaginary part of the instance
	 * 
	 * @returns The imaginary part as a double
	 */
	double getImag() const;
	
	/**
	 * Getter to the module of the instance
	 * 
	 * @returns The module as a double
	 */
	double getMag() const;
	
	/**
	 * Getter to the square of the module of the instance (more efficient than using getMag())
	 * 
	 * @returns The module squared as a double
	 */
	double getMagSquared() const;
	
	/**
	 * Getter of the argument of the instance
	 * 
	 * @returns The argument in radians as a double
	 * 
	 * @throws std::invalid_argument Thrown if the real or the imaginary part is NaN
	 */
	double getArg() const;

	/**
	 * Setter of the real part of the instance
	 * 
	 * @param x The new value for the real part
	 */
	void setReal(double x);
	
	/**
	 * Setter of the imaginary part of the instance
	 * 
	 * @param y The new value for the imaginary part
	 */
	void setImag(double y);

	/**
	 * Method returning a printable image of the complex number
	 * 
	 * @returns A std::string showing the current state of the instance
	 */
	std::string toString() const;
	
	/**
	 * Method returning a condensed and minimal image of the complex number to save in a file
	 * 
	 * @returns A std::string in the format "<real part> <imaginary part>"
	 */
	std::string toFileFormat() const;

	/**
	 * Operator checking the equality of the instance to another
	 * 
	 * @returns The boolean true if both real and imaginary parts are equal, and false otherwise
	 */
	bool operator==(Complex const &z) const;
	
	/**
	 * Operator checking the non-equality of the instance to another
	 * 
	 * @returns The boolean opposite to the operator ==
	 */
	bool operator!=(Complex const &z) const;

	/**
	 * Method computing the inverse of the instance
	 * 
	 * @returns Another instance representing the inverse of the instance calling the method
	 * 
	 * @throws std::invalid_argument Thrown if the instance calling the method is equal to zero
	 */
	Complex inverse() const;
	
	/**
	 * Method computing the complex conjugate of the instance
	 * 
	 * @returns Another instance representing the complex conjugate of the instance calling the method
	 */
	Complex conj() const;

	/**
	 * Operator adding to the instance the value of another
	 * 
	 * @param z Another instance of the class representing the value to add
	 */
	void operator+=(Complex const &z);
	
	/**
	 * Operator substracting to the instance the value of another
	 * 
	 * @param z Another instance of the class representing the value to substract
	 */
	void operator-=(Complex const &z);
	
	/**
	 * Operator multiplying the instance by the value of another
	 * 
	 * @param z Another instance of the class representing the value to multiply by
	 */
	void operator*=(Complex const &z);
	
	/**
	 * Operator dividing the instance the value of another
	 * 
	 * @param z Another instance of the class representing the value to divide by
	 * 
	 * @throws std::invalid_argument Thrown if the parameter z is equal to zero
	 */
	void operator/=(Complex const &z);

	/**
	 * Operator returning the result of the addition of the instance with another
	 * 
	 * @param z Another instance of the class representing the value to add
	 */
	Complex operator+(Complex const &z) const;
	
	/**
	 * Operator returning the result of the substraction of the instance with another
	 * 
	 * @param z Another instance of the class representing the value to substract
	 */
	Complex operator-(Complex const &z) const;
	
	/**
	 * Operator returning the result of the multiplication of the instance by another
	 * 
	 * @param z Another instance of the class representing the value to multiply by
	 */
	Complex operator*(Complex const &z) const;
	
	/**
	 * Operator returning the result of the division of the instance by another
	 * 
	 * @param z Another instance of the class representing the value to divide by
	 */
	Complex operator/(Complex const &z) const;

	/**
	 * Operator returning the opposite of the instance
	 */
	Complex operator-() const;
	
	/**
	 * Operator returning the result of the multiplication of the instance with a double
	 * 
	 * @param x The value to multiply by
	 */
	Complex operator*(double x) const;

	/**
	 * Static method returning the complex exponential of an instance
	 * 
	 * @param z The instance on which computation will be done
	 * 
	 * @returns Another instance representing the value of exp(z)
	 */
	static Complex exp(Complex const &z);
	
	/**
	 * Static method returning the complex logarithm of an instance
	 * 
	 * @param z The instance on which computation will be done
	 * 
	 * @returns Another instance representing the value of ln(z)
	 */
	static Complex ln(Complex const &z);

	/**
	 * Static method returning the dot product of two instances representing vectors of the plane
	 * 
	 * @param z1 The first vector used to compute the dot product
	 * 
	 * @param z1 The second vector used to compute the dot product
	 * 
	 * @returns The dot product as a double
	 */
	static double dotProduct(Complex const &z1, Complex const &z2);
	
	/**
	 * Static method returning the cross product of two instances representing vectors of the plane
	 * 
	 * @param z1 The first vector used to compute the cross product
	 * 
	 * @param z1 The second vector used to compute the cross product
	 * 
	 * @returns The norm of the cross product as a double
	 */
	static double crossProduct(Complex const &z1, Complex const &z2);
};