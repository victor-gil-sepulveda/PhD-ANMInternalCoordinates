///////////////////////////////////////////////////////////
/// Point.cpp
///
/// Implementation of class Point
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////

#include "Point.h"

#include <sstream>
#include "MathTools.h"
#include "../../Tools/vectorTools.h"
#include <iomanip>

using namespace std;


Point::Point()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}
///////////////////////////////////////////////////////////////
/// \remarks
/// Class constructor.
///
/// \param x [In] Coordinate x
/// \param y [In] Coordinate y
/// \param z [In] Coordinate z
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
Point::Point(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Class constructor.
///
/// \param array [In] A 3 elements array (at least) where pos. 0 is coordinate x, pos. 1 is y
/// and pos 2 is z.
///
/// \author vgil
/// \date 10/9/2013
///////////////////////////////////////////////////////////////
Point::Point (double* array){
	setCoordinates(array);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Class copy constructor.
///
/// \param other [In] Another Point
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
Point::Point(const Point& other)
{
	*this = other;
}

Point::~Point()
{}

///////////////////////////////////////////////////////////////
/// \remarks
/// It sets the x, y and z coordinates of a Point object.
///
/// \param coords [In] Pointer to the x coordinate
///
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
void Point::setCoordinates(double *coords)
{
	this->x = *coords;
	this->y = *(coords + 1);
	this->z = *(coords + 2);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It sums the coordinates of this Point and another Point
/// and stores them in this Point, this += p
///
/// \param other [In] Another Point object.
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
void Point::add(const Point & other)
{
	double result[3];
	Math::addVectors(3, &getCoordinates()[0], &(other.getCoordinates()[0]), result);
	this->setCoordinates(result);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It substracts the coordinates of two Point objects
///
/// \param other [In] a Point object.
///
/// \author vgil
/// \author mrivero
/// \date 26/07/2011
///////////////////////////////////////////////////////////////
void Point::subtract(const Point & other)
{
	double result[3];
	Math::subtractVectors(3, &getCoordinates()[0], &(other.getCoordinates()[0]), result);
	this->setCoordinates(result);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// In place normalization to 1 of the vector represented by this point.
///
/// \return Itself.
///
/// \author vgil
///////////////////////////////////////////////////////////////
Point* Point::normalize(){
	vector<double> norm_coords = getCoordinates();
	Math::normalizeVector(&(norm_coords[0]),norm_coords.size());
	setCoordinates(&(norm_coords[0]));
	return this;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It sums the coordinates of two Point objects and produce a new Point
///
/// \param p1 [In] Point object.
/// \param p2 [In] Point object.
///
/// \return Point with the coordinates resulting from the sum of two points
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
Point Point::add(const Point & p1, const Point & p2)
{
	Point p(p1.x, p1.y, p1.z);
	p.add(p2);
	return p;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It subtracts the coordinates of two Point objects and produce a new one
///
/// \param p1 [In] Point object.
/// \param p2 [In] Point object.
///
/// \return Point with the coordinates resulting from the subtraction of two points (p-p2)
///
/// \author vgil
/// \author mrivero
/// \date 26/07/2011
///////////////////////////////////////////////////////////////
Point Point::subtract(const Point & p1, const Point & p2)
{
	Point p(p1.x, p1.y, p1.z);
	p.subtract(p2);
	return p;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It multiplies the coordinates of this Point by a scalar
/// and stores them in this Point, this *= k
///
/// \param k [In] A scalar
///
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
void Point::multiplyByScalar(double k)
{
	double result[3];
	Math::multiplyVectorByScalar(result, &getCoordinates()[0], k, 3);
	this->setCoordinates(result);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It multiplies the coordinates of a Point object by a scalar and produce a new Point
///
/// \param p [In] Point object.
/// \param k [In] A scalr
///
/// \return The resulting point object
///
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
Point Point::multiplyByScalar(const Point & p, double k)
{
	Point point(p.x, p.y, p.z);
	point.multiplyByScalar(k);
	return point;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It makes the dot product of two points
///
/// \param p1 [In] Point object.
/// \param p2 [In] Point object.
///
/// \return The dot product
///
/// \author icabeza
/// \author mrivero
/// \date 12/06/2013
///////////////////////////////////////////////////////////////
double Point::dotProduct(const Point & p1, const Point & p2)
{
	return Math::dotProduct(p1.getCoordinates(), p2.getCoordinates());
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Assignment operator
///
/// \param other [In] Another Point
///
/// \return A copy of the Point
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
Point & Point::operator=(const Point & other) {

	this->setCoordinates(&(other.getCoordinates()[0]));

	return *(this);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Equal operator.
/// Two Point objects are considered equal if they have the same coordinates
/// with a precision of 1.0e-12
///
/// \param other [In] Other Point
///
/// \return True if they have the same coordinates, false otherwise
///
/// \author mrivero
/// \date 16/01/2013
///////////////////////////////////////////////////////////////
bool Point::operator==(const Point& other) const
{
	const double EQUAL_PRECISION = 1.0e-12;

	return Math::areEqual(this->x, other.x, EQUAL_PRECISION) and
			Math::areEqual(this->y, other.y, EQUAL_PRECISION) and
			Math::areEqual(this->z, other.z, EQUAL_PRECISION);
}

double Point::getX() const
{
	return x;
}

double Point::getY() const
{
	return y;
}

double Point::getZ() const
{
	return z;
}

string Point::toString() const
{
	ostringstream s;
	s <<setprecision(16)<<"("<< x<<string(", ")<<y<<string(", ")<<z<<")";
	return s.str();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It returns the distance between this Point object and another Point object
///
/// \param other [In] Another point object
///
/// \return The distance between this Point object and another Point object
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////////
double Point::distance(const Point & other) const
{
	return Math::distance(this->x, this->y, this->z, other.x, other.y, other.z);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It returns the squared distance between this Point object and another Point object
///
/// \param other [In] Another point object
///
/// \return The squared distance
///
/// \author vgil
/// \date 25/06/2011
///////////////////////////////////////////////////////////////
#include<iostream>
double Point::squaredDistance(const Point& other) const
{
//	std::cout << "DBG center: "<< this->x << " " << this->y << " " << this->z << " atom:" << other.x << " " << other.y << " "<< other.z ;
	double distance =  Math::squaredDistance(this->x, this->y, this->z, other.x, other.y, other.z);
//	cout<< " distance: "<< distance <<endl;
	return distance;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It returns a copy of the coordinates of the Point
///
/// \return A vector with a copy of the coordinates of the Point
///
/// \author mrivero
/// \date 29/06/2012
///////////////////////////////////////////////////////////////
std::vector<double> Point::getCoordinates() const
{
	vector<double> coords;

	coords.push_back(this->x);
	coords.push_back(this->y);
	coords.push_back(this->z);

	return coords;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It computes the Points module (considering it as a vector of three components)
///
/// \return The module
///
/// \author mrivero
/// \date 27/05/2013
///////////////////////////////////////////////////////////////
double Point::module()
{
	return Math::euclideanNorm(3, &(this->getCoordinates()[0]));
}
