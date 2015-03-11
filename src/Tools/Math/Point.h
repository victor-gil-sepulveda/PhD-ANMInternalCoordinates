///////////////////////////////////////////////////////////
/// \file Point.h
///
/// \brief This class models a 3D point
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
///////////////////////////////////////////////////////////

#pragma once

#ifndef POINT_H_
#define POINT_H_

#include <iosfwd>
#include <vector>

/////////////////////////////////////////////////////////////////////////////
/// \brief This class models a 3D point
///
/// \author vgil
/// \author mrivero
/// \date 19/10/2010
//////////////////////////////////////////////////////////////////////////////
class Point
{
	public:
		Point();
		Point(double x, double y, double z);
		Point (double*);
		Point(const Point& pt);
		virtual ~Point();

		// Getters
		double getX() const;
		double getY() const;
		double getZ() const;
		std::vector<double> getCoordinates() const;

		// Other methods
		void add(const Point & p);
		void subtract(const Point & p);
		Point* normalize();
		void multiplyByScalar(double k);
		double distance(const Point & p) const;
		double squaredDistance(const Point& p) const;
		std::string toString() const;
		double module();

		// Operators
		Point & operator=(const Point & that);
		bool operator==(const Point& other) const;

		// Static methods
		static Point add(const Point & p1, const Point & p2);
		static Point subtract(const Point & p1, const Point & p2);
		static Point multiplyByScalar(const Point & p, double k);
		static double dotProduct(const Point & p1, const Point & p2);


	private:
		double x, y, z;
		void setCoordinates(double *startOfAtomInCoordinates);

};

#endif /* POINT_H_ */
