/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 15/10/2014
/////////////////////////////////////////////////////////////////////////////

#include "InverseExponentialElasticConstant.h"
#include "../../../Tools/Math/Point.h"
#include <cmath>

InverseExponentialElasticConstant::InverseExponentialElasticConstant(double k,
		double x0, int power):ElasticConstantCalculator(k){
	this->x0 = x0;
	this->power = power;
}

InverseExponentialElasticConstant::~InverseExponentialElasticConstant() {

}

double InverseExponentialElasticConstant::calculateK(Point& r_i, Point& r_j){
	// initial values: k = 1, x0 = 3.8, power = 6, cutoff = 10

	double d = r_j.distance(r_i);
	double k0 = initial_constant;

	return k0 / ( 1.0 + pow( d / x0, power ) );
//	return (k0 / ( 1.0 + pow( d / x0, power ) ))/(d*d);
}
