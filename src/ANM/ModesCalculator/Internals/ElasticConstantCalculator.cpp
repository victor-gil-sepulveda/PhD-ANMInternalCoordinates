/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 25/09/2014
/////////////////////////////////////////////////////////////////////////////

#include "ElasticConstantCalculator.h"
#include "../../../Tools/Math/Point.h"

ElasticConstantCalculator::ElasticConstantCalculator(double k) {
	this->initial_constant = k;
}

ElasticConstantCalculator::~ElasticConstantCalculator() {

}

double ElasticConstantCalculator::calculateK(Point& r_i, Point& r_j){
	return this->initial_constant / r_j.squaredDistance(r_i);
}


