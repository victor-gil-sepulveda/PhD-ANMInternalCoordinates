/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 01/10/2014
/////////////////////////////////////////////////////////////////////////////

#include "ElasticCalculatorMock.h"
#include <iostream>
#include "../../../../../Tools/Math/Point.h"
using namespace std;

ElasticCalculatorMock::ElasticCalculatorMock(double k):ElasticConstantCalculator(k) {
}

ElasticCalculatorMock::~ElasticCalculatorMock() {
}

double ElasticCalculatorMock::calculateK(Point& r_i, Point& r_j){
	return this->initial_constant * r_i.squaredDistance(r_j);//It will be divided later by sq dist
}
