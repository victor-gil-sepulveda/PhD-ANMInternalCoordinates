/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 01/10/2014
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ELASTICCALCULATORMOCK_H_
#define ELASTICCALCULATORMOCK_H_

#include "../../ElasticConstantCalculator.h"

class Point;

class ElasticCalculatorMock: public ElasticConstantCalculator {
	public:
		ElasticCalculatorMock(double k);
		virtual ~ElasticCalculatorMock();

		double calculateK(Point& r_i, Point& r_j);
};

#endif /* ELASTICCALCULATORMOCK_H_ */
