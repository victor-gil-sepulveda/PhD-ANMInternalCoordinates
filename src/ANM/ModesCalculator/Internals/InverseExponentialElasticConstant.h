/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 15/10/2014
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef INVERSEEXPONENTIALELASTICCONSTANT_H_
#define INVERSEEXPONENTIALELASTICCONSTANT_H_
#include "ElasticConstantCalculator.h"

class InverseExponentialElasticConstant :public ElasticConstantCalculator {
	public:
		InverseExponentialElasticConstant(double k, double x0, int power);
		virtual ~InverseExponentialElasticConstant();

		double calculateK(Point& r_i, Point& r_j);

	double x0;
	int power;

};

#endif /* INVERSEEXPONENTIALELASTICCONSTANT_H_ */
