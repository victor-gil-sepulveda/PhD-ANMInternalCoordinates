/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 25/09/2014
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ELASTICCONSTANTCALCULATOR_H_
#define ELASTICCONSTANTCALCULATOR_H_

class Point;

class ElasticConstantCalculator {
public:
	ElasticConstantCalculator(double k);
	virtual ~ElasticConstantCalculator();

	virtual double calculateK(Point& r_i, Point& r_j);

	double initial_constant;
};

#endif /* ELASTICCONSTANTCALCULATOR_H_ */
