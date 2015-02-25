/////////////////////////////////////////////////////////////////////////////
/// \file AnmNormalizer.h
///
/// \brief Different functions to normalize eigenvectors
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author atarraco
/// \author mrivero
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMNORMALIZER_H_
#define ANMNORMALIZER_H_
#include <vector>

/////////////////////////////////////////////////////////////////////////////
/// \brief Different functions to normalize eigenvectors
///
/// \author atarraco
/// \author mrivero
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////
namespace AnmNormalizer {
	void normalizeByInverseLargestNorm(std::vector<double> & vector);
	void normalizeByInverseLargestNormWithThermalScaling(std::vector<double> & vector, double eigenValue, double constantForHessian);
	void normalizeByLargestValue(std::vector<double> & vector);
}

#endif /* ANMNORMALIZER_H_ */
