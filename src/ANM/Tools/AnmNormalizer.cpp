/////////////////////////////////////////////////////////////////////////////
/// AnmNormalizer.cpp
///
/// Implementation of AnmNormalizer class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////

#include "AnmNormalizer.h"
#include "../../Tools/Math/MathTools.h"
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// This function normalizes the vector by its inverse largest norm
///
/// \param eigenVector [In/Out] Eigenvector to be normalized
///
/// \author atarraco
/// \author mrivero
/// \date 03/09/2012
///////////////////////////////////////////////////////////////
void AnmNormalizer::normalizeByInverseLargestNorm(std::vector<double> & eigenVector)
{
	unsigned int eigenVectorsDimension = eigenVector.size();
	double factorNormalization = Math::inverseLargestNorm(eigenVectorsDimension, &eigenVector[0]);
	Math::multiplyVectorByScalar(&eigenVector[0], factorNormalization, eigenVectorsDimension);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function normalizes the eigen vector by its inverse largest norm
/// with thermal scaling
///
/// \param eigenValue [In] Corresponding eigenValue
/// \param constantForHessian [In] Parameter used to create the hessian
///
/// \param eigenVector [In/Out] Eigenvector to be normalized
///
/// \author atarraco
/// \author mrivero
/// \date 03/09/2012
///////////////////////////////////////////////////////////////
void AnmNormalizer::normalizeByInverseLargestNormWithThermalScaling(std::vector<double> & eigenVector, double eigenValue, double constantForHessian)
{
	unsigned int eigenVectorsDimension = eigenVector.size();

	double factorNormalization = Math::inverseEuclideanNorm(eigenVectorsDimension, &eigenVector[0])
																* 1.0 / sqrt(constantForHessian * eigenValue);

	Math::multiplyVectorByScalar(&eigenVector[0], factorNormalization, eigenVectorsDimension);
}


///////////////////////////////////////////////////////////////
/// \remarks
/// Scales a vector by dividing all its elements by its biggest element value.
///
/// \param eigenVector [In/Out] vector (eigenvector) to be scaled
///
/// \author vgil
/// \date 7/01/2015
///////////////////////////////////////////////////////////////
void AnmNormalizer::normalizeByLargestValue(std::vector<double> & eigenVector){
	unsigned int eigenVectorsDimension = eigenVector.size();
	double max_element = std::abs(*(std::max_element(eigenVector.begin(), eigenVector.end())));
	double min_element = std::abs(*(std::min_element(eigenVector.begin(), eigenVector.end())));
	double largest = std::max(min_element, max_element);
	double factorNormalization = 1. / largest;
	cout<<"DBG: largest value for scaling:"<<largest<<endl;
	Math::multiplyVectorByScalar(&eigenVector[0], factorNormalization, eigenVectorsDimension);
}
