/////////////////////////////////////////////////////////////////////////////
/// \file AnmEigen.h
///
/// \brief This class holds the eigenvalues and eigenvectors
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \author atarraco
/// \author xoro
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMEIGEN_H_
#define ANMEIGEN_H_

#include <vector>
#include <iosfwd>

/////////////////////////////////////////////////////////////////////////////
/// \brief This class holds the eigenvalues and eigenvectors
///
/// \author mrivero
/// \author atarraco
/// \author xoro
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////
class AnmEigen
{
	public:
		AnmEigen();
		virtual ~AnmEigen();

		void computeAverageEigenVector(const std::vector<double> & weights, std::vector<double> & directions);
		void initialize(const double * const eigenValues, const double * const eigenVectors,
				unsigned int numberOfModes, unsigned int numberOfNodes, bool usingCartesian=true);
		void initialize(std::vector<double>& eigenValues,
						std::vector<std::vector<double> >& eigenVectors,
						bool usingCartesian);

		void normalizeByInverseLargestNorm();
		void normalizeByInverseNormWithThermal(double constantForHessian);
		void normalizeByLargestValue();

		unsigned int getNumberOfModes() const;
		unsigned int getEigenVectorsDimension() const;

		std::vector<double> & getEigenVectorOfMode(unsigned int mode);
		double getEigenValueOfMode(unsigned int mode);

		std::string toString() const;

		AnmEigen& operator=(const AnmEigen&);

		unsigned int numberOfNodes;
		unsigned int numberOfModes;
		std::vector<double> values;
		std::vector<std::vector<double> > vectors;
};

#endif /* ANMEIGEN_H_ */
