/////////////////////////////////////////////////////////////////////////////
/// \file HarmonicDihedralConstraintTerm.h
/// \brief Definition of the HarmonicDihedralConstraintTerm class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author arincon
/// \date 16/08/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef HARMONICDIHEDRALCONSTRAINTTERM_H_
#define HARMONICDIHEDRALCONSTRAINTTERM_H_

#include "../ConstraintTerm.h"
#include <vector>

/////////////////////////////////////////////////////////////////////////////
/// \brief HarmonicDihedralConstraintTerm class
///
/// This class models a harmonic dihedral constraint term.
///
/// \author arincon
/// \date 16/08/2012
/////////////////////////////////////////////////////////////////////////////
class HarmonicDihedralConstraintTerm : public ConstraintTerm
{
	public:
		HarmonicDihedralConstraintTerm(double springConstant, double equilibriumAngle);
		HarmonicDihedralConstraintTerm(double springConstant, double equilibriumAngle, std::vector<Atom*> &atoms);
		virtual ~HarmonicDihedralConstraintTerm();

		// Implementation of the ConstraintTerm interface
		double getCurrentAngle(const double * const coords);
		double getEnergy(const double * const coords);
		double getEnergyGradient(const double * const coords, double * grad);
		double getEnergyGradientHess(const double * const coords, double *grad, SparseMatrixCRS *hessian);
		void 	getHess(const double * const coords, SparseMatrixCRS *hessian);

		std::string generateReport() const { return ""; };

	private:

		double springConstant;
		double equilibriumAngle;

		std::vector<Atom*> atoms;

		static void generateFullHessianMatrixFromTriangularMatrix(std::vector<double> &triangularMatrix, double fullMatrix[][12]);

		friend class TestHarmonicDihedralConstraintTerm;
};

#endif /* HARMONICDIHEDRALCONSTRAINTTERM_H_ */
