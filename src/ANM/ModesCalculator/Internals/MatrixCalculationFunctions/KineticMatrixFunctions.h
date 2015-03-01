/////////////////////////////////////////////////////////////////////////////
/// \file KineticMatrixFunctions.h
///
/// \brief Functions used in the computation of the kinetic matrix
///
/// \author vgil
/// \author arincon
/// \date 11/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef KINETICMATRIXFUNCTIONS_H_
#define KINETICMATRIXFUNCTIONS_H_
#include <vector>
#include <utility>
#include "TriangularMatrices.h"
#include "TensorCalcTypes.h"
#include "../../../../Tools/Math/Point.h"

class Atom;
class Unit;
class CenterOfMass;

namespace ANMICKineticMatrixCalculator {

	// I Calculation
	void calculateI(double I[3][3],
			std::vector<Unit*>& units,
			const std::pair<unsigned int, unsigned int>& range,
			ICalcType calc_type);

	void calculateC(const Point& p, double result[3][3]);
	void calculateP(const Point& p, double result[3][3]);

	// K calculation
	double calculateM(std::vector<Unit*>& units, const std::pair<int,int>& range,bool skip_OXT);

	CenterOfMass calculateMobileBodyCOM(std::vector<Unit*>& units,
			const std::pair<int,int>& range,
			bool skip_OXT);

	double calculateKab( int alpha_dihedral,
							int beta_dihedral,
							Point* r_a,
							Point* r_b,
							Point* e_a,
							Point* e_b,
							std::vector<Unit*>& units,
							double M,
							double const I[3][3],
							double const I_inv[3][3],
							ICalcType tensor_calc_type);

	TriangularMatrix* calculateK( std::vector<Unit*>& units,
							ICalcType tensor_calc_type);

	// Jacobian
	void Jacobi(std::vector<Unit*>& units);

	// Derivatives
	void  dri_dq(std::vector<Unit*>& units,
			unsigned int alpha,
			double const I_inv[3][3],
			std::vector<Point>& term1_v,
			std::vector<Point>& term2_v);

	void  dri_dq_2(std::vector<Unit*>& units,
				unsigned int beta,
				double const I_inv[3][3],
				std::vector<Point>& term1_v,
				std::vector<Point>& term2_v);
};


#endif /* KINETICMATRIXFUNCTIONS_H_ */
