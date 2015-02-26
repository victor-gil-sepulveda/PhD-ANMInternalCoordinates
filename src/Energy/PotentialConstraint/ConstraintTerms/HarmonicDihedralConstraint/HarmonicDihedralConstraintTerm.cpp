///////////////////////////////////////////////////////////
/// HarmonicDihedralConstraintTerm.cpp
///
/// Implementation of the HarmonicDihedralConstraintTerm class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author arincon
/// \date 16/08/2012
///////////////////////////////////////////////////////////

#include "HarmonicDihedralConstraintTerm.h"

#include "../../../../Molecules/Atom.h"
#include "../../../../Tools/Math/SparseMatrixes/SparseMatrixCRS.h"
#include "HarmonicDihedralConstraintFunctions/HarmonicDihedralConstraintFunctions.h"
#include <vector>
#include "../../../TermCalculators/Utils/GradientHessianDihedralUpdater.h"
#include "../../../../Tools/vectorTools.h"
using namespace std;

HarmonicDihedralConstraintTerm::HarmonicDihedralConstraintTerm(double springConstant,
		double equilibriumAngle, std::vector<Atom*> &nodes)
{
	this->springConstant = springConstant;
	this->equilibriumAngle = equilibriumAngle;

	VectorTools::copy(this->atoms, nodes);
}

HarmonicDihedralConstraintTerm::~HarmonicDihedralConstraintTerm(){}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function returns the angle formed by the four term atoms
/// using calculateDihedralAngleWithArcTanFunction
///
/// \param coords [In] vector containing the space coordinates of
/// the four term atoms
///
/// \return dihedral angle
///
/// \author arincon
/// \date
///////////////////////////////////////////////////////////////
double HarmonicDihedralConstraintTerm::getCurrentAngle(const double * const coords){
	Atom *a1 = atoms[0];
	Atom *a2 = atoms[1];
	Atom *a3 = atoms[2];
	Atom *a4 = atoms[3];

	return HarmonicDihedralConstraintFunctions::calculateDihedralAngleWithArcTanFunction(
			coords[a1->ix], coords[a1->iy], coords[a1->iz],
			coords[a2->ix], coords[a2->iy], coords[a2->iz],
			coords[a3->ix], coords[a3->iy], coords[a3->iz],
			coords[a4->ix], coords[a4->iy], coords[a4->iz]);
}
///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the energy of the constraint term
/// formed by the four term atoms using calculateEnergyWithArcTanFunction
///
/// \param coords [In] vector containing the space coordinates of
/// the four term atoms
///
/// \return energy of the constraint term
///
/// \author arincon
/// \date 20/08/2012
///////////////////////////////////////////////////////////////
double HarmonicDihedralConstraintTerm::getEnergy(const double * const coords)
{
	// Get atoms from list
	Atom *a1 = atoms[0];
	Atom *a2 = atoms[1];
	Atom *a3 = atoms[2];
	Atom *a4 = atoms[3];

	// Total energy
	return HarmonicDihedralConstraintFunctions::calculateEnergyWithArcTanFunction(springConstant,
																coords[a1->ix], coords[a1->iy], coords[a1->iz],
																coords[a2->ix], coords[a2->iy], coords[a2->iz],
																coords[a3->ix], coords[a3->iy], coords[a3->iz],
																coords[a4->ix], coords[a4->iy], coords[a4->iz],
																equilibriumAngle);
}


///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the gradient using calculateGradient function
/// and updates the grad vector using updateGradient
///
/// \param coords [In] vector containing the space coordinates of
/// the four term atoms
/// \param grad [Out] Gradient vector
///
/// \return energy of the constraint term
///
/// \author arincon
/// \date 20/08/2012
///////////////////////////////////////////////////////////////
double HarmonicDihedralConstraintTerm::getEnergyGradient(const double * const coords, double * grad)
{
	Atom *a1 = atoms[0];
	Atom *a2 = atoms[1];
	Atom *a3 = atoms[2];
	Atom *a4 = atoms[3];

	vector<double> gradient = HarmonicDihedralConstraintFunctions::calculateGradient(springConstant, coords[a1->ix],
																					 coords[a1->iy], coords[a1->iz],
																					 coords[a2->ix], coords[a2->iy],
																					 coords[a2->iz], coords[a3->ix],
																					 coords[a3->iy], coords[a3->iz],
																					 coords[a4->ix], coords[a4->iy],
																					 coords[a4->iz], equilibriumAngle);

	GradientHessianDihedralUpdater::updateGradient(grad, a1, a2, a3, a4, 1, &(gradient[0]));

	// Total energy
	return getEnergy(coords);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates saves into a full matrix the value
/// of the hessian matrix stored in triangularMatrix
///
/// \param triangularMatrix [In] hessian triangular matrix
/// \param fullMatrix [Out] hessian full matrix (12x12 size)
///
/// \author arincon
/// \author vgil
/// \date 20/08/2012
///////////////////////////////////////////////////////////////
void HarmonicDihedralConstraintTerm::generateFullHessianMatrixFromTriangularMatrix(vector<double> &triangularMatrix,
																					double fullMatrix[][12]) {
	// The upper part of submatrix is initialized to zero
	for(int i=0; i<12; ++i) {
			for(int j=0; j<12; ++j) {
			fullMatrix[i][j] = 0.0;
		}
	}

	int vector_index = 0;
	for(int i=0; i<12; ++i) {
		for(int j=i; j<12; ++j) {
			fullMatrix[i][j] = triangularMatrix[vector_index];
			vector_index++;
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the constrained energy, gradient and hessian of the constraint
///
/// \param coords [In] vector containing the space coordinates of
/// the four term atoms
/// \param grad [Out] gradient vector
/// \param hessian [Out] hessian matrix
///
/// \return energy of the constraint
///
/// \author arincon
/// \date 21/08/2012
///////////////////////////////////////////////////////////////
double HarmonicDihedralConstraintTerm::getEnergyGradientHess(const double * const coords, double *grad, SparseMatrixCRS *hessian)
{
	HarmonicDihedralConstraintTerm::getHess(coords, hessian);

	double energy = getEnergyGradient(coords, grad);

	return energy;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the hessian matrix using calculateHessian function
/// and updates the value using updateHessian
///
/// \param coords [In] vector containing the space coordinates of
/// the four term atoms
/// \param hessian [Out] hessian matrix
///
/// \return energy of the constraint
///
/// \author arincon
/// \date 21/08/2012
///////////////////////////////////////////////////////////////
void HarmonicDihedralConstraintTerm::getHess(const double * const coords, SparseMatrixCRS *hessian){
	Atom *a1 = atoms[0];
	Atom *a2 = atoms[1];
	Atom *a3 = atoms[2];
	Atom *a4 = atoms[3];

	// Compute hessian
	vector<double> triangularHessian= HarmonicDihedralConstraintFunctions::calculateHessian(springConstant, coords[a1->ix],
																							coords[a1->iy], coords[a1->iz],
																							coords[a2->ix], coords[a2->iy],
																							coords[a2->iz], coords[a3->ix],
																							coords[a3->iy], coords[a3->iz],
																							coords[a4->ix], coords[a4->iy],
																							coords[a4->iz], equilibriumAngle);

	// Transform it to 2D matrix
	double submatrix[12][12];
	HarmonicDihedralConstraintTerm::generateFullHessianMatrixFromTriangularMatrix(triangularHessian, submatrix);

	// Update hessian
	GradientHessianDihedralUpdater::updateHessian(a1, a2, a3, a4, submatrix, hessian);
}

