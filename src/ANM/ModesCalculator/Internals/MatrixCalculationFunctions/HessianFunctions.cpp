/////////////////////////////////////////////////////////////////////////////
/// ANMICHessianCalculation.cpp
///
/// Implementation of ANMICHessianCalculation
///
/// \author arincon
/// \date 17/11/2012
/////////////////////////////////////////////////////////////////////////////

#include "HessianFunctions.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <iosfwd>
#include "../../../../Tools/Math/Point.h"
#include "../../../../Molecules/Atom.h"
#include "../CoarseGrainModel/Unit.h"
#include <iomanip>
#include "ANMICMath.h"
#include <cstring>
#include <algorithm>
#include "TriangularMatrices.h"
#include "../ElasticConstantCalculator.h"
#include <string>
#include "../../../../Molecules/AtomNames.h"
using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates Dij given two points.
//	Dij is independent of global translations
///
/// \param fij [In] force constant of the spring between
//					the two atoms in the points
/// \param r_i [In] point i
/// \param r_j [In] point j
///
/// \param res [Out] Dij matrix
///
/// \author arincon
/// \author vgil
/// \date 20/11/2013
///////////////////////////////////////////////////////////////
void ANMICHessianCalculator::calculateDij(double k, Point& r_i, Point& r_j, double *res){
	//---------------------
	//!!!!!!!!!!!!!!!!
	//double fij = 1;
	//double term = fij / r_j.squaredDistance(r_i);
	//term = 1; // !!!
	//-------------------

	Point Rij = Point::subtract(r_i, r_j);

	Point multRiRj = ANMICMath::crossProduct(r_i, r_j);

	double A[6] = {multRiRj.getX(), multRiRj.getY(), multRiRj.getZ(), Rij.getX(), Rij.getY(), Rij.getZ()};
	double B[6] = {multRiRj.getX(), multRiRj.getY(), multRiRj.getZ(), Rij.getX(), Rij.getY(), Rij.getZ()};

	ANMICMath::multiplyColumnByRow(A, B, res, 6);
	double d = r_j.squaredDistance(r_i);
	ANMICMath::multiplyMatrixByScalar(k/d, res, 6);

}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates Tab (the sum of the Dij matrices)
///	given two units
///
/// \param unit_a [In] a unit
/// \param unit_b [In] b unit
/// \param sq_cutoff [In] squared cutoff distance
///
/// \param result [Out] Tab
///
/// \author arincon
/// \author vgil
/// \date 01/12/2013
///////////////////////////////////////////////////////////////
void ANMICHessianCalculator::calculateTab(Unit* unit_a, Unit* unit_b,
												double sq_cutoff,
												ElasticConstantCalculator* force_constant_calculator,
												double *Tab,
												bool skip_OXT) {
	double Dij[6*6];

	for(unsigned int i = 0; i < unit_a->atoms.size(); i++) {
		Atom* atom_i = unit_a->atoms[i];
		Point r_i = atom_i->toPoint();
		for(unsigned int j = 0; j < unit_b->atoms.size(); j++) {
			Atom* atom_j = unit_b->atoms[j];
			Point r_j = atom_j->toPoint();
//			double distance = r_i.distance(r_j);
			if(r_i.squaredDistance(r_j) < sq_cutoff
					&& !(skip_OXT && (atom_i->name == AtomNames::OXT
					|| atom_j->name == AtomNames::OXT))){

//				cout<<"DBG: "<<atom_i->name<<" ("<<atom_i->serial<<")"<<" "
//						<<atom_j->name<<" ("<<atom_j->serial<<")"<<" i "
//						<< i <<" j " << j <<" C "
//						<<force_constant_calculator->calculateK(r_i, r_j)<<" d "<<distance<<endl;

				ANMICHessianCalculator::calculateDij(force_constant_calculator->calculateK(r_i, r_j),
														r_i, r_j, Dij);
				ANMICMath::addMatrices(Tab, Dij, Tab, 6);

			}
			else{
//				cout<<"DBG: Out of scope "<<atom_i->name<<" "<<atom_j->name<<endl;
			}
		}
	}
//	ANMICMath::printMatrix(Tab,6);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the U matrix containing all the
///	interaction information between each unit
///
/// \param units [In] Units of the molecule
///
/// \return double matrix
///			The interaction between each pair of units is
///			represented by a 6x6 matrix, so the U has a
///			size of (6*M)x(6*M)
///
/// \author arincon
/// \author vgil
/// \date 11/12/2013
///////////////////////////////////////////////////////////////
// Si no caben, U se tiene que transformar en una matriz triangular
std::vector< std::vector<double> > ANMICHessianCalculator::calculateU(double sq_cutoff,
															ElasticConstantCalculator* force_constant_calculator,
															std::vector<Unit*>& units,
															bool skip_OXT) {
	int M = units.size();
	int colLength = M * 6;

	// initialize U matrix
	std::vector< std::vector<double> > UMatrix(colLength, std::vector<double>(colLength));

	double U_tmp_buffer[6*6];
	double Tab[6*6];
	// units a
	for(int a = 0; a < M; ++a) {
		for(int b = M-1; b > a; --b) {
			// Tab = {0, 0, .... , 0}
			memset(Tab, 0, 6*6*sizeof(double));

			// First row and first column are treated differently than the rest of position
			// If b is not in the last column (A)
			if(b != M-1) {
				ANMICMath::getABMatrix(UMatrix, U_tmp_buffer, a, b+1);
				ANMICMath::addMatrixToU(UMatrix, U_tmp_buffer, a, b);
			}

			// If a is not in the first row (B)
			if(a != 0) {
				//U[a][b] += U[a-1][b];
				ANMICMath::getABMatrix(UMatrix, U_tmp_buffer, a-1, b);
				ANMICMath::addMatrixToU(UMatrix, U_tmp_buffer, a, b);
			}

			// All the other cases
			if(b+1 < M && a-1 >= 0) {
				//U[a][b] -= U[a-1][b+1];
				ANMICMath::getABMatrix(UMatrix, U_tmp_buffer, a-1, b+1);
				ANMICMath::subsMatrixFromU(UMatrix, U_tmp_buffer, a, b);
			}

			//U[a][b] += Tab;
//			cout<<"DBG: calc "<<a <<" "<<b <<endl;
			ANMICHessianCalculator::calculateTab(units[a], units[b],
													sq_cutoff,
													force_constant_calculator,
													Tab,
													skip_OXT);

			ANMICMath::addMatrixToU(UMatrix, Tab, a, b);

			ANMICMath::getABMatrix(UMatrix, U_tmp_buffer, a, b);
		}
	}

	return UMatrix;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates Rab (Ua,b+1), a matrix
///	containing the interaction information between all
///	the atoms of the units "a" and "b"
///
/// \param a [In] a unit
/// \param b [In] b unit
/// \param U [In] U matrix
///
/// \param Rab [Out] 6x6 matrix
///
/// \author arincon
/// \date 11/12/2013
///////////////////////////////////////////////////////////////
void ANMICHessianCalculator::calculateRab(int a, int b, std::vector< std::vector<double> > &U, double *Rab) {
	ANMICMath::getABMatrix(U, Rab, a, b+1);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This functions calculates an element of the hessian matrix
///
/// \param e_a [In] e_a point
/// \param e_b [In] e_b point
/// \param r_a [In] r_a point
/// \param r_b [In] r_b point
/// \param a_dihedral [In] a dihedral angle
/// \param b_dihedral [In] b dihedral angle
/// \param U [In] U matrix
///
/// \return double Hab
///
/// \author arincon
/// \date 11/12/2013
///////////////////////////////////////////////////////////////
double ANMICHessianCalculator::calculateHab(Point* e_a, Point* e_b, Point* r_a, Point* r_b,
		int a_dihedral, int b_dihedral, std::vector< std::vector<double> > &U) {

	// Calculate cross products
	Point crossEaRa = ANMICMath::crossProduct(*e_a, *r_a);
	Point crossEbRb = ANMICMath::crossProduct(*e_b, *r_b);

	double A[6] = {e_a->getX(), e_a->getY(), e_a->getZ(), crossEaRa.getX(), crossEaRa.getY(), crossEaRa.getZ()};
	double B[6] = {e_b->getX(), e_b->getY(), e_b->getZ(), crossEbRb.getX(), crossEbRb.getY(), crossEbRb.getZ()};

	double Rab[6*6];
	double result[6];

	ANMICHessianCalculator::calculateRab(a_dihedral, b_dihedral, U, Rab);

	ANMICMath::multitplyRowByMatrix(A, Rab, result, 6);
	double hab = ANMICMath::multiplyRowByColumn(result, B, 6);

	return hab;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the hessian matrix
///
/// \param units [In] vector containing all the molecule units
/// \param U [In] U matrix
///
/// \return double matrix
///			Hessian matrix of M*M values
///
/// \author arincon
/// \author vgil
/// \date 11/12/2013
///////////////////////////////////////////////////////////////
TriangularMatrix* ANMICHessianCalculator::calculateH(std::vector<Unit*>& units, std::vector< std::vector<double> > &U){
	unsigned int M = units.size()-1;
	TriangularMatrix* H = new TriangularMatrix(M, M);

	// Iterates over the dihedrals of M units (psi0, phi1, psi1... phiM-1, psiM-1)
	for(unsigned int alpha_dihedral = 0; alpha_dihedral < M; ++alpha_dihedral){
		Unit* unit_a = units[alpha_dihedral];
		for(unsigned int beta_dihedral = alpha_dihedral; beta_dihedral < M; ++beta_dihedral){
			Unit* unit_b = units[beta_dihedral];
			(*H)(alpha_dihedral,beta_dihedral) = ANMICHessianCalculator::calculateHab(unit_a->e_right, unit_b->e_right,
					unit_a->r_right, unit_b->r_right,
					alpha_dihedral, beta_dihedral, U);
		}
	}

	return H;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It tries to lower the tip effect adding extra dihedral constraints to the hessian.
///
/// \param H [In] Hessian matrix
///
/// \author vgil
/// \date 11/12/2013
///////////////////////////////////////////////////////////////

void ANMICHessianCalculator::modifyHessianWithExtraTorsion(TriangularMatrix* H){
	vector<double> diagonal;
	vector<double> unsigned_diagonal;
	for (unsigned int i = 0; i < H->size1(); ++i) {
		//ANM TODO:
		diagonal.push_back((*H)(i,i));
		unsigned_diagonal.push_back((*H)(i,i));
	}

	double min = *min_element(unsigned_diagonal.begin(), unsigned_diagonal.end());
	unsigned int index_of_min = 0;
	for(unsigned int i =0; i < unsigned_diagonal.size(); ++i){
		if(min == unsigned_diagonal[i])
			index_of_min = i;
	}

	min = diagonal[index_of_min];
	double omega = 3*(min);
//	cout<<"DBG: Omega "<<omega<<endl;
	for (unsigned int i = 0; i < H->size1(); ++i) {
		(*H)(i,i) += omega;
	}
}
