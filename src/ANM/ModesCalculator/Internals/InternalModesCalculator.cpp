/////////////////////////////////////////////////////////////////////////////
/// InternalModesCalculator.cpp
///
/// Implementation of InternalModesCalculator class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////

#include "InternalModesCalculator.h"
#include "MatrixCalculationFunctions/HessianFunctions.h"
#include "MatrixCalculationFunctions/KineticMatrixFunctions.h"
#include "../../Parameters/AnmParameters.h"
#include "../../ModesCalculator/AnmEigen.h"
#include "../../AnmNodeList.h"
#include "../../AnmUnitNodeList.h"
#include "CoarseGrainModel/Unit.h"
#include <vector>
#include "../../../Tools/Utils.h"
#include "../../../Tools/Math/Point.h"
#include "MatrixCalculationFunctions/ANMICMath.h"
#include <string>
#include "../../../Tools/TestTools.h"
#include "ElasticConstantCalculator.h"
#include <iostream>
#include "../../Tools/Inout/ModesWriter.h"
#include "InverseExponentialElasticConstant.h"
#include "CoarseGrainModel/UnitTools.h"
#include <map>
using namespace std;

//-------------------------------------
// C definition of the needed functions from LAPACK. This definitions must be moved to the Math package upon
// final merge.
//-------------------------------------
extern "C" {
	void dspgvx_(int* itype, char* jobz, char* range, char* uplo, int* N, double* A, double* B, double* vl, double* vu,
			int* il, int* iu, double* abstol, int* M, double*W, double* Z, int* ldz, double* work,	int* iwork, int* ifail,
			int* info);
	double dlamch_(char* cmach);
}
//--------------------------------


InternalModesCalculator::InternalModesCalculator() : ModesCalculator() {
}

InternalModesCalculator::~InternalModesCalculator(){

}


///////////////////////////////////////////////////////////////
/// \remarks
/// Calculation of the eigenvalues and vectors of a given structure defined by its coarse grained model.
///
/// \param anmParameters [In ] Parameters of this ANM calculation
/// \param node_list [In ] Stores the coarse grain representation of the structure in Unit objects.
/// \param eigen [In / Out] The object that will hold the eigenvectors and values once calculated.
///
/// \author vgil
/// \date 07/01/2015
///////////////////////////////////////////////////////////////
void InternalModesCalculator::calculateEigenValuesAndVectors(AnmParameters * anmParameters,
																const AnmNodeList & node_list,
																AnmEigen * eigen)
{
	// const std::vector<Atom*> & node_list viene de  NodeListGenerator::selectNodes
	// Build coarse grain model
	const AnmUnitNodeList &unitNodeList = dynamic_cast<const AnmUnitNodeList &>(node_list);
	std::vector<Unit*> units = unitNodeList.getNodeList(); // TODO: RELLENAR modelo
	cout<<"DBG: Calculating Eig. for "<< units.size() <<" units."<<endl;

	// Calculate K and H
	double cutoff = anmParameters->getCutoff();
	double k = anmParameters->getConstantForHessian();
	cout<<"DBG: InternalModesCalculator::calculateEigenValuesAndVectors K "<<k<<" cutoff "<<cutoff<<endl;
	//ElasticConstantCalculator ecc(k);
	InverseExponentialElasticConstant ecc(k, 3.8, 6);
	std::vector< std::vector<double> > U = ANMICHessianCalculator::calculateU(cutoff*cutoff,
			&ecc,
			units,
			false); // Skip OXT

	TriangularMatrix* H = ANMICHessianCalculator::calculateH(units,  U);
	//ANM TODO: This must be an option !!
	//ANMICHessianCalculator::modifyHessianWithExtraTorsion(H);

	TriangularMatrix* K = ANMICKineticMatrixCalculator::calculateK(units, INMA);

	// And then the modes
	calculate_modes(anmParameters, H, K, eigen);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Does the actual eigencalculation using K and H matrices.
///
/// \param anmParameters [In ] Parameters of this ANM calculation
/// \param H [In / Out] H matrix (hessian matrix). May be modified by LAPACK.
/// \param K [In / Out] K matrix (kinetic matrix). May be modified by LAPACK.
/// \param eigen [In / Out] The object that will hold the calculated eigenvectors and values.
///
/// \author vgil
/// \date 07/01/2015
///////////////////////////////////////////////////////////////
void InternalModesCalculator::calculate_modes(AnmParameters * anmParameters,
													TriangularMatrix* H,
													TriangularMatrix* K,
													AnmEigen* eigen){

	TestTools::save_vector( H->data().begin(), "H.txt", 0, H->data().size()-1, 4);
	TestTools::save_vector( K->data().begin(), "K.txt", 0, K->data().size()-1, 4);

	int number_of_eigen = anmParameters->getNumberOfModes();
	int itype = 1;
	char jobz = 'V';
	char range = 'I';
	char uplo = 'U';
	int N = K->size1();
	int il = 1;
	int iu = number_of_eigen;
	char cmach = 'S';
	double abstol = 2*dlamch_(&cmach);
	//override
	abstol = 0.0001;
	int M;
	double* W = new double[number_of_eigen];
	double* Z = new double[N*number_of_eigen];
	double * WORK = new double[N*8];
	int * IWORK = new int[N*5];
	int * IFAIL = new int[N];
	int info;
	dspgvx_(&itype,	// ITYPE A*x = (lambda)*B*x
			&jobz,	// JOBZ Compute eigenvalues and eigenvectors.
			&range,	// RANGE
			&uplo,
			&N,
			H->data().begin(),
			K->data().begin(),
			NULL,	// vl vu
			NULL,
			&il,	// il iu
			&iu,
			&abstol,// ABSTOL
			&M,		// M The total number of eigenvalues found
			W,		// W (vector) eigenvalues in ascending order
			Z,		// Z (Matrix) eigenvectors
			&N,		// LDZ leading dimension de Z
			WORK,
			IWORK,
			IFAIL,	// IFAIL  contains the indices of the eigenvectors that failed to converge
			&info);

	delete [] WORK;
	delete [] IWORK;
	delete [] IFAIL;

	eigen->initialize(W, Z, number_of_eigen, H->size1(), false);

	delete [] W;
	delete [] Z;

}

///////////////////////////////////////////////////////////////
/// \remarks
/// Converts the modes expressed in internal coordinates to modes expressed in
/// cartesian coordinates (all atom). See code for details on calculations.
///
/// \param units [In ] Coarse grain representation of the structure in Unit objects.
/// \param in [In / Out] AnmEigen object containing modes in internal coordinates.
///
/// \return The AnmEigen object holding the converted modes.
///
/// \author vgil
/// \date 07/01/2015
///////////////////////////////////////////////////////////////
AnmEigen* InternalModesCalculator::internalToCartesian(vector<Unit*>& units, AnmEigen* in, bool onlyHeavy) {

	AnmEigen* out = new AnmEigen;
	unsigned int number_of_modes = in->getNumberOfModes();
	unsigned int number_of_dihedrals = units.size() - 1;

	vector<double> eigenvalues_cc;
	vector< vector<double> > eigenvectors_cc;

	vector<Atom*> unit_atoms;
	UnitTools::getAllAtomsFromUnits(units, unit_atoms, onlyHeavy);
	cout<<"DBG: Conversion from IC to CC. N.Units: "<<units.size()<<" N.Atoms "<<(onlyHeavy? "(Heavy): ": "(All): ")<<unit_atoms.size()<<endl;

	// Create a dictionary to index them using the pointer
	map<Atom*,int> atom2index;
	for (unsigned int i =0; i< unit_atoms.size();++i){
		atom2index[unit_atoms[i]] = i;
	}

	// Precalculate terms of ecs 4.1.12 and .13 in J.R. Lopez Blanco Thesis (dr/dq)
	vector<Point> term1, term2, term1b, term2b;
	double I[3][3], I_inv[3][3];
	ANMICKineticMatrixCalculator::calculateI(I, units, pair<int,int>(0, units.size()-1), INMA); // calculate I for all units i in [0,number_of_units-1]
	ANMICMath::invertMatrix(I,I_inv);
	for(unsigned int dihedral = 0; dihedral < number_of_dihedrals; ++dihedral){
		ANMICKineticMatrixCalculator::dri_dq( units, dihedral, I_inv, term1, term2);
		ANMICKineticMatrixCalculator::dri_dq_2( units, dihedral, I_inv, term1b, term2b);
	}

	// Calculate per-atom contributions of each dihedral turn
	for (unsigned int k = 0; k < number_of_modes; k++){
		// Add eigenvalue (no need to convert this one)
		eigenvalues_cc.push_back( in->getEigenValueOfMode(k));

		//Calc eigenvector for mode
		vector<double> eigenvector_ic = in->getEigenVectorOfMode(k); // Mode k , eigenvector in internal coords
		Point p_eigenvector_cc[unit_atoms.size()];
		vector<Atom*> all_atoms;
		UnitTools::getAllAtomsFromUnitRange(units, all_atoms, 0, units.size()-1, onlyHeavy);

		for(unsigned int dihedral = 0; dihedral < number_of_dihedrals; ++dihedral){
			vector<Atom*> left_atoms, right_atoms;
			UnitTools::getAllAtomsFromUnitRange(units, left_atoms, 0, dihedral, onlyHeavy);
			UnitTools::getAllAtomsFromUnitRange(units, right_atoms, dihedral+1, units.size()-1, onlyHeavy);

			// Contributions of dihedral i to left atoms
			// term1 - term2 x r1
			for(unsigned int i = 0; i < left_atoms.size(); ++i){
				Point r_i =  left_atoms[i]->toPoint();
				Point dr_dqa = Point::subtract(term1[dihedral],
						ANMICMath::crossProduct(term2[dihedral], r_i));
				dr_dqa.multiplyByScalar(eigenvector_ic[dihedral]);
				p_eigenvector_cc[atom2index[left_atoms[i]]].add(dr_dqa);
			}

			// Contributions of dihedral i to right atoms
			// -term1 + term2 x r1
			for(unsigned int i = 0; i < right_atoms.size(); ++i){
				Point r_i =  right_atoms[i]->toPoint();
				Point dr_dqb = Point::subtract(ANMICMath::crossProduct(term2b[dihedral], r_i),
									term1b[dihedral]);
				dr_dqb.multiplyByScalar(eigenvector_ic[dihedral]);
				p_eigenvector_cc[atom2index[right_atoms[i]]].add(dr_dqb);

			}
		}

		// Convert from point cc to flat vector
		vector<double> mode_cc;
		for (unsigned int i = 0; i < unit_atoms.size(); ++i){
			mode_cc.push_back(p_eigenvector_cc[i].getX());
			mode_cc.push_back(p_eigenvector_cc[i].getY());
			mode_cc.push_back(p_eigenvector_cc[i].getZ());
		}
		eigenvectors_cc.push_back(mode_cc);
	}

	out->initialize(eigenvalues_cc, eigenvectors_cc, false);

	return out;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Converts a collection of modes expressed in cartesian coordinates in torsional rotations (can be
/// seen as their conversion to IC.
///
/// \param units [In] Coarse grain representation of the structure in Unit objects.
/// \param in [In / Out] AnmEigen object containing modes (or forces or whatever) in cartesian coordinates.
///
/// \return The AnmEigen object holding the converted modes.
///
/// \author vgil
/// \date 01/03/2015
///////////////////////////////////////////////////////////////
AnmEigen* InternalModesCalculator::cartesianToInternal(vector<Unit*>& units, AnmEigen* in){
	AnmEigen* out = new AnmEigen;
	unsigned int number_of_modes = in->getNumberOfModes();
	unsigned int number_of_coordinates = in->getEigenVectorsDimension();

	// Calculate Jacobi
	vector<vector<double> > J;
	ANMICKineticMatrixCalculator::Jacobi(units, J, M);

	// Calculate K
	TriangularMatrix* K = ANMICKineticMatrixCalculator::calculateK(units, INMA);

	// Invert K
	TriangularMatrix* K_inv;

	// Get all masses
	bool onlyHeavyAtoms = true;
	vector<Atom*> all_atoms;
	UnitTools::getAllAtomsFromUnits(units, all_atoms, onlyHeavyAtoms);
	vector<double> M;
	for (unsigned int i = 0; i < all_atoms.size(); ++i){
		M.push_back(all_atoms[i]->getMass());
		M.push_back(all_atoms[i]->getMass());
		M.push_back(all_atoms[i]->getMass());
	}

	// We do this for every of the modes we have in cartesian coordinates
	for (unsigned int i = 0;i < number_of_modes; ++i){
		vector<double> torsion_rotations;
		vector<double>& original_mode = in->vectors[i];
		// Multiply M by r


	}



	return out;
}

std::string InternalModesCalculator::generateReport() const {
	return "Modes calculated using internal coordinates";
}
