/////////////////////////////////////////////////////////////////////////////
/// TestANMKineticMatrixFunctions.cpp
///
/// Implementation of TestANMKineticMatrixFunctions class
///
/// \author vgil
/// \author arincon
/// \date 12/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "TestANMKineticMatrixFunctions.h"
#include "TestANMICTools.h"

#include "../../../../../System/System.h"
#include "../../CoarseGrainModel/Unit.h"
#include "../../../../../Molecules/Complex.h"

#include <iostream>
#include <vector>
#include "../KineticMatrixFunctions.h"
#include "../../../../../Tools/Utils.h"
#include "../ANMICMath.h"
#include "../../../../../Tools/TestTools.h"
#include "../../../../../Tools/Assertion.h"
#include "../TriangularMatrices.h"
#include "../../CoarseGrainModel/UnitTools.h"
using namespace std;

TestANMKineticMatrixFunctions::TestANMKineticMatrixFunctions(string name)
{
    test_name = name;
}

TestANMKineticMatrixFunctions::~TestANMKineticMatrixFunctions()
{}

void TestANMKineticMatrixFunctions::init()
{}

void TestANMKineticMatrixFunctions::run()
{
    Test::run();

    bool SKIP_OXT = true;
    bool DO_NOT_SKIP_OXT = false;

//    TEST_REGRESSION_FUNCTION(testCalculateI,
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/I_INMA.txt",
//		INMA,
//		DO_NOT_SKIP_OXT,
//		1e-6);
//
//    TEST_FUNCTION(testCalculateI,
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/I_INMA.txt",
//    	INMA_NO_OXT,
//    	SKIP_OXT,
//    	1e-6);
//
//    TEST_FUNCTION(testCalculateI,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/I.txt",
//		INMA,
//		DO_NOT_SKIP_OXT,
//		1e-6);
//
//    TEST_REGRESSION_FUNCTION(testCalculateI,
//       	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
//       	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/I_BRAUN.txt",
//       	BRAUN,
//    	SKIP_OXT,
//    	1e-6);
//
//    TEST_FUNCTION(testCalculateI,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/I.txt",
//		INMA,
//		DO_NOT_SKIP_OXT,
//		1e-7);
//
//    // Uses OXT
//    TEST_FUNCTION(testCalculateI,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/I.txt",
//		INMA,
//		DO_NOT_SKIP_OXT,
//		1e-7);
//
//    // K  wo skipping OXT
//    TEST_FUNCTION(testCalculateK,
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/K.txt",
//		INMA,
//		1e-4);
//
//    // K skipping OXT
//    TEST_FUNCTION(testCalculateK,
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
//    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/K.txt",
//    	INMA_NO_OXT,
//		1e-3);
//
//    // K  wo skipping OXT
//    TEST_FUNCTION(testCalculateK,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/K.txt",
//		INMA,
//		1e-4);
//
//    //This test fails because of unknown reasons
//    // Relative errors > 0.1: 1.3 (+/- scenario),  0.169, 0.1046, 0.178, 0.111, 0.907
//    // K  wo skipping OXT
//	TEST_FUNCTION(testCalculateK,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/K.txt",
//		INMA,
//		0.1);

//	TEST_FUNCTION(testCalculateK,
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
//		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/K.txt",
//		INMA,
//		0.1);
    TEST_FUNCTION(testJacobi,
    		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
    		1e-12)

    finish();
}

void TestANMKineticMatrixFunctions::finish()
{
    Test::finish();
}

bool TestANMKineticMatrixFunctions::testCalculateI(
		const char* complex_pdb,
		const char* goldenI,
		ICalcType calc_type,
		bool skip_OXT,
		double tolerance){

	cout<<"Calculate I for "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;
	TestANMICTools::createUnitsFromFile(complex_pdb, units, complex, skip_OXT);

	vector<vector<double> > expectedI;
	TestTools::load_vector_of_vectors(expectedI, goldenI);

	double I[3][3];

	ANMICKineticMatrixCalculator::calculateI(I,
			units,
			pair<int,int>(0,units.size()-1),
			calc_type);

	bool ok = true;
	for(unsigned int  i = 0; i < expectedI.size(); ++i){
		for(unsigned int j = 0 ; j < expectedI[i].size(); ++j){
			// Low precission as it is read from a file with decs.
			ok &= Assertion::expectedEqualsCalculatedWithTolerance(expectedI[j][i], I[i][j], tolerance);
		}
	}

	Utils::clearVector<Unit>(units);
	delete complex;

	return ok;
}

bool TestANMKineticMatrixFunctions::testCalculateK(const char* complex_pdb,
		const char* goldenK,
		ICalcType tensor_calc_type,
		double tolerance){

	cout<<"Calculate K for "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;

	bool skip_OXT = tensor_calc_type == INMA_NO_OXT? true: false;
	TestANMICTools::createUnitsFromFile(complex_pdb, units, complex, skip_OXT);

	// Load expected K
	vector<vector<double> > expected_K;
	TestTools::load_vector_of_vectors(expected_K, goldenK);

	// Calculate K
	TriangularMatrix* K = ANMICKineticMatrixCalculator::calculateK(units, tensor_calc_type);

//	cout<<"DBG: K size ("<<K->size1()<<", "<<K->size2()<<")"<<endl;
//	for (unsigned int i = 0; i < K->size1(); ++i){
//		TriangularMatrixRow K_r (*K,i);
//		for (unsigned int j = i; j < K->size2(); ++j){
//			cout<<K_r(j)<<", ";
//		}
//		cout<<endl;
//	}

	bool ok = true;
	for (unsigned int i = 0; i < expected_K.size(); ++i){
		TriangularMatrixRow K_ith_row (*K,i);
		for (unsigned int j = i; j < expected_K[i].size(); ++j){
			bool one_check = Assertion::expectedEqualsCalculatedWithTolerance(expected_K[i][j], K_ith_row(j), tolerance);
			ok &= one_check;
			if(one_check == false ){
				cout<<"Failed at "<<i<<" "<<j<<endl;
			}
		}
	}

	Utils::clearVector<Unit>(units);
	delete complex;
	delete K;
	return ok;
}


bool TestANMKineticMatrixFunctions::testJacobi(const char* prot_path, double tolerance){

	System sys;
	vector<Unit*> units;
	Complex* complex;
	TestANMICTools::createUnitsFromFile(prot_path, units, complex, false);

	// Load expected K
	vector<vector<double> > expected_K;

	// Calculate K
	TriangularMatrix* K = ANMICKineticMatrixCalculator::calculateK(units, INMA);

	// If the Jacobi is ok, K => T = JtMJ
	vector<vector<double> > J;
	ANMICKineticMatrixCalculator::Jacobi(units, J);// Transposed Jacobi

	// Get M matrix
	bool onlyHeavyAtoms = true;
	vector<Atom*> all_atoms;
	UnitTools::getAllAtomsFromUnits(units, all_atoms, onlyHeavyAtoms);
	vector<double> M;
	for (unsigned int i = 0; i < all_atoms.size(); ++i){
		M.push_back(all_atoms[i]->getMass());
		M.push_back(all_atoms[i]->getMass());
		M.push_back(all_atoms[i]->getMass());
	}

	// M is diagonal matrix of the masses
	vector<vector<double> > JtM;
	for (unsigned int i =0; i< J.size(); ++i){ // d
		vector<double> row;
		for (unsigned int j =0; j< M.size(); ++j){ // n
			row.push_back(J[i][j]*M[j]);
		}
		JtM.push_back(row);
	}

	// Now get JtMJ (regular multiplication, remember that J is already trasposed,
	// so we have to retranspose during the mult.)
	vector<vector<double> > JtMJ;
	for (unsigned int i =0; i< JtM.size(); ++i){ //n
		vector<double> row;
		for (unsigned int j =0; j< J.size(); ++j){//d
			// perform dot between row and column of J (row too)
			double dot = 0;
			for (unsigned int k=0; k< J[0].size(); ++k){
				dot += JtM[i][k]*J[j][k];
			}
			row.push_back(dot);
		}
		JtMJ.push_back(row);
	}

		cout<<"DBG: K size ("<<K->size1()<<", "<<K->size2()<<")"<<endl;
		for (unsigned int i = 0; i < K->size1(); ++i){
			TriangularMatrixRow K_r (*K,i);
			for (unsigned int j = i; j < K->size2(); ++j){
				cout<<K_r(j)<<", ";
			}
			cout<<endl;
		}

		cout <<"--------"<<endl;
	for (unsigned int i = 0; i< JtMJ.size(); ++i){
		for (unsigned int j = 0; j< JtMJ[0].size(); ++j){
			cout<<JtMJ[i][j]<<", ";
		}
		cout<<endl;
	}

}

