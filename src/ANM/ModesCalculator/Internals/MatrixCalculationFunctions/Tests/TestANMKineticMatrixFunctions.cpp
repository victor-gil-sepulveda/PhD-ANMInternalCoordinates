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

    TEST_REGRESSION_FUNCTION(testCalculateI,
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/I_INMA.txt",
		INMA,
		DO_NOT_SKIP_OXT,
		1e-6);

    TEST_FUNCTION(testCalculateI,
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/I_INMA.txt",
    	INMA_NO_OXT,
    	SKIP_OXT,
    	1e-6);

    TEST_FUNCTION(testCalculateI,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/I.txt",
		INMA,
		DO_NOT_SKIP_OXT,
		1e-6);

    TEST_REGRESSION_FUNCTION(testCalculateI,
       	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
       	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/I_BRAUN.txt",
       	BRAUN,
    	SKIP_OXT,
    	1e-6);

    TEST_FUNCTION(testCalculateI,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/I.txt",
		INMA,
		DO_NOT_SKIP_OXT,
		1e-7);

    // Uses OXT
    TEST_FUNCTION(testCalculateI,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/I.txt",
		INMA,
		DO_NOT_SKIP_OXT,
		1e-7);

    // K  wo skipping OXT
    TEST_FUNCTION(testCalculateK,
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/K.txt",
		INMA,
		1e-6);

    // K skipping OXT
    TEST_FUNCTION(testCalculateK,
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
    	"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/K.txt",
    	INMA_NO_OXT,
    	1e-5);

    // K  wo skipping OXT
    TEST_FUNCTION(testCalculateK,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/K.txt",
		INMA,
		1e-6);

    // K matrices are very sensitive to small changes
    // K  wo skipping OXT
	TEST_FUNCTION(testCalculateK,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/K.txt",
		INMA,
		4e-3);

	TEST_FUNCTION(testCalculateK,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/K.txt",
		INMA,
		6e-4);

//    TEST_FUNCTION(testJacobi,
//    		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
//    		1e-12)

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

	// Check
	bool ok = true;
	for (unsigned int i = 0; i < expected_K.size(); ++i){

		vector<double> k_row;
		TriangularMatrixRow K_ith_row (*K,i);
		for(unsigned int j = 0; j < expected_K[i].size();++j){
			if(j>=i){
				k_row.push_back(K_ith_row(j));
			}
			else{
				k_row.push_back(0.0);
			}
		}

		ok &= Assertion::expectedVectorHasRelativeErrorNotBiggerThan(
							expected_K[i],
							k_row,
							tolerance);
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
	vector<vector<double> > J,Jt,JtM,JtMJ;
	ANMICKineticMatrixCalculator::Jacobi2(units, J);// Transposed Jacobi
	ANMICMath::transpose(J,Jt);

	// M matrix
	bool onlyHeavyAtoms = true;
	vector<Atom*> all_atoms;
	UnitTools::getAllAtomsFromUnits(units, all_atoms, onlyHeavyAtoms);
	vector<vector<double> > M(all_atoms.size()*3, vector<double>(all_atoms.size()*3,0));
	for (unsigned int i = 0; i < all_atoms.size(); ++i){
		unsigned int offset = i*3;
		M[offset][offset] = all_atoms[i]->getMass();
		M[offset+1][offset+1] = all_atoms[i]->getMass();
		M[offset+2][offset+2] = all_atoms[i]->getMass();
	}


	ANMICMath::multiplyMatrixByMatrix(Jt,M,JtM);
	ANMICMath::multiplyMatrixByMatrix(JtM,J,JtMJ);

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

	// Also, sum(MiJia) = 0
	for (unsigned int alpha = 0; alpha < Jt.size(); ++alpha){
		double total = 0;
		for  (unsigned int i = 0; i< Jt[alpha].size(); ++i){
			total += M[i][i]*Jt[alpha][i];
		}
		cout<<"** "<<total<<endl;
	}

	// And sum(mi rixJia) = 0
	for (unsigned int alpha = 0; alpha < Jt.size(); ++alpha){
		double total = 0;
		for  (unsigned int i = 0; i< Jt[alpha].size()/3; ++i){
			Point ri = all_atoms[i]->toPoint();
			Point Jia(Jt[alpha][i*3],Jt[alpha][i*3+1],Jt[alpha][i*3+2]);
			Point rixJia = ANMICMath::crossProduct(ri, Jia);
			total += M[i][i]* rixJia.getX();
			total += M[i][i]* rixJia.getY();
			total += M[i][i]* rixJia.getZ();
		}
		cout<<"xx "<<total<<endl;
	}
}
