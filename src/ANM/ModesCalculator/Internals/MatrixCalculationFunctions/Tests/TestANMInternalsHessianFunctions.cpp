/////////////////////////////////////////////////////////////////////////////
/// TestANMICHessianCalculation.cpp
///
/// Implementation of TestANMICHessianCalculation class
///
/// \author arincon
/// \date 18/11/2012
/////////////////////////////////////////////////////////////////////////////

#include "TestANMInternalsHessianFunctions.h"
#include "../HessianFunctions.h"
#include <cmath>
#include <vector>
#include "../../../../../Tools/Math/Point.h"
#include "../../../../../Tools/Assertion.h"
#include "../../../../../Tools/stringTools.h"
#include "../../../../../Molecules/ComplexBuilder.h"
#include "../../../../../Molecules/Complex.h"
#include "../../../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../../../Molecules/Selection/Selector.h"
#include "../../CoarseGrainModel/UnitsBuilder.h"
#include "../../CoarseGrainModel/Unit.h"
#include "../../../../../System/System.h"
#include "../../../../../Tools/Math/MathTools.h"
#include "../../../../../Tools/TestTools.h"
#include "../../../../../Tools/Utils.h"
#include <iomanip>
#include "../ANMICMath.h"
#include "../../../../../Molecules/ElementParametersLoader.h"
#include "TestANMICTools.h"
#include "ElasticCalculatorMock.h"
#include "../../InverseExponentialElasticConstant.h"

using namespace std;
using StringTools::toDouble;
using StringTools::split;

TestANMInternalHessianFunctions::TestANMInternalHessianFunctions(string name) {
	test_name = name;
}

TestANMInternalHessianFunctions::~TestANMInternalHessianFunctions() {
}

void TestANMInternalHessianFunctions::init() {
}

void TestANMInternalHessianFunctions::run() {
	Test::run();

	bool SKIP_OXT = true;

	TEST_FUNCTION(testCalculateD,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/D.txt");

	TEST_FUNCTION(testCalculateT,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/T_iexp.txt",
		not SKIP_OXT,
		1e-5);

	TEST_FUNCTION(testCalculateT,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/T_iexp.txt",
		SKIP_OXT,
		1e-5);

	TEST_FUNCTION(testCalculateT,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/T_iexp.txt",
		not SKIP_OXT,
		1e-5);

	TEST_FUNCTION(testCalculateU,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/U.txt",
		not SKIP_OXT,
		1e-4);

	// These tests use an error check which is relative to the maximum absolute magnitude
	// of the array/ matrix. This is done because H matrix (as well as K) is very sensitive
	// to small changes in masses and compilation flags.
	//-----------------
	// Note that changing coordinates by those at ala3/plop_ala.coords using
	// createUnitsFromFileChangingCoordinates can increase precision to 10-12
	TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/H.txt",
		not SKIP_OXT,
		1e-5);// allowed relative error (with respect to the maximum absolute value)

	TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/H.txt",
		not SKIP_OXT,
		2e-5);

	TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/H.txt",
		not SKIP_OXT,
		1e-5);

	// With prior test checks this one was passing, but i will not, as INMA was using OXT in this case
	FAILING_TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/H.txt",
		SKIP_OXT,
		1e-5);

	TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/H.txt",
		not SKIP_OXT,
		1e-5);

	TEST_FUNCTION(testCalculateH,
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/H.txt",
		not SKIP_OXT,
		2e-5);

	finish();
}

void TestANMInternalHessianFunctions::finish() {
	Test::finish();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of Dij matrix
///
/// \param golde_file [In] file containing the correct results
/// of the algorithm
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \author vgil
/// \date 20/11/2012
///////////////////////////////////////////////////////////////
bool TestANMInternalHessianFunctions::testCalculateD(const char* golden_file){
	vector<vector<double> > data;
	TestTools::load_vector_of_vectors(data, golden_file);
	// Golden data has, for each Dij submatrix:
	// Point i coordinates
	// Point j coordinates
	// Each Dij, as a 6x6 matrix
	bool ok = true;
	unsigned int number_of_dijs = data.size() / 8;
	double dij[6*6];
	for (unsigned int i = 0; i < number_of_dijs; ++i){
		int current_start = i*8;
		Point r_i(Utils::vectorToPointer(data[current_start]));
		Point r_j(Utils::vectorToPointer(data[current_start+1]));
		// We multiply the constant by sq. matrix as it is going to be divided
		// later (it wasn't at the time of getting the golden data)
		ANMICHessianCalculator::calculateDij(1*r_i.squaredDistance(r_j), r_i, r_j, dij);
		for (unsigned int j = 0; j < 6; j++) {
			ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(Utils::vectorToPointer(data[current_start+2+j]),
					&dij[j*6], 6, 1e-12);
		}
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of Tab matrix
///
/// \param complex_pdb [In]
/// \param unit_atoms_coords [In]
/// \param goldenT [In] file containing the correct results
/// of the algorithm
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \author vgil
/// \date 01/12/2012
///////////////////////////////////////////////////////////////
bool TestANMInternalHessianFunctions::testCalculateT(const char* complex_pdb,
		const char* goldenT,
		bool skip_OXT,
		double test_precision){

	System sys;
	vector<Unit*> units;
	Complex* complex;

	cout <<"Testing "<<complex_pdb<<endl;

	TestANMICTools::createUnitsFromFile(complex_pdb,
			units,
			complex,
			skip_OXT);

	vector<vector<double> > data;
	TestTools::load_vector_of_vectors(data, goldenT);
	// Golden data has, for each of the Tab submatrices:
	// int(a) int(b) -> index of the a and b units
	// Each Tab, as a 6x6 matrix
	bool ok = true;
	unsigned int number_of_Tabs = data.size() / 7;

	ElasticConstantCalculator* ecc = new InverseExponentialElasticConstant(1, 3.8, 6);

	for (unsigned int i = 0; i < number_of_Tabs; ++i){
		double tab[6*6] = {
			0 ,0 ,0, 0, 0, 0,
			0 ,0 ,0, 0, 0, 0,
			0 ,0 ,0, 0, 0, 0,
			0 ,0 ,0, 0, 0, 0,
			0 ,0 ,0, 0, 0, 0,
			0 ,0 ,0, 0, 0, 0,
		};
		int current_start = i*7;
		int a = (int) data[current_start][0];
		int b = (int) data[current_start][1]+1; // Correction from golden file to ours (dihedrals vs unit number)
		cout<<"DBG Number of units "<<units.size()<<endl;
		cout<<"DBG a "<<a<<" b "<<b<<endl;

		ANMICHessianCalculator::calculateTab(units[a], units[b], 100, ecc, tab, skip_OXT);
		for (unsigned int j = 0; j < 6; j++) {
			ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(Utils::vectorToPointer(data[current_start+1+j]),
					&tab[j*6], 6, test_precision);
		}
	}

	Utils::clearVector<Unit>(units);
	delete complex;
	delete ecc;

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of the U matrix
///
/// \param complex_pdb [In]
/// \param unit_atoms_coords [In]
/// \param goldenU [In] file containing the correct results
/// of the algorithm
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \author vgil
/// \date 11/12/2012
///////////////////////////////////////////////////////////////
bool TestANMInternalHessianFunctions::testCalculateU(const char* complex_pdb,
		const char* goldenU,
		bool skip_OXT,
		double test_precision) {

	System sys;
	vector<Unit*> units;
	Complex* complex;

	TestANMICTools::createUnitsFromFile(complex_pdb,
				units,
				complex,
				skip_OXT);

	vector<vector<double> > data;
	TestTools::load_vector_of_vectors(data, goldenU);
	// Golden data has, for each of the Uab submatrices:
	// int(a) int(b) -> index of the a and b units
	// Each Uab, as a 6x6 matrix

	ElasticCalculatorMock ecc(1);
	std::vector<std::vector<double> > U = ANMICHessianCalculator::calculateU(100, &ecc, units, false);

	double uab[6*6];
	bool ok = true;
	unsigned int number_of_Uabs = data.size() / 7;
	for (unsigned int i = 0; i < number_of_Uabs; ++i){
		int current_start = i*7;
		int a = (int) data[current_start][0];
		int b = (int) data[current_start][1];
		ANMICMath::getABMatrix(U, uab, a, b);
		for (unsigned int j = 0; j < 6; j++) {
			ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(Utils::vectorToPointer(data[current_start+1+j]),
					&uab[j*6], 6, test_precision);
		}
	}

	Utils::clearVector<Unit>(units);
	delete complex;

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of the H matrix
///
/// \param complex_pdb [In]
/// \param unit_atoms_coords [In]
/// \param goldenH [In] file containing the correct results
/// of the algorithm
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \author vgil
/// \date 11/12/2012
///////////////////////////////////////////////////////////////
bool TestANMInternalHessianFunctions::testCalculateH(
		const char* complex_pdb,
		const char* goldenH,
		bool skip_OXT,
		double test_tolerance){

	cout << "Testing H with "<<complex_pdb<<endl;
	System sys;
	vector<Unit*> units;
	Complex* complex;

	cout <<"Testing "<<complex_pdb<<endl;

	TestANMICTools::createUnitsFromFile(complex_pdb,
					units,
					complex,
					skip_OXT);

	vector<vector<double> > expected;
	TestTools::load_vector_of_vectors(expected, goldenH);

	InverseExponentialElasticConstant ecc(1, 3.8, 6);
	std::vector<std::vector<double> > U = ANMICHessianCalculator::calculateU(100, &ecc, units, skip_OXT);

	TriangularMatrix* H = ANMICHessianCalculator::calculateH(units, U);

	bool ok = true;
	for (unsigned int i = 0; i < expected.size(); ++i){
		vector<double> h_row;
		for(unsigned int j = 0; j < expected.size();++j){
			if(j>=i){
				h_row.push_back((*H)(i,j));
			}
			else{
				h_row.push_back(0.0);
			}
		}
		ok &= Assertion::expectedVectorHasRelativeErrorNotBiggerThan(
							expected[i],
							h_row,
							test_tolerance);
	}

	Utils::clearVector<Unit>(units);
	delete complex;
	delete H;
	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the...
///
/// \param initialH [In] calculated H matrix
/// \param goldenHmodif [In] file containing the correct results
/// of the algorithm
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \author vgil
/// \date 11/12/2012
///////////////////////////////////////////////////////////////
bool TestANMInternalHessianFunctions::testModifyHessianWithExtraTorsion(const char* initialH, const char* goldenHmodif) {

	vector<vector<double> > Hv, expected;
	TestTools::load_vector_of_vectors(Hv, initialH);
	TestTools::load_vector_of_vectors(expected, goldenHmodif);

	TriangularMatrix H(expected.size(),expected.size());
	for(unsigned int i = 0; i< expected.size();++i){
		for(unsigned int j = i; j < expected.size();++j){
			H(i,j) = Hv[i][j];
		}
	}

	ANMICHessianCalculator::modifyHessianWithExtraTorsion(&H);

	bool ok = true;
	for (unsigned int i = 0; i < expected.size(); ++i){
		for(unsigned int j = i; j < expected.size();++j){
			ok &= Assertion::expectedEqualsCalculatedWithinPrecision(
					expected[i][j],
					H(i,j),
					1e-14);
		}
	}

	return ok;
}
