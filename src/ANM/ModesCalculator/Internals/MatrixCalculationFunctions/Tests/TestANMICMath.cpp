/////////////////////////////////////////////////////////////////////////////
/// TestANMICHessianCalculation.cpp
///
/// Implementation of TestANMICHessianCalculation class
///
/// \author arincon
/// \date 18/11/2012
/////////////////////////////////////////////////////////////////////////////

#include "TestANMICMath.h"
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

using namespace std;

TestANMICMath::TestANMICMath(string name) {
	test_name = name;
}

TestANMICMath::~TestANMICMath() {
}

void TestANMICMath::init() {
}

void TestANMICMath::run() {
	Test::run();

	TEST_FUNCTION(testMultiplyColumnByRow);
	TEST_FUNCTION(testMultiplyMatrixByScalar);
  	TEST_FUNCTION(testMultiplyMatrixByScalarGeneral);
  	TEST_FUNCTION(testAddMatrices);
	TEST_FUNCTION(testAddMatricesGeneral);
	TEST_FUNCTION(testAddMatricesInPlace);
	TEST_FUNCTION(testAddMatricesInPlaceGeneral);
	TEST_FUNCTION(testSubtractMatrices);
	TEST_FUNCTION(testSubtractMatricesGeneral);
  	TEST_FUNCTION(testMultitplyRowByMatrix);
  	TEST_FUNCTION(testMultiplyRowByColumn);
  	TEST_FUNCTION(testMultiplyMatrixByMatrix);
  	TEST_FUNCTION(testMultiplyMatrixByMatrixGeneral);
  	TEST_FUNCTION(testMultiplyIMatrixByEVector);
  	TEST_FUNCTION(testInvertMatrix);
	TEST_FUNCTION(testInvertMatrixGeneral);
	TEST_FUNCTION(testTransposeMatrix);
	TEST_FUNCTION(testTransposeMatrixGeneral);
	TEST_FUNCTION(testGetABMatrix);
	TEST_FUNCTION(testAddMatrixToU);
	TEST_FUNCTION(testConvertMatrixToArray);
	TEST_FUNCTION(testArbitrarySizeMatrixMult);
	TEST_FUNCTION(testTranspose)

	finish();
}

void TestANMICMath::finish() {
	Test::finish();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of a column matrix by a row matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyColumnByRow() {
	int size = 6;
	double a[6] = { 2, 3, 1, 9, 4, 2 };
	double b[6] = { 1, 2, 3, 4, 3, 1 };
	double result[6*6];

	double expected[6*6] = {
		2.000,  4.000,  6.000,  8.000,  6.000, 2.000,
		3.000,  6.000,  9.000, 12.000,  9.000, 3.000,
		1.000,  2.000,  3.000,  4.000,  3.000, 1.000,
		9.000, 18.000, 27.000, 36.000, 27.000, 9.000,
		4.000,  8.000, 12.000, 16.000, 12.000, 4.000,
		2.000,  4.000,  6.000,  8.000,  6.000, 2.000
	};

	ANMICMath::multiplyColumnByRow(a, b, result, size);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 6, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyMatrixByMatrix() {
	double first[3][3] = {
		{1, 2, 9},
		{3, 8, 2},
		{4, 3, 7},
	};
	double second[3][3] = {
		{2, 2, 3},
		{0, 1, 0},
		{0, 1, 5},
	};

	double expected[3][3] = {
		{2, 13, 48},
		{6, 16, 19},
		{8 ,18, 47},
	};

	double result[3][3];

	ANMICMath::multiplyMatrixByMatrix(first, second, result);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 3; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], result[i], 3, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyMatrixByMatrixGeneral() {
	double first[3*3] = {
		1, 2, 9,
		3, 8, 2,
		4, 3, 7,
	};
	double second[3*3] = {
		2, 2, 3,
		0, 1, 0,
		0, 1, 5,
	};

	double expected[3*3] = {
		2, 13, 48,
		6, 16, 19,
		8 ,18, 47,
	};

	double result[3*3];

	ANMICMath::multiplyMatrixByMatrix(first, second, result, 3);

	// compares the result with the expected matrix
	bool ok = true;
	ok = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 3*3, 1e-13);

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of a matrix by a scalar
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyMatrixByScalar() {
	double result[6][6] = {
		{ 2.000,  4.000,  6.000,  8.000,  6.000, 2.000 },
		{ 3.000,  6.000,  9.000, 12.000,  9.000, 3.000 },
		{ 1.000,  2.000,  3.000,  4.000,  3.000, 1.000 },
		{ 9.000, 18.000, 27.000, 36.000, 27.000, 9.000 },
		{ 4.000,  8.000, 12.000, 16.000, 12.000, 4.000 },
		{ 2.000,  4.000,  6.000,  8.000,  6.000, 2.000 }
	};

	double expected[6][6] = {
		{  6.000, 12.000, 18.000,  24.000, 18.000,  6.000 },
		{  9.000, 18.000, 27.000,  36.000, 27.000,  9.000 },
		{  3.000,  6.000,  9.000,  12.000,  9.000,  3.000 },
		{ 27.000, 54.000, 81.000, 108.000, 81.000, 27.000 },
		{ 12.000, 24.000, 36.000,  48.000, 36.000, 12.000 },
		{  6.000, 12.000, 18.000,  24.000, 18.000,  6.000 }
	};

	ANMICMath::multiplyMatrixByScalar(3, result);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 2; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], result[i], 6, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of a matrix by a scalar
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyMatrixByScalarGeneral() {
	int length = 6;
	int scalar = 3;
	double result[6*6] = {
		2.000,  4.000,  6.000,  8.000,  6.000, 2.000,
		3.000,  6.000,  9.000, 12.000,  9.000, 3.000,
		1.000,  2.000,  3.000,  4.000,  3.000, 1.000,
		9.000, 18.000, 27.000, 36.000, 27.000, 9.000,
		4.000,  8.000, 12.000, 16.000, 12.000, 4.000,
		2.000,  4.000,  6.000,  8.000,  6.000, 2.000
	};

	double expected[6*6] = {
	 	 6.000, 12.000, 18.000,  24.000, 18.000,  6.000,
		 9.000, 18.000, 27.000,  36.000, 27.000,  9.000,
		 3.000,  6.000,  9.000,  12.000,  9.000,  3.000,
		27.000, 54.000, 81.000, 108.000, 81.000, 27.000,
		12.000, 24.000, 36.000,  48.000, 36.000, 12.000,
		 6.000, 12.000, 18.000,  24.000, 18.000,  6.000
	};

	ANMICMath::multiplyMatrixByScalar(scalar, result, length);

	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, length*length, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the addition of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 02/12/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testAddMatrices() {
	double result[6][6];
	double a[6][6] = {
		{ 1, 1, 1, 1, 1, 1 },
		{ 1, 1, 1, 1, 1, 1 },
		{ 1, 1, 1, 1, 1, 1 },
		{ 1, 1, 1, 1, 1, 1 },
		{ 1, 1, 1, 1, 1, 1 },
		{ 1, 1, 1, 1, 1, 1 }
	};

	double expected[6][6] = {
		{ 2, 2, 2, 2, 2, 2 },
		{ 2, 2, 2, 2, 2, 2 },
		{ 2, 2, 2, 2, 2, 2 },
		{ 2, 2, 2, 2, 2, 2 },
		{ 2, 2, 2, 2, 2, 2 },
		{ 2, 2,	2, 2, 2, 2 }
	};

	ANMICMath::addMatrices(a, a, result);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 6; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], result[i], 6, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the addition of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testAddMatricesGeneral() {
	int length = 6;
	double result[6*6];
	double a[6*6] = {
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1
	};

	double expected[6*6] = {
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 2, 2
	};

	ANMICMath::addMatrices(a, a, result, length);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, length*length, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the addition of two matrices and store the result
/// in the first one
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testAddMatricesInPlace() {
	double a[3][3] = {
		{ 1, 1, 1 },
		{ 1, 1, 1 },
		{ 1, 1, 1 },
	};

	double expected[3][3] = {
		{ 2, 2, 2 },
		{ 2, 2, 2 },
		{ 2, 2, 2 },
	};

	ANMICMath::addMatricesInPlace(a, a);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 3; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], a[i], 3, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the addition of two matrices and store the result
/// in the first one
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testAddMatricesInPlaceGeneral() {
	int length = 3;
	double a[3*3] = {
		1, 1, 1,
		1, 1, 1,
		1, 1, 1,
	};

	double expected[3*3] = {
		2, 2, 2,
		2, 2, 2,
		2, 2, 2,
	};

	ANMICMath::addMatricesInPlace(a, a, 3);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, a, length*length, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the addition of a 6x6 matrix to the indicated
/// positions of U matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testAddMatrixToU() {
	int a = 2;
	int b = 3;
	double expected[6*6];
	double result[6*6];
	double matrix[6*6] = {
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1
	};

	// initialization of U matrix
	std::vector<std::vector<double> > UMatrix ;
	for (int i =0; i < 36; i++){
		UMatrix.push_back(vector<double>(36,0));
	}

	ANMICMath::getABMatrix(UMatrix, expected, a, b);

	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			expected[i+j*6] += 1;
		}
	}

	ANMICMath::addMatrixToU(UMatrix, matrix, a, b);

	ANMICMath::getABMatrix(UMatrix, result, a, b);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 6*6, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the subtraction of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testSubtractMatrices() {
	double result[3][3];
	double a[3][3] = {
		{ 1, 1, 1 },
		{ 1, 1, 1 },
		{ 1, 1, 1 }
	};

	double expected[3][3] = {
		{ 0, 0, 0 },
		{ 0, 0, 0 },
		{ 0, 0, 0 }
	};

	ANMICMath::subtractMatrices(a, a, result);

	// compares the result with the expected matrix
	bool ok;
	for (unsigned int i = 0; i < 3; i++) {
		ok = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], result[i], 3, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the subtraction of two matrices
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testSubtractMatricesGeneral() {
	int length = 6;
	double result[6*6];
	double a[6*6] = {
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1
	};

	double expected[6*6] = {
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
	};

	ANMICMath::subtractMatrices(a, a, result, length);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, length*length, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests if the indicated submatrix of U matrix is correct
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testGetABMatrix() {
	int a = 0;
	int b = 5;
	double result[6*6];
	double expected[6*6] = {
		142.010526355199,  31.609879942883, 132.811573223004,   39.386759703363,    4.512907899362,  -60.216261721861,
		 31.609879942883,  55.025110625176,  42.550012308292,  -55.224205531400,   -5.361755577275,   53.502132648078,
		132.811573223004,  42.550012308292, 131.412353089313,   16.369528954957,   -0.502612532484,  -34.025004126087,
		 39.386759703363, -55.224205531400,  16.369528954957,  293.252204768955,   91.428813686419, -325.352979906807,
		  4.512907899362,  -5.361755577275,  -0.502612532484,   91.428813686419,   38.058212139223, -104.022541630676,
		-60.216261721861,  53.502132648078, -34.025004126087, -325.352979906807, -104.022541630676,  364.957468651242
	};

	// loads the uab data into the U matrix
	std::vector<std::vector<double> > U_data;
	TestTools::load_vector_of_vectors(U_data, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/U.txt");
	std::vector<std::vector<double> > U(36,std::vector<double>(36,0));
	for (unsigned int i = 0; i< U_data.size(); i+=7){
		int subindex_a = U_data[i][0];
		int subindex_b = U_data[i][1];
		for (unsigned int j = 0; j< 6; ++j ){
			for (unsigned int k = 0; k < 6; ++k ){
				U[subindex_a*6+j][subindex_b*6+k] = U_data[i+1+j][k];
			}
		}
	}

	ANMICMath::getABMatrix(U, result, a, b);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 6*6, 1e-12);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of row matrix by a matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultitplyRowByMatrix() {
	double a[6] = { 1.000, 2.000, 3.000, 1.000, 2.000, 3.000 };

	double Rab[6*6] = {
		2.000, 0.000, 4.000, 1.000, 2.000, 2.000,
		1.000, 3.000, 4.000, 0.000, 2.000, 1.000,
		2.000, 0.000, 4.000, 1.000, 2.000, 2.000,
		1.000, 3.000, 4.000, 0.000, 2.000, 1.000,
		1.000, 3.000, 4.000, 0.000, 2.000, 1.000,
		5.000, 3.000, 4.000, 1.000, 2.000, 4.000
	};

	double result[6];
	double expected[6] = { 28, 24, 48, 7, 24, 25 };

	ANMICMath::multitplyRowByMatrix(a, Rab, result, 6);

	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 6, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of row matrix by a matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyRowByColumn() {
	double result;
	double a[6] = { 1, 2, 3, 1, 0, 1 };
	double b[6] = { 2, 1, 0, 0, 1, 1 };
	double expected = 5;
	int size = 6;

	result = ANMICMath::multiplyRowByColumn(a, b, size);

	return Assertion::expectedEqualsCalculatedWithinPrecision(expected, result, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the multiplication of the I matrix by the E vector
///
/// \return True if the resulting matrix is correct
///
/// \author vgil
/// \date 10/08/2012
///////////////////////////////////////////////////////////////
bool TestANMICMath::testMultiplyIMatrixByEVector() {
	double e[3] = { 1.000, 0.000, 1.000 };

	double I[3][3] = {
		{ 1.000, 0.000, 2.000 },
		{ 0.000, 1.000, 3.000 },
		{ 2.000, 2.000, 3.000 }
	};

	double result[3];
	double expected[3] = { 3, 3, 5 };

	Point res = ANMICMath::multiplyIMatrixByEVector(I, e);
	result[0] = res.getX();
	result[1] = res.getY();
	result[2] = res.getZ();

	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 3, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the inversion of a matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 10/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testInvertMatrix() {
	double ident[3][3];
	double matrix[3][3] = {
		{1, 2, 3},
		{0, 1, 4},
		{5, 6, 0},
	};

	double result[3][3];

	double expected[3][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1},
	};

	ANMICMath::invertIMatrix(matrix, result);

	ANMICMath::multiplyMatrixByMatrix(matrix, result, ident);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 3; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], ident[i], 3, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the inversion of a matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 10/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testInvertMatrixGeneral() {
	double ident[3*3];
	double matrix[3*3] = {
		1, 2, 3,
		0, 1, 4,
		5, 6, 0,
	};

	double result[3*3];

	double expected[3*3] = {
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};

	ANMICMath::invertMatrix(matrix, result, 3);

	ANMICMath::multiplyMatrixByMatrix(matrix, result, ident, 3);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, ident, 3*3, 1e-13);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of the transpose matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testTransposeMatrix() {
	double matrix[3][3] = {
		{1, 2, 3},
		{0, 1, 4},
		{5, 6, 0},
	};

	double expected[3][3] = {
		{1, 0, 5},
		{2, 1, 6},
		{3, 4, 0},
	};

	double result[3][3];

	ANMICMath::transposeMatrix(matrix, result);

	// compares the result with the expected matrix
	bool ok = true;
	for (unsigned int i = 0; i < 3; i++) {
		ok &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected[i], result[i], 3, 1e-13);
	}

	return ok;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Tests the calculation of the transpose matrix
///
/// \return True if the resulting matrix is correct
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
bool TestANMICMath::testTransposeMatrixGeneral() {
	double matrix[3*3] = {
		1, 2, 3,
		0, 1, 4,
		5, 6, 0,
	};

	double expected[3*3] = {
		1, 0, 5,
		2, 1, 6,
		3, 4, 0,
	};

	double result[3*3];

	ANMICMath::transposeMatrix(matrix, result, 3);

	// compares the result with the expected matrix
	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 3*3, 1e-13);
}

bool TestANMICMath::testConvertMatrixToArray() {
	double *result = new double[9];

	double **in = new double*[3];
	in[0] = new double[3];
	in[1] = new double[3];
	in[2] = new double[3];

	int k = 0;
	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			in[i][j] = k;
			k++;
		}
	}

	double *expected = new double[9];
	for(int i=0; i<9; i++) {
		expected[i] = i;
	}

	ANMICMath::convertMatrixToArray(in, result, 3, 3);
	bool ok = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected, result, 9, 1e-13);

	delete [] result;
	delete [] in[0];
	delete [] in[1];
	delete [] in[2];
	delete [] in;
	delete [] expected;

	return ok;
}

bool TestANMICMath::testArbitrarySizeMatrixMult(){
	vector<vector<double> > matrix1,matrix2, out1,out2, expected_out1, expected_out2;

	TestTools::load_vector_of_vectors(matrix1, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/math_mult_1.txt");
	TestTools::load_vector_of_vectors(matrix2, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/math_mult_2.txt");
	TestTools::load_vector_of_vectors(expected_out1, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/math_mult_result.txt");
	TestTools::load_vector_of_vectors(expected_out2, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/math_mult_result2.txt");

	ANMICMath::multiplyMatrixByMatrix(matrix1,matrix2,out1);
	ANMICMath::multiplyMatrixByMatrix(matrix2,matrix1,out2);

	bool ok = true;
	for (unsigned int i =0; i<expected_out1.size();++i){
		ok = ok && Assertion::expectedVectorEqualsCalculatedWithinPrecision(
				expected_out1[i],
				out1[i],
				1e-16);
	}

	for (unsigned int i =0; i<expected_out2.size();++i){
		ok = ok && Assertion::expectedVectorEqualsCalculatedWithinPrecision(
				expected_out2[i],
				out2[i],
				1e-16);
	}

	return ok;
}

bool TestANMICMath::testTranspose(){
	vector<vector<double> > matrix, transposed, expected_transposed;
	TestTools::load_vector_of_vectors(matrix, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/math_mult_1.txt");
	TestTools::load_vector_of_vectors(expected_transposed, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/math/transposed.txt");
	ANMICMath::transpose(matrix, transposed);
	bool ok = true;
	for (unsigned int i =0; i<expected_transposed.size();++i){
		ok = ok && Assertion::expectedVectorEqualsCalculatedWithinPrecision(
				expected_transposed[i],
				transposed[i],
				1e-16);
	}
	return ok;
}
