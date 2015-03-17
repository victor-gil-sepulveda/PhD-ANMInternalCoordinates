/////////////////////////////////////////////////////////////////////////////
/// TestDihedralConstraintFunctions.cpp
///
/// Implementation of TestTemplate class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \\author arincon
/// \date 25/05/2012
/////////////////////////////////////////////////////////////////////////////

#include "TestDihedralConstraintFunctions.h"

#include <string>
#include "../HarmonicDihedralConstraintFunctions.h"
#include "../../../../../../Tools/stringTools.h"
#include "../../../../../../Tools/Assertion.h"
#include "../../../../../../Tools/Math/MathTools.h"

using namespace std;
using StringTools::split;
using StringTools::toDouble;

TestDihedralConstraintFunctions::TestDihedralConstraintFunctions(string name)
{
    test_name = name;
}

TestDihedralConstraintFunctions::~TestDihedralConstraintFunctions()
{}

void TestDihedralConstraintFunctions::init()
{}

void TestDihedralConstraintFunctions::run()
{
    Test::run();

    TEST_FUNCTION(testDihedralAngleWithArcTanFunction, 1e-12)
	TEST_REGRESSION_FUNCTION(testEnergyWithArcTanFunction, 1.0e-12)
	TEST_REGRESSION_FUNCTION(testDerivatives, 1.0e-4) // precission lowered (golden data precission was also lowered)
	TEST_REGRESSION_FUNCTION(testHessian, 1.0e-3) // Idem


	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(30), Math::degToRad(30))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(150), Math::degToRad(150))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(200), Math::degToRad(200))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(330), Math::degToRad(330))

	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(-30), Math::degToRad(330))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(-150), Math::degToRad(210))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(-200), Math::degToRad(160))
	TEST_FUNCTION(testPIto2PIRange, Math::degToRad(-330), Math::degToRad(30))

	TEST_FUNCTION(testAngleFixing, Math::degToRad(30), Math::degToRad(30))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(150), Math::degToRad(150))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(200), Math::degToRad(-160))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(330), Math::degToRad(-30))

	TEST_FUNCTION(testAngleFixing, Math::degToRad(-30), Math::degToRad(-30))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(-150), Math::degToRad(-150))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(-200), Math::degToRad(160))
	TEST_FUNCTION(testAngleFixing, Math::degToRad(-330), Math::degToRad(30))
    TEST_FUNCTION(testAngleFixing, 11.09, -1.476370614359173)

	TEST_FUNCTION(testRadSubtraction, Math::degToRad(40), Math::degToRad(30), Math::degToRad(10))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(30), Math::degToRad(40), Math::degToRad(-10))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(180), Math::degToRad(30), Math::degToRad(150))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(30), Math::degToRad(180), Math::degToRad(-150))

	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-40), Math::degToRad(-30), Math::degToRad(-10))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-30), Math::degToRad(-40), Math::degToRad(10))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-180), Math::degToRad(-30), Math::degToRad(-150))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-30), Math::degToRad(-180), Math::degToRad(150))

	TEST_FUNCTION(testRadSubtraction, Math::degToRad(170), Math::degToRad(-30), Math::degToRad(-160))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-30), Math::degToRad(170), Math::degToRad(160))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-170), Math::degToRad(30), Math::degToRad(160))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(30), Math::degToRad(-170), Math::degToRad(-160))

	//Final difference is > 180 or < -180 (we must pick the smaller angle)
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(80), Math::degToRad(-80), Math::degToRad(160))
	TEST_FUNCTION(testRadSubtraction, Math::degToRad(-80), Math::degToRad(80), Math::degToRad(-160))

    finish();
}

void TestDihedralConstraintFunctions::finish()
{
    Test::finish();
}

bool TestDihedralConstraintFunctions::testRadSubtraction(double a, double b, double expected_increment){
	return Assertion::expectedEqualsCalculatedWithinPrecision(expected_increment,
				HarmonicDihedralConstraintFunctions::rad_subtraction(a,b),
				1e-12);
}


bool TestDihedralConstraintFunctions::testPIto2PIRange(double angle, double expected){
	return Assertion::expectedEqualsCalculatedWithinPrecision(expected,
			HarmonicDihedralConstraintFunctions::to_0_2PI_range(angle),
			1e-12);
}

bool TestDihedralConstraintFunctions::testAngleFixing(double angle, double expected){
	return Assertion::expectedEqualsCalculatedWithinPrecision(expected,
			HarmonicDihedralConstraintFunctions::put_in_pi_minus_pi_range(angle),
			1e-12);
}
/////////////////////////////////////////////////////////////////////////////
/// \remarks
/// This test checks that the dihedral angle is correct
///
/// \param precision [In] To consider two doubles as equal, its difference must be less than this value.
///
/// \return True if the test passed, false if not
///
/// \author arincon
/// \date 13/10/2012
/////////////////////////////////////////////////////////////////////////////
bool TestDihedralConstraintFunctions::testDihedralAngleWithArcTanFunction(double precision) {
	bool passed = true;

	// A file with test data is loaded in lines vector
	vector<string> lines;
	TestTools::load_vector(lines, "src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/Tests/data/test_angle.txt");

	// For each element of the line vector, a test for each interaction is performed
	unsigned int size = lines.size();
	for(unsigned int i=0; i<size; i++) {
		vector<string> test = split(lines[i], ",");

		double expected = toDouble(test[0]);
		double result = HarmonicDihedralConstraintFunctions::calculateDihedralAngleWithArcTanFunction(
				toDouble(test[1]), toDouble(test[2]),
				toDouble(test[3]), toDouble(test[4]),
				toDouble(test[5]), toDouble(test[6]),
				toDouble(test[7]), toDouble(test[8]),
				toDouble(test[9]), toDouble(test[10]),
				toDouble(test[11]), toDouble(test[12]));


		passed &= Assertion::expectedEqualsCalculatedWithinPrecision(expected, result, precision);
	}

	return passed;
}

/////////////////////////////////////////////////////////////////////////////
/// \remarks
/// This test checks that the dihedral constraint energy functions are correct
///
/// \param precision [In] To consider two doubles as equal, its difference must be less than this value.
///
/// \return True if the test passed, false if not
///
/// \\author arincon
/// \date 18/10/2012
/////////////////////////////////////////////////////////////////////////////
bool TestDihedralConstraintFunctions::testEnergyWithArcTanFunction(double precision) {
	bool passed = true;

	// A file with test data is loaded in lines vector
	vector<string> lines;
	TestTools::load_vector(lines, "src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/Tests/data/test_energia.txt");

	// For each element of the line vector, a test for each interaction is performed
	unsigned int size = lines.size();
	for(unsigned int i=0; i<size; i++) {
		vector<string> test = split(lines[i], ",");
		double expected = toDouble(test[0]);
		double result = HarmonicDihedralConstraintFunctions::calculateEnergyWithArcTanFunction(toDouble(test[1]), toDouble(test[2]),
																		toDouble(test[3]), toDouble(test[4]),
																		toDouble(test[5]), toDouble(test[6]),
																		toDouble(test[7]), toDouble(test[8]),
																		toDouble(test[9]), toDouble(test[10]),
																		toDouble(test[11]), toDouble(test[12]),
																		toDouble(test[13]), toDouble(test[14]));

		passed &= Assertion::expectedEqualsCalculatedWithinPrecision(expected, result, precision);
		//cout << expected << " == " << result << " passed? " << passed << endl;
	}

	return passed;
}

/////////////////////////////////////////////////////////////////////////////
/// \remarks
/// This test checks that the dihedral constraint derivative functions are correct
///
/// \param precision [In] To consider two doubles as equal, its difference must be less than this value.
///
/// \return True if the test passed, false if not
///
/// \\author arincon
/// \date 25/05/2012
/////////////////////////////////////////////////////////////////////////////
bool TestDihedralConstraintFunctions::testDerivatives(double precision) {

	typedef double (*fptr)(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);

	bool passed = true;

	// A vector containing all the derivatives functions is created
	std::vector<fptr> f;
	f.push_back(HarmonicDihedralConstraintFunctions::derivativeXa);
	f.push_back(HarmonicDihedralConstraintFunctions::derivativeYa);
	f.push_back(HarmonicDihedralConstraintFunctions::derivativeZa);
	f.push_back(HarmonicDihedralConstraintFunctions::derivativeXb);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeYb);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeZb);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeXc);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeYc);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeZc);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeXd);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeYd);
  	f.push_back(HarmonicDihedralConstraintFunctions::derivativeZd);

  	// A file with test data is loaded in lines vector
	vector<string> lines;
  	TestTools::load_vector(lines, "src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/Tests/data/test_derivatives.reg.txt");

  	// For each element of the line vector, a test for each interaction is performed
  	unsigned int size = lines.size();
  	if(lines.empty()) passed = false;

  	for(unsigned int i=0; i<size; i++) {
  		// Each partial derivative needs 12 lines.
  		vector<string> test = split(lines[i], " ");
  		double expected = toDouble(test[0]);
  		double result = f[i%12](toDouble(test[1]), toDouble(test[2]), toDouble(test[3]),
  						toDouble(test[4]), toDouble(test[5]), toDouble(test[6]),
  						toDouble(test[7]), toDouble(test[8]), toDouble(test[9]),
  						toDouble(test[10]), toDouble(test[11]), toDouble(test[12]),
  						toDouble(test[13]), toDouble(test[14]));
  		passed = Assertion::expectedEqualsCalculatedWithinPrecision(expected, result, precision) and passed;

  	}

  	if(not passed){
  		cout<<"Problem with individual interactions!"<<endl;
  	}
  	else{
  		cout<<"Individual interactions OK"<<flush<<endl;
  	}


  	// Test cases for calculateGradient are run
  	unsigned int lines_per_test_case = 12;
  	unsigned int num_test_cases = lines.size() / lines_per_test_case;

  	for(unsigned int i=0; i<num_test_cases; i++) {
		vector<double> expected_values;
		vector<double> calculated_values;
		vector<string> test;

  		for(unsigned int j=0; j<lines_per_test_case; j++) {
  			test = split(lines[i*lines_per_test_case+j], " ");
  			double expected = toDouble(test[0]);
			expected_values.push_back(expected);
  		}

  		calculated_values = HarmonicDihedralConstraintFunctions::calculateGradient(toDouble(test[1]), toDouble(test[2]),
																				toDouble(test[3]), toDouble(test[4]),
																				toDouble(test[5]), toDouble(test[6]),
																				toDouble(test[7]), toDouble(test[8]),
																				toDouble(test[9]), toDouble(test[10]),
																				toDouble(test[11]), toDouble(test[12]),
																				toDouble(test[13]), toDouble(test[14]));
  		passed = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_values, calculated_values, precision)
  							and passed;
  	}

  	if(not passed){
  		cout<<"Problem with calculateGradient!"<<endl;
  	}

  	return passed;
}

/////////////////////////////////////////////////////////////////////////////
/// \remarks
/// This test checks that the dihedral constraint hessian functions are correct
///
/// \param precision [In] To consider two doubles as equal, its difference must be less than this value.
///
/// \return True if the test passed, false if not
///
/// \author arincon
/// \date 25/05/2012
/////////////////////////////////////////////////////////////////////////////
bool TestDihedralConstraintFunctions::testHessian(double precision) {
	bool passed = true;

	// A file with test data is loaded in lines vector
	vector<string> lines;
	TestTools::load_vector(lines, "src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/Tests/data/test_hessian.reg.txt");

	// For each element of the line vector, a test for each interaction is performed
	unsigned int size = lines.size();
	for(unsigned int i=0; i<size; i++) {
		vector<string> test = split(lines[i], " ");
		vector<double> hessian = HarmonicDihedralConstraintFunctions::calculateHessian(toDouble(test[0]), toDouble(test[1]), toDouble(test[2]), toDouble(test[3]), toDouble(test[4]), toDouble(test[5]), toDouble(test[6]), toDouble(test[7]), toDouble(test[8]), toDouble(test[9]), toDouble(test[10]), toDouble(test[11]), toDouble(test[12]), toDouble(test[13]));
		unsigned int k = 0;
		for(unsigned int j=14; j<hessian.size(); j++) {
			passed &= Assertion::expectedEqualsCalculatedWithinPrecision(toDouble(test[j]), hessian[k], precision);
			k++;
		}
	}

	return passed;
}

