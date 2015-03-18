/////////////////////////////////////////////////////////////////////////////
/// TestUnitsBuilder.cpp
///
/// Implementation of TestUnitsBuilder class
///
/// \author vgil
/// \date 06/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "TestUnitsBuilder.h"
#include "../Unit.h"
#include "../UnitsBuilder.h"
#include "../../../../../Tools/Math/Point.h"
#include "../../../../../Molecules/ComplexBuilder.h"
#include "../../../../../Molecules/Selection/Selector.h"
#include "../../../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../../../Tools/Assertion.h"
#include "../../../../../Tools/TestTools.h"
#include "../../../../../Tools/Utils.h"
#include "../../../../../System/System.h"
#include "../../../../../Molecules/Complex.h"
#include "../../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"
#include "../../MatrixCalculationFunctions/Tests/TestANMICTools.h"

#include <vector>
using namespace std;

TestUnitsBuilder::TestUnitsBuilder(string name){
    test_name = name;
}

TestUnitsBuilder::~TestUnitsBuilder(){}

void TestUnitsBuilder::init(){
}

void TestUnitsBuilder::run(){
    Test::run();

    TEST_FUNCTION(testUnitaryBondVectors,
			"src/Molecules/Test/Data/ProteinResiduesSamples/plop_ala.pdb",
			"src/ANM/ModesCalculator/Internals/CoarseGrainModel/Test/data/unitary_trialanine.txt")

    TEST_FUNCTION(testNunitsCenterAndMass,
    		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
    		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/mass_center_units.txt")

	TEST_FUNCTION(testNunitsCenterAndMass,
				"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
				"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/mass_center_units.txt")

    finish();
}

void TestUnitsBuilder::finish(){
    Test::finish();
}

bool TestUnitsBuilder::testUnitaryBondVectors(const char* complex_path, const char* vectors_file_path){

	System sys;
	vector<Unit*> units;
	CenterOfMass* com;
	Complex* complex = TestANMICTools::createUnitsFromFilePickingCOM(complex_path, units, com, false);

	vector<vector<double> > unitary_bond_vectors;
	TestTools::load_vector_of_vectors(unitary_bond_vectors, vectors_file_path);

	bool all_vectors_are_equal = true;
	unsigned int number_of_torsions = units.size()-1;
	for (unsigned int i = 0; i < number_of_torsions; ++i){
		all_vectors_are_equal &= Assertion::expectedVectorEqualsCalculatedWithinPrecision(
				Utils::vectorToPointer(unitary_bond_vectors[2*i+1]),
				&(units[i]->e_right->getCoordinates()[0]), 3, 1e-7);
	}

	delete complex;
	delete com;
	Utils::clearVector<Unit>(units);

	return all_vectors_are_equal;
}

bool TestUnitsBuilder::testNunitsCenterAndMass(const char* complex_path, const char* expected_data_file){
	// vector layout [N Units (1xi), Mass (1xf), Center (3xf)]
	vector<double> values;
	TestTools::load_vector(values, expected_data_file);
	int number_of_units = (int) values[0];
	double total_mass = values[1];
	Point calculated_center(values[2], values[3], values[4]);

	System sys;
	vector<Unit*> units;
	CenterOfMass* com;

	bool USE_INMA_ATOM_WEIGHTS = true;
	TestANMICTools::createUnitsFromFilePickingCOM(complex_path, units, com, USE_INMA_ATOM_WEIGHTS);

	bool all_ok = Assertion::expectedEqualsCalculated(number_of_units, units.size());

	all_ok = all_ok && Assertion::expectedEqualsCalculatedWithinPrecision(total_mass, com->getMass(), 1e-3);

	all_ok = all_ok && Assertion::expectedEqualsCalculatedWithinPrecision(calculated_center.getX(),
																	com->getX(), 1e-7);

	all_ok = all_ok && Assertion::expectedEqualsCalculatedWithinPrecision(calculated_center.getY(),
																	com->getY(), 1e-7);

	all_ok = all_ok && Assertion::expectedEqualsCalculatedWithinPrecision(calculated_center.getZ(),
																	com->getZ(), 1e-7);

	delete com;

	return all_ok;
}
