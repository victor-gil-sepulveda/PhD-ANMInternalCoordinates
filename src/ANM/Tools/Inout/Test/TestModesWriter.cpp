/////////////////////////////////////////////////////////////////////////////
/// \file TestModesLoader.cpp
///
/// \brief bla
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author priera
/// \date May 2, 2014
/////////////////////////////////////////////////////////////////////////////

#include "TestModesWriter.h"

#include <string>
#include <vector>
#include <iostream>
#include "../../../ModesCalculator/Internals/MatrixCalculationFunctions/Tests/TestANMICTools.h"
#include "../../../../Molecules/Atom.h"
#include "../ModesWriter.h"
#include "../../../../System/System.h"
#include "../../../../Tools/Assertion.h"
#include "../../../ModesCalculator/Internals/MatrixCalculationFunctions/ANMICMath.h"
#include "../../../ModesCalculator/AnmEigen.h"
using namespace std;

TestModesWriter::TestModesWriter(string name) {
    test_name = name;
}

TestModesWriter::~TestModesWriter(){
}

void TestModesWriter::init(){
}

void TestModesWriter::run()
{
    Test::run();

    TEST_FUNCTION(testCAFiltering,
        		"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala_pro_ala/ala_pro_ala.pdb",
        		"src/ANM/Tools/Inout/Test/writer_data/ala_pro_ala_ca_indices.txt")

    TEST_FUNCTION(testCAFiltering,
				"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb",
				"src/ANM/Tools/Inout/Test/writer_data/ala5_indices.txt")

	TEST_FUNCTION(testEigenFiltering,
			"src/ANM/Tools/Inout/Test/writer_data/ala_pro_ala_ca_indices.txt",
			"src/ANM/Tools/Inout/Test/writer_data/apa_eigenvectors.txt",
			"src/ANM/Tools/Inout/Test/writer_data/apa_filtered_eigenvectors.txt")

	TEST_FUNCTION(testEigenFiltering,
			"src/ANM/Tools/Inout/Test/writer_data/ala5_indices.txt",
			"src/ANM/Tools/Inout/Test/writer_data/ala5_eigenvectors.txt",
			"src/ANM/Tools/Inout/Test/writer_data/ala5_filtered_eigenvectors.txt")

    finish();
}

bool TestModesWriter::testCAFiltering(const char* prot_path,
		const char* expected_indices_path){

	// Load units from a file
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model"<<endl;
	TestANMICTools::createUnitsFromFile(prot_path, units, complex, false);

	// Get the indices
	vector<int> indices;
	vector<Atom*> ca_atoms;
	ModesWriter::getCAsFromUnits(units,indices,ca_atoms);

	// Load the golden indices
	vector<int> expected_indices;
	TestTools::loadVectorOfIntegerValues(expected_indices, expected_indices_path);

	// Check the indices
	return Assertion::expectedVectorEqualsCalculated(expected_indices, indices);
}

bool TestModesWriter::testEigenFiltering(const char* indices_to_keep,
				const char* eigenvectors_file,
				const char* expected_eigenvectors_file){

	vector<int> indices;
	TestTools::loadVectorOfIntegerValues(indices, indices_to_keep);

	vector<vector<double> > eigenvectors, expected_eigenvectors;
	TestTools::load_vector_of_vectors(eigenvectors, eigenvectors_file);
	vector<double> eigenvalues(0, eigenvectors.size());
	TestTools::load_vector_of_vectors(expected_eigenvectors, expected_eigenvectors_file);

	AnmEigen eigen, filtered;
	bool isCartesian = true;
	eigen.initialize(eigenvalues, eigenvectors, isCartesian);

	ModesWriter::filterEigenvectors(&eigen, &filtered, indices);

	//ANMICMath::printMatrix(filtered.vectors);

	bool ok = true;
	for (unsigned i = 0; i < expected_eigenvectors.size(); ++i){
		ok = ok && Assertion::expectedVectorEqualsCalculatedWithTolerance(expected_eigenvectors[i],
				filtered.vectors[i],
				1e-12);
	}
	return ok;
}

void TestModesWriter::finish()
{
    Test::finish();
}

