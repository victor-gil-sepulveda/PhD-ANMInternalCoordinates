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

#include "TestModesLoader.h"

#include <string>
#include <map>
#include <iosfwd>
#include <json/value.h>

#include "../ModesLoader.h"
#include "../../../NodeListGenerator/NodeListGeneratorBuilder.h"
#include "../../../NodeListGenerator/NodeListGenerator.h"
#include "../../../Algorithms/AnmAlgorithmTypes.h"
#include "../../../ModesCalculator/AnmEigen.h"
#include "../../../../Tools/TestTools.h"
#include "../../../../Tools/Assertion.h"
#include "../../../../Molecules/ComplexBuilder.h"
#include "../../../../Molecules/Complex.h"
#include "../../../../Molecules/Atom.h"
#include "../../../../Molecules/AtomPrinter.h"

using namespace std;

#include <iostream>

TestModesLoader::TestModesLoader(string name) {
    test_name = name;
}

TestModesLoader::~TestModesLoader() {
}

void TestModesLoader::init(){
}

void TestModesLoader::run()
{
    Test::run();

    TEST_FUNCTION(testParsingAtomInfo,
       		"src/ANM/Tools/Inout/Test/Data/inputfile",
       		"src/ANM/Tools/Inout/Test/Data/atomNames",
       		"src/ANM/Tools/Inout/Test/Data/residuesIds",
       		"src/ANM/Tools/Inout/Test/Data/residueNames",
       		"src/ANM/Tools/Inout/Test/Data/chainIds",
       		6)

    TEST_FUNCTION(testNodeListIndexes,
    		"src/ANM/Tools/Inout/Test/Data/plop_isr_out.nmd",
    		"src/Molecules/Test/Data/plop_isr_out.pdb",
    		"src/ANM/Tools/Inout/Test/Data/nodeListIndexes",
    		6)

    TEST_FUNCTION(testModesLoading,
    		"src/ANM/Tools/Inout/Test/Data/plop_isr_out.nmd",
    		"src/Molecules/Test/Data/plop_isr_out.pdb",
    		"src/ANM/Tools/Inout/Test/Data/eigenVectorsGoldenData",
    		6)

    finish();
}

bool TestModesLoader::testParsingAtomInfo(const char * modesFilePath, const char * atomNamesFile,
			const char * residuesIdsFile, const char * residueNamesFile, const char * chainIdsFile, unsigned int numberOfModes){
	//Set up
	unsigned int notUsedValue = 0;
	vector<Atom *> * nodeListPointer = TestTools::NOT_USED<vector<Atom*> >();
	ModesLoader loader(numberOfModes, notUsedValue, *nodeListPointer);

	//Code to test
	ModesLoader::ParsingData data;
	ifstream * stream = loader.getStream(string(modesFilePath));
	loader.parseLines(*stream, data);

	//Assert
	vector<string> obtainedAtomNames;
	vector<unsigned int> obtainedResiduesIds;
	vector<string> obtainedResidueNames;
	vector<char> obtainedChainIds;

	const std::map<unsigned int, ModesLoader::AtomInfo> & mapOfAtomInfo = data.mapOfAtomInfo;
	for (std::map<unsigned int, ModesLoader::AtomInfo>::const_iterator it = mapOfAtomInfo.begin(); it != mapOfAtomInfo.end(); ++it){
		const ModesLoader::AtomInfo & info = it->second;

		obtainedAtomNames.push_back(info.atomName);
		obtainedResiduesIds.push_back(info.residueId);
		obtainedChainIds.push_back(info.chainId);
		obtainedResidueNames.push_back(info.residueName);
	}

	vector<string> expectedAtomNames;
	vector<unsigned int> expectedResiduesIds;
	vector<string> expectedResidueNames;
	vector<char> expectedChainIds;

	TestTools::load_vector(expectedAtomNames, atomNamesFile);
	TestTools::load_vector(expectedResiduesIds, residuesIdsFile);
	TestTools::load_vector(expectedResidueNames, residueNamesFile);
	TestTools::load_vector(expectedChainIds, chainIdsFile);

	bool allOk = Assertion::expectedVectorEqualsCalculated(expectedAtomNames, obtainedAtomNames);
	allOk &= Assertion::expectedVectorEqualsCalculated(expectedResiduesIds, obtainedResiduesIds);
	allOk &= Assertion::expectedVectorEqualsCalculated(expectedResidueNames, obtainedResidueNames);
	allOk &= Assertion::expectedVectorEqualsCalculated(expectedChainIds, obtainedChainIds);

	//Clean-up
	loader.closeAndDeleteStream(stream);

	return allOk;
}

bool TestModesLoader::testNodeListIndexes(const char * modesFilePath, const char * dataPath,
		const char * expectedIndexes, unsigned int numberOfModes){

	//Set up
	ComplexBuilder complexBuilder(OPLS2005);
	Complex * complex = complexBuilder.build(dataPath);

	NodeListGeneratorBuilder nodeListGeneratorBuilder;
	NodeListGenerator * nodeListGenerator = nodeListGeneratorBuilder.createFromConfiguration(Json::nullValue, ANM_CARTESIAN_ALPHACARBONS);

	std::vector<Atom*> nodeList;
	nodeListGenerator->generateNodeList(complex, nodeList);

	ModesLoader loader(numberOfModes, nodeList.size(), nodeList);
	ifstream * stream = loader.getStream(string(modesFilePath));
	ModesLoader::ParsingData data;
	loader.parseLines(*stream, data);

	//Code to test
	loader.findNodeListIndexes(data);

	//Assert
	vector<unsigned int> expectedDataVector;
	TestTools::loadVectorOfIntegerValues(expectedDataVector, expectedIndexes);

	vector<unsigned int> obtainedVector;
	map<unsigned int, ModesLoader::AtomInfo> & infoMap = data.mapOfAtomInfo;
	for (map<unsigned int, ModesLoader::AtomInfo>::iterator it = infoMap.begin(); it != infoMap.end(); ++it) {
		obtainedVector.push_back(it->second.nodeListIndex);
	}

	bool ok = Assertion::expectedVectorEqualsCalculated(expectedDataVector, obtainedVector);

	//Clean-up
	loader.closeAndDeleteStream(stream);

	delete nodeListGenerator;
	delete complex;

	return ok;
}

bool TestModesLoader::testModesLoading(const char * modesFilePath, const char * dataPath, const char * goldenDataPath, unsigned int numberOfModes){

	//Set up
	ComplexBuilder complexBuilder(OPLS2005);
	Complex * complex = complexBuilder.build(dataPath);

	NodeListGeneratorBuilder nodeListGeneratorBuilder;
	NodeListGenerator * nodeListGenerator = nodeListGeneratorBuilder.createFromConfiguration(Json::nullValue, ANM_CARTESIAN_ALPHACARBONS);

	std::vector<Atom*> nodeList;
	nodeListGenerator->generateNodeList(complex, nodeList);
	AnmEigen * eigen = new AnmEigen();

	//Code to test
	ModesLoader loader(string(modesFilePath), numberOfModes, nodeList.size(), nodeList);
	loader.loadIntoEigen(eigen);

	//Assert
	double expectedEigenValuesArray[] = {7.48121e-06,0.000294294,0.00050402,0.0345915,0.193013,1.1663352009};
	vector<double> expectedEigenValues(expectedEigenValuesArray, expectedEigenValuesArray + sizeof(expectedEigenValuesArray) / sizeof(double) );

	vector<double> expectedComponents;
	TestTools::load_vector(expectedComponents, goldenDataPath);

	vector<double> obtainedComponents;
	vector<double> obtainedValues;
	for (unsigned int i = 0; i < numberOfModes; i++){
		obtainedValues.push_back(eigen->getEigenValueOfMode(i));

		const vector<double> & currentMode = eigen->getEigenVectorOfMode(i);
		unsigned int numberOfComponents = currentMode.size();
		for (unsigned int j = 0; j < numberOfComponents; j++){
			obtainedComponents.push_back(currentMode[j]);
		}
	}

	bool vectorsAreEqual = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expectedComponents, obtainedComponents, 1.e-6);
	bool valuesAreEqual = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expectedEigenValues, obtainedValues, 1.e-6);

	//Clean up
	delete eigen;
	delete nodeListGenerator;
	delete complex;

	return vectorsAreEqual and valuesAreEqual;
}

void TestModesLoader::finish()
{
    Test::finish();
}

