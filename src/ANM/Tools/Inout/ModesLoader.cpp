/////////////////////////////////////////////////////////////////////////////
/// \file ModesLoader.cpp
///
/// \brief Implementations of the functions of the ModesLoader class.
/// TODO: In order to maximize encapsulation, this class can be splitted in several ones. This has not been done since currently is not necessary
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author priera
/// \date 28/04/2014
/////////////////////////////////////////////////////////////////////////////

#include "ModesLoader.h"

#include <stdlib.h>
#include <fstream>
#include <algorithm>

#include "../../ModesCalculator/AnmEigen.h"
#include "../../../Exceptions/ReadFileException.h"
#include "../../../Tools/stringTools.h"
#include "../../../Molecules/Atom.h"
#include "../../../Molecules/AtomSet/AtomSet.h"
#include "../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../Molecules/Selection/Selector.h"
#include "../../../Molecules/AtomPrinter.h"
#include "../../../Tasks/Configuration/BlockNames.h"

using namespace std;

//Contains the identifier that appears at the beginning of each line
namespace LineNames {

	const string NMWIZ_LOAD = "nmwiz_load";
	const string TITLE = "name";
	const string NAMES = "atomnames";
	const string RESNAMES = "resnames";
	const string CHIDS = "chainids";
	const string RESIDS = "resids";
	const string BETAS = "bfactors";
	const string COORDINATES = "coordinates";
	const string MODE = "mode";
	const string SEGNAMES = "segnames";

}

///////////////////////////////////////////////////////////////
/// \param fileName [In] Name of the file to be parsed
/// \param numberOfModes [In] Number of modes to be loaded. No more modes than this will be used, despite the file can contain more.
/// \param numberOfNodes [In] Number of nodes
/// \param nodeList [In] List of nodes
///
/// \author priera
/// \date 20/05/2014
///////////////////////////////////////////////////////////////
ModesLoader::ModesLoader(string fileName, unsigned int numberOfModes, unsigned int numberOfNodes,
		const std::vector<Atom*> & nodesList) :
		numberOfModes(numberOfModes), numberOfNodes(numberOfNodes), nodesList(&nodesList) {

	eigenValues = (double *)calloc(numberOfModes, sizeof(double));
	eigenVectors = (double *)calloc(numberOfModes * numberOfNodes * 3, sizeof(double));

	generateParsingDataStructures();

	load(fileName);

}

///////////////////////////////////////////////////////////////
/// \param numberOfModes [In] Number of modes to be loaded. No more modes than this will be used, despite the file can contain more.
/// \param numberOfNodes [In] Number of nodes
/// \param nodeList [In] List of nodes
///
/// This constructor is meant for testing purposes only, so it won't do the whole parsing procedure. It only will create the basic data structures.
///
/// \author priera
/// \date 20/05/2014
///////////////////////////////////////////////////////////////
ModesLoader::ModesLoader(unsigned int numberOfModes, unsigned int numberOfNodes, const std::vector<Atom*> & nodesList) :
	numberOfModes(numberOfModes), numberOfNodes(numberOfNodes), nodesList(&nodesList) {

	eigenValues = (double *)calloc(numberOfModes, sizeof(double));
	eigenVectors = (double *)calloc(numberOfModes * numberOfNodes * 3, sizeof(double));

	generateParsingDataStructures();
}

ModesLoader::~ModesLoader() {
	free(eigenValues);
	free(eigenVectors);
}

void ModesLoader::generateParsingDataStructures(){

	typesOfLine[LineNames::NMWIZ_LOAD] = NOT_RELEVANT_LINE;
	typesOfLine[LineNames::TITLE] = NOT_RELEVANT_LINE;
	typesOfLine[LineNames::BETAS] = NOT_RELEVANT_LINE;
	typesOfLine[LineNames::SEGNAMES] = NOT_RELEVANT_LINE;
	typesOfLine[LineNames::MODE] = MODES;
	typesOfLine[LineNames::NAMES] = NAMES;
	typesOfLine[LineNames::CHIDS] = CHAIN_IDS;
	typesOfLine[LineNames::RESIDS] = RESIDUE_NUMBERS;
	typesOfLine[LineNames::RESNAMES] = RESIDUE_NAMES;
	typesOfLine[LineNames::COORDINATES] = COORDINATES;

	linesHaveBeenSeen[NAMES] = false;
	linesHaveBeenSeen[CHAIN_IDS] = false;
	linesHaveBeenSeen[RESIDUE_NUMBERS] = false;
	linesHaveBeenSeen[COORDINATES] = false;
}

///////////////////////////////////////////////////////////////
/// \param fileName [In] Name of the nmd file
///
/// Does the actual parsing of the nmd file
///
/// \author priera
/// \date 21/05/2014
///////////////////////////////////////////////////////////////
void ModesLoader::load(const string & fileName){

	try {
		ifstream * stream = getStream(fileName);

		ParsingData data;

		parseLines(*stream, data);

		processModes(data);

		closeAndDeleteStream(stream);

	} catch(exception &e){
		throw ReadFileException(string("Error when processing file ") + fileName, "ModesLoader::load");
	}

}

void ModesLoader::parseLines(istream & stream, ParsingData & data){

	string currentLine;

	getline(stream, currentLine);

	while (not stream.eof()){

		parseLine(currentLine, data);

		getline(stream, currentLine);
	}

	if (not linesHaveBeenSeen[COORDINATES]){
		throw ReadFileException("There is no coordinates line in the .nmd file ", \
			"ModesLoader::parseLines");
	}

}

///////////////////////////////////////////////////////////////
/// \param data [In] Data collected during the parsing phase
///
/// Process the modes that have been stored during the parsing of the file
///
/// \author priera
/// \date 21/05/2014
///////////////////////////////////////////////////////////////
void ModesLoader::processModes(ParsingData & data){

	findNodeListIndexes(data);

	parseModesLines(data);
}

///////////////////////////////////////////////////////////////
/// \param data [In] Data collected during the parsing phase
///
/// Finds for each node collected during the parsing which is its index in the nodeList passed as argument to the constructor. This is because the
/// final vectors will be ordered according the nodeList created previously during build time, and the order may differ from what is specified in the
/// nmd file
///
/// The fact that we read information about the nodes to which apply normal modes from the nmd file introduces a redundancy: nodes
/// are specified in both the control file and the nmd file. Due this reason, this class wants as input the list of nodes that have been specified
/// in the control file, and will complain if what's in the nmd file does not match what's in the control file
///
/// \author priera
/// \date 21/05/2014
///////////////////////////////////////////////////////////////
void ModesLoader::findNodeListIndexes(ParsingData & data){

	map<unsigned int, AtomInfo> & infoMap = data.mapOfAtomInfo;
	unsigned int atomIndex;
	for (map<unsigned int, AtomInfo>::iterator it = infoMap.begin(); it != infoMap.end(); ++it) {
		atomIndex = findAtomIndexInNodeList(it->second);
		it->second.nodeListIndex = atomIndex;
	}

	//TODO: check whether an atom in the node list has not been processed by this function
}

unsigned int ModesLoader::findAtomIndexInNodeList(const AtomInfo & referenceAtom){

	unsigned int ret;
	unsigned int nodesListSize = nodesList->size();

	AtomPrinter printer(referenceAtom.atomName, referenceAtom.residueId, referenceAtom.residueName, referenceAtom.chainId);
	const string & referenceIdentifier = printer.getPDBIdentifier(false);
	for (ret = 0; ret < nodesListSize; ret++){

		//TODO: make tests for getPDBIdentifier
		if ((*nodesList)[ret]->getPrinter()->getPDBIdentifier(false) == referenceIdentifier){
			break;
		}
	}

	if (ret == nodesListSize){
		throw runtime_error("Atom with identifier " + referenceIdentifier + " has not been found in the list of atoms selected to apply ANM");
	}

	return ret;

}

void ModesLoader::parseLine(const string & line, ParsingData & data){

	vector<string> lineElements;
	StringTools::split(lineElements, line);

	LineType type = typesOfLine.at(lineElements[0]);

	if (not validLineType(type)){
		throw runtime_error("Line identifier not valid"); // In case the same line has been seen two times
	}

	switch(type){
		case NAMES:
			collectAtomNames(lineElements, data);
			break;
		case CHAIN_IDS:
			collectChainIds(lineElements, data);
			break;
		case RESIDUE_NUMBERS:
			collectResiduesIds(lineElements, data);
			break;
		case RESIDUE_NAMES:
			collectResidueNames(lineElements, data);
			break;
		case COORDINATES:
			collectCoordinates(lineElements, data);
			break;
		case MODES:
			data.modesLines.push_back(lineElements);
			break;
		case NOT_RELEVANT_LINE:
			break;
		default:
			throw runtime_error("Not recognized line identifier");
	}
}

void ModesLoader::collectAtomNames(const vector<string> & lineElements, ParsingData & data){
	static const string ALPHACARBON = "CA";
	static const string STANDARIZED_PELE_ALPHACARBON = " CA "; //For PELE atom names have four characters

	string atomName;
	unsigned int numberOfElements = lineElements.size();
	for (unsigned int i = 1; i < numberOfElements; i++){
		atomName = lineElements[i];
		if (lineElements[i] == ALPHACARBON){
			atomName = STANDARIZED_PELE_ALPHACARBON;
		}
		data.mapOfAtomInfo[i - 1].atomName = atomName;
	}

}

void ModesLoader::collectChainIds(const vector<string> & lineElements, ParsingData & data){

	unsigned int numberOfElements = lineElements.size();
	for (unsigned int i = 1; i < numberOfElements; i++){
		data.mapOfAtomInfo[i - 1].chainId = lineElements[i][0];
	}

}

void ModesLoader::collectResiduesIds(const vector<string> & lineElements, ParsingData & data){

	unsigned int numberOfElements = lineElements.size();
	for (unsigned int i = 1; i < numberOfElements; i++){
		data.mapOfAtomInfo[i - 1].residueId = StringTools::fromString<unsigned int>(lineElements[i]);
	}

}

void ModesLoader::collectResidueNames(const vector<string> & lineElements, ParsingData & data){
	unsigned int numberOfElements = lineElements.size();
	for (unsigned int i = 1; i < numberOfElements; i++){
		data.mapOfAtomInfo[i - 1].residueName = lineElements[i];
	}
}

///////////////////////////////////////////////////////////////
/// \param lineElements [In] Vector containing the coordinates of each  node
/// \param data [In Out] Container for the information which is collected during the parsing
///
/// Not only collects the coordinates. According to the nmd specification the number of atoms is computed
/// from the information which is in this line, so this function also does so.
///
/// \author priera
/// \date 21/05/2014
///////////////////////////////////////////////////////////////
void ModesLoader::collectCoordinates(const vector<string> & lineElements, ParsingData & data){

	unsigned int numberOfElementsSoFar = 0;
	unsigned int numberOfElements = lineElements.size();
	for (unsigned int i = 1; i < numberOfElements; i++){
		for (int j = 0; j < 3; j++){
			stringstream stream(lineElements[i]);
			stream >> data.mapOfAtomInfo[numberOfElementsSoFar].coords[j];
		}

		numberOfElementsSoFar = (i - 1) / 3;
	}

	data.numberOfAtoms = numberOfElementsSoFar + 1; //numberOfElementsSoFar is a 0-based index, while numberOfAtoms is 1-based
}

void ModesLoader::parseModesLines(const ParsingData & data){

	if (data.modesLines.size() == 0){
		throw runtime_error("Input .nmd file does not have information about normal modes");
	}

	unsigned int numberOfElements = data.modesLines.size();
	if (numberOfModes < numberOfElements ) numberOfElements = numberOfModes;

	for (unsigned int i = 0; i < numberOfElements; i++){
		parseModeLine(data.modesLines[i], i, data);
	}

}

///////////////////////////////////////////////////////////////
/// \param line [In] Mode line from the nmd file
/// \param normalModeIndex [In] Index of this normal mode (0 - based)
/// \param data [In] Data collected during parsing
///
/// Each modes line contains both the eigenvalue and the eigenvector components. This function stores both. Since the components of the eigenvector
/// may be in a different order than what's in the nodeList, the previously found indexes are used to save them in the appropriate positions.
///
/// \author priera
/// \date 21/05/2014
///////////////////////////////////////////////////////////////
void ModesLoader::parseModeLine(const vector<string> & line, unsigned int normalModeIndex, const ParsingData & data){
	unsigned int firstComponentIndex;
	checkLineSizeAndStoreEigenValue(line, normalModeIndex, data, firstComponentIndex);

	parseEigenVectorComponents(line, data, normalModeIndex, firstComponentIndex);
}

void ModesLoader::parseEigenVectorComponents(const vector<string> & line, const ParsingData & data, unsigned int normalModeIndex,
		unsigned int firstComponentIndex){

	unsigned int elementsInLine = line.size();
	unsigned int vectorOffset = (3 * data.numberOfAtoms) * normalModeIndex;

	const map<unsigned int, AtomInfo> & mapInfo = data.mapOfAtomInfo;

	for (unsigned int i = firstComponentIndex; i < elementsInLine; i++){
		double component = StringTools::toDouble(line[i]);
		unsigned int currentAtomIndex = (i - firstComponentIndex) / 3;
		unsigned int currentComponent = (i - firstComponentIndex) % 3;
		unsigned int baseIndexInRow = mapInfo.at(currentAtomIndex).nodeListIndex * 3;
		unsigned int rowOffset = baseIndexInRow + currentComponent;

		eigenVectors[vectorOffset + rowOffset]  = component;
	}
}

void ModesLoader::checkLineSizeAndStoreEigenValue(const vector<string> & line, unsigned int & normalModeIndex, const ParsingData & data,
		unsigned int & firstComponentIndex){

	const unsigned int firstIndexWhenEigenValueIsPresent = 2;
	const unsigned int firstIndexWhenEigenValueAndIndexArePresent = 3;

	unsigned int lineShouldHaveAtLeast = data.numberOfAtoms * 3 + 2;
	unsigned int lineShouldHaveAtMost = data.numberOfAtoms * 3 + 3;
	unsigned int sizeOfLine = line.size();
	double scalingFactor;
	double eigenValue;

	if (sizeOfLine < lineShouldHaveAtLeast or sizeOfLine > lineShouldHaveAtMost){
		throw runtime_error("Wrongly formated normal mode");
	} else	if (sizeOfLine == lineShouldHaveAtLeast){
		checkLineHasEigenValue(line, scalingFactor);
		firstComponentIndex = firstIndexWhenEigenValueIsPresent;
	} else {
		normalModeIndex = StringTools::fromString<unsigned int>(line[1]) - 1;
		scalingFactor = StringTools::fromString<double>(line[2]);
		firstComponentIndex = firstIndexWhenEigenValueAndIndexArePresent;
	}

	/*According to the file format documentation (http://prody.csb.pitt.edu/manual/reference/dynamics/nmdfile.html), what is provided is a scaling factor
	which is the square-root of the eigenvalue*/
	eigenValue = scalingFactor * scalingFactor;
	eigenValues[normalModeIndex] = eigenValue;

}

void ModesLoader::checkLineHasEigenValue(const vector<string> & line, double & scalingFactor){

	if (line[1].find(".") == string::npos){
		throw runtime_error("Normal mode index found in .nmd file, but eigenvalue is required too");
	}

	scalingFactor = StringTools::fromString<double>(line[1]);
}

bool ModesLoader::validLineType(const ModesLoader::LineType & typeOfLine) {
	if (typeOfLine == MODES or typeOfLine == NOT_RELEVANT_LINE) {
		// We don't mind if some line not used by PELE is repeated. Also, modes lines can be repeated
		return true;
	}

	bool lineHasAlreadyBeenSeen = linesHaveBeenSeen[typeOfLine];
	linesHaveBeenSeen[typeOfLine] = true;
	return not lineHasAlreadyBeenSeen;
}

void ModesLoader::loadIntoEigen(AnmEigen * eigen){
	eigen->initialize(eigenValues, eigenVectors, numberOfModes, numberOfNodes);
}

std::ifstream * ModesLoader::getStream(string fileName){
	ifstream * stream = new ifstream(fileName.c_str());

	if (not *stream) {
		throw ReadFileException(string("Error when opening file " ) + fileName, "ModesLoader::getStream");
	}

	return stream;

}

void ModesLoader::closeAndDeleteStream(std::ifstream * stream){
	stream->close();
	delete stream;
}

