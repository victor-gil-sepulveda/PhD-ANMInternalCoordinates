/////////////////////////////////////////////////////////////////////////////
/// \file ModesLoader.h
///
/// \brief Definition of the ModesLoader class
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author priera
/// \date 28/04/2014
/////////////////////////////////////////////////////////////////////////////

#ifndef MODESLOADER_H_
#define MODESLOADER_H_

#include <string>
#include <vector>
#include <map>

class Atom;
class AtomSet;
class AnmEigen;

/////////////////////////////////////////////////////////////////////////////
/// \brief Given a file containing the eigenvalues and eigenvectors for a given system, this class parses it, and stores the data in a format able to
/// fill an AnmEigen object with it
///
/// The input file has to be in nmd format (http://prody.csb.pitt.edu/manual/reference/dynamics/nmdfile.html)
///
/// \author priera
/// \date 20/05/2012
/////////////////////////////////////////////////////////////////////////////
class ModesLoader {
public:
	ModesLoader(std::string fileName, unsigned int numberOfModes, unsigned int numberOfNodes, const std::vector<Atom*> & nodesList);
	virtual ~ModesLoader();

	void loadIntoEigen(AnmEigen * eigen);

private:

	//Constructor only for test-purposes
	ModesLoader(unsigned int numberOfModes, unsigned int numberOfNodes, const std::vector<Atom*> & nodesList);

	typedef enum {
		NOT_RELEVANT_LINE,
		COORDINATES,
		MODES,
		NAMES,
		CHAIN_IDS,
		RESIDUE_NUMBERS,
		RESIDUE_NAMES
	} LineType;

	typedef struct {
		std::string atomName;
		unsigned int residueId;
		std::string residueName;
		char chainId;
		double coords[3];
		unsigned int nodeListIndex;
	} AtomInfo;

	typedef struct {
		std::map<unsigned int, AtomInfo> mapOfAtomInfo;
		std::vector<std::vector<std::string> > modesLines;
		unsigned int numberOfAtoms;
	} ParsingData;

	void generateParsingDataStructures();
	void load(const std::string & fileName);
	void parseLines(std::istream & stream, ParsingData & data);
	void parseLine(const std::string & line, ParsingData & data);
	bool validLineType(const ModesLoader::LineType & typeOfLine);
	void collectAtomNames(const std::vector<std::string> & lineElements, ParsingData & data);
	void collectChainIds(const std::vector<std::string> & lineElements, ParsingData & data);
	void collectResiduesIds(const std::vector<std::string> & lineElements, ParsingData & data);
	void collectResidueNames(const std::vector<std::string> & lineElements, ParsingData & data);
	void collectCoordinates(const std::vector<std::string> & lineElements, ParsingData & data);
	unsigned int findAtomIndexInNodeList(const AtomInfo & referenceAtom);
	void parseModeLine(const std::vector<std::string> & line, unsigned int normalModeIndex, const ParsingData & data);
	void processModes(ParsingData & data);
	void findNodeListIndexes(ParsingData & data);
	void parseModesLines(const ParsingData & data);
	void checkLineSizeAndStoreEigenValue(const std::vector<std::string> & line, unsigned int & normalModeIndex, const ParsingData & data,
			unsigned int & firstComponentIndex);
	void checkLineHasEigenValue(const std::vector<std::string> & line, double & eigenValue);
	void parseEigenVectorComponents(const std::vector<std::string> & line, const ParsingData & data, unsigned int normalModeIndex,
			unsigned int firstComponentIndex);
	void loadEigenvectorsAndEigenvalues(std::istream * file);

	std::ifstream * getStream(std::string fileName);
	void closeAndDeleteStream(std::ifstream * stream);

	unsigned int numberOfModes;
	unsigned int numberOfNodes;
	const std::vector<Atom*> * nodesList;

	std::map<std::string, LineType> typesOfLine; //Map used to decide which function has to be called to parse a line
	std::map<LineType, bool> linesHaveBeenSeen; //Some lines can appear only once, so this maps allows to see that

	double * eigenVectors;
	double * eigenValues;

	friend class TestModesLoader;

};

#endif /* MODESLOADER_H_ */
