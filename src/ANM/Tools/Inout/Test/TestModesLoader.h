/////////////////////////////////////////////////////////////////////////////
/// \file TestModesLoader.h
///
/// \brief bla
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author priera
/// \date May 2, 2014
/////////////////////////////////////////////////////////////////////////////

#ifndef TESTMODESLOADER_H_
#define TESTMODESLOADER_H_

#include <iosfwd>
#include "../../../../Tools/TestTools.h"

class ModesLoader;

class TestModesLoader : public Test {
public:
	TestModesLoader(std::string name);
	virtual ~TestModesLoader();

    void run();
    void init();

  private:
	void finish();

	bool testParsingAtomInfo(const char * modesFilePath, const char * atomNamesFile,
			const char * residuesIdsFile, const char * residueNamesFile, const char * chainIdsFile, unsigned int numberOfModes);

	bool testNodeListIndexes(const char * modesFilePath, const char * dataPath, const char * goldenDataIndexes,
			unsigned int numberOfModes);

	bool testModesLoading(const char * modesFilePath, const char * dataPath,  const char * goldenDataPath, unsigned int numberOfModes);

	void printDict(const ModesLoader & loader);

};

#endif /* TESTMODESLOADER_H_ */
