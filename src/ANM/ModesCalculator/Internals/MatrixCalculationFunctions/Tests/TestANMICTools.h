/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 13/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTANMICTOOLS_H_
#define TESTANMICTOOLS_H_
#include <vector>


class Complex;
class Chain;
class Unit;
class Atom;
class CenterOfMass;

namespace TestANMICTools {
	bool pointsLookSimilar(double*a,
							double*b);

	void changeCoordsByThoseInTheFile(Chain* chain,
									const char* coordinates_file_name);

	Chain* createChain(const char* file_path,
					Complex*& r_complex);

	Chain* createUnitsFromFile(const char* file_path,
								std::vector<Unit*>& units,
								Complex*& r_complex,
								bool skip_OXT);

	Complex* createUnitsFromFilePickingCOM(const char* file_path,
									std::vector<Unit*>& units,
									CenterOfMass*& com, bool injectINMAMasses);


	void createUnitsFromFileChangingCoordinates(const char* file_path,
								const char* coordinates_file_name,
								std::vector<Unit*>& units,
								Complex*& r_complex);

	void changeAtomMasses(std::vector<Atom*>& atoms);
};

#endif /* TESTANMICTOOLS_H_ */
