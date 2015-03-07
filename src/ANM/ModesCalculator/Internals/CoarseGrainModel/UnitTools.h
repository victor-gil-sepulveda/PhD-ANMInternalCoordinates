/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author vgil
/// \date 15/12/2014
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef UNITTOOLS_H_
#define UNITTOOLS_H_

#include <vector>
class Unit;
class Atom;

namespace UnitTools {
	void getAllAtomsFromUnits(std::vector<Unit*>& units, std::vector<Atom*>& atoms, bool onlyHeavy);
	void getAllAtomsFromUnitRange(std::vector<Unit*>& units, std::vector<Atom*>& atoms, int i, int j, bool onlyHeavy);
	int  getNumberOfAtomsOfUnitRange(std::vector<Unit*>& units, int start, int end, bool onlyHeavy);
	void printUnits(std::vector<Unit*>& units);
}
#endif /* UNITTOOLS_H_ */
