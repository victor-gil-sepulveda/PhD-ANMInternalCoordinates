/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author vgil
/// \date 15/12/2014
/////////////////////////////////////////////////////////////////////////////

#include "UnitTools.h"
#include <vector>
#include "../../../../Molecules/Atom.h"
#include "../CoarseGrainModel/Unit.h"
#include "../../../../Tools/vectorTools.h"
#include <iostream>
#include <iomanip>
#include "../../../../Tools/Math/Point.h"
#include "../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"

using namespace std;

void UnitTools::getAllAtomsFromUnits(vector<Unit*>& units, vector<Atom*>& atoms, bool onlyHeavy){

	for (unsigned int i = 0; i < units.size(); ++i){
		VectorTools::add(atoms, units[i]->atoms);
		if (!onlyHeavy){
			VectorTools::add(atoms, units[i]->hydrogens);
		}
	}
}

void UnitTools::getAllAtomsFromUnitRange(vector<Unit*>& units, vector<Atom*>& atoms, int start, int end, bool onlyHeavy){

	for (unsigned int i = start; i <= end; ++i){
		VectorTools::add(atoms, units[i]->atoms);
		if (!onlyHeavy){
			VectorTools::add(atoms, units[i]->hydrogens);
		}
	}
}

int UnitTools::getNumberOfAtomsOfUnitRange(vector<Unit*>& units, int start, int end, bool onlyHeavy){
	unsigned int num_atoms = 0;
	for (unsigned int i = start; i <= end; ++i){
		num_atoms += units[i]->atoms.size();
		if (!onlyHeavy){
			num_atoms += units[i]->hydrogens.size();
		}
	}
	return num_atoms;
}


void UnitTools::printUnits(std::vector<Unit*>& units){
	for (unsigned int i = 0; i < units.size(); ++i){
		Unit* unit = units[i];
		Point* com = dynamic_cast<Point*>(unit->com);
		//Point e_left = *(unit->e_left);
		Point e_right = *(unit->e_right);
		//Point r_left = *(unit->r_left);
		Point r_right = *(unit->r_right);
		cout << setprecision(16)
				<< "Unit resname:"<<unit->resname
				<<" description:"<<unit->description
				<<" com: "<< com->toString()
				//<< "e_left:"<<e_left.toString()
				<< "e_right"<<e_right.toString()
				//<< "r_left"<< r_left.toString()
				<< "r_right"<<r_right.toString()
				<< "mass: "<<unit->getMass()
				<< endl;
	}
}

