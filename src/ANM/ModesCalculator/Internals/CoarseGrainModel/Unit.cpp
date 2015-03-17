#include "Unit.h"
#include <iostream>
#include <vector>
#include <string>
#include "../../../../Tools/Math/Point.h"
#include "../../../../Tools/vectorTools.h"
#include "../../../../Molecules/ForceFieldTerms/Dihedral.h"
#include "../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"

using namespace std;



///////////////////////////////////////////////////////////////
/// \remarks
/// Creates a unit with all its parameters.
///     [Unit i -1] <---e_left----r_left[atom 2, atom 3.... ] r_rigth ------e_right----> [Unit i+1]
///                              [            Unit i                 ]

/// \param atoms [In] A list with all the heavy atoms the unit is composed of.
/// \param hydrogens [In] The hydrogens bonded to the first list of atoms.
/// \param left, right [In] Left and right dihedral for the unit.
/// \param e_left, e_right [In] The axis of rotation (unit vectors).
/// \param resname [In] Name of the residue the unit is representing.
/// \param description [In] String description of that unit.
///
/// \author vgil
/// \date 1/07/2015
///////////////////////////////////////////////////////////////
Unit::Unit(vector<Atom*>& atoms,
			vector<Atom*>& hydrogens,
			Dihedral* left, Dihedral* right,
			Point* e_left, Point* e_right,
			Point* r_left, Point* r_right,
			Atom* main_atom,
			string resname,
			string description ){

	VectorTools::copy<Atom*>(this->atoms, atoms);
	VectorTools::copy<Atom*>(this->hydrogens, hydrogens);

	this->left_dihedral = left;
	this->right_dihedral = right;

	com = new CenterOfMass(CenterOfMass::compute(atoms));

	this->e_left = e_left;
	this->e_right = e_right;
	this->r_left = r_left;
	this->r_right = r_right;

	this->resname = resname;

	this->description = description;

	this->main_atom = main_atom;
}

Unit::~Unit() {
	delete com;
	delete e_left;
	delete e_right;
	delete r_left;
	delete r_right;
}

Point Unit::getCenterOfMass() const{
	return Point(com->getX(),com->getY(), com->getZ());
}

std::vector<Atom*> Unit::get_all_atoms(){
	// TODO: improve this with bookkeeping! Units must be inmutable but its atoms are not
	vector<Atom*> all;

	VectorTools::add(all, hydrogens);
	VectorTools::add(all, atoms);

	return all;
}

double Unit::getMass() const{
	return com->getMass();
}

Dihedral* Unit::getRightDihedral() {
	return right_dihedral;
}

Dihedral* Unit::getLeftDihedral() {
	return left_dihedral;
}

string Unit::toString(){
	return this->resname + " ("+this->description+")";
}
