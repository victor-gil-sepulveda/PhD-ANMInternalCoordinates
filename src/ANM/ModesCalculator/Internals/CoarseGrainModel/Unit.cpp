#include "Unit.h"
#include <iostream>
#include <vector>
#include <string>
#include "../../../../Tools/Math/Point.h"
#include "../../../../Tools/vectorTools.h"
#include "../../../../Molecules/ForceFieldTerms/Dihedral.h"
#include "../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"
#include "../../../../Tools/Utils.h"
#include "../../../../Molecules/Atom.h"

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
Unit::Unit(std::vector<Atom*>& atoms,
			std::vector<Atom*>& hydrogens,
			Dihedral* right,
			Atom* left_torsion_bond_atom,
			Atom* right_torsion_bond_atom,
			Atom* main_atom,
			std::string resname,
			std::string description ){

	VectorTools::copy<Atom*>(this->atoms, atoms);
	VectorTools::copy<Atom*>(this->hydrogens, hydrogens);

	this->right_dihedral = right;

	this->left_torsion_bond_atom = left_torsion_bond_atom;
	this->right_torsion_bond_atom = right_torsion_bond_atom;

	this->resname = resname;
	this->description = description;
	this->main_atom = main_atom;

	com = NULL;
	e_right = r_right = NULL;

	update();
}

Unit::~Unit() {
	delete com;
	delete e_right;
	delete r_right;
}

void Unit::update(){
	// compute center of mass
	Utils::deleteAndSetPointerToNull(com);
	com = new CenterOfMass(CenterOfMass::compute(atoms));

	// Calculate the axis only if the bond exists
	if (left_torsion_bond_atom!=NULL && right_torsion_bond_atom!=NULL){
		// compute ra
		Utils::deleteAndSetPointerToNull(r_right);
		r_right = new Point(left_torsion_bond_atom->toPoint());

		// compute ea
		Utils::deleteAndSetPointerToNull(e_right);
		e_right = (new Point(Point::subtract(right_torsion_bond_atom->toPoint(),left_torsion_bond_atom->toPoint())))->normalize();
	}
}


std::vector<Atom*> Unit::getAllAtoms(){
	// TODO: improve this with bookkeeping! Units must be inmutable but its atoms are not
	vector<Atom*> all;

	VectorTools::add(all, hydrogens);
	VectorTools::add(all, atoms);

	return all;
}

double Unit::getMass() const{
	return com->getMass();
}

string Unit::toString(){
	return this->resname + " ("+this->description+")";
}
