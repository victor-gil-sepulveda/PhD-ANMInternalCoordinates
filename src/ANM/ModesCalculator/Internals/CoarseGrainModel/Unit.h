/////////////////////////////////////////////////////////////////////////////
/// \brief Description of each of the units in the coarse grain model used for
/// the ANM algorithm in internal coordinates.
///
/// \author vgil
/// \date 06/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef UNIT_H_
#define UNIT_H_

///////////////////////////////////////////////////////////
/// \brief Models a coarse grain element composed of many atoms. Units may be generalized so that they can
/// represent any coarse grain model. This Unit implementation is too specific, as it is used only in the
/// IC ANM calculations.
///////////////////////////////////////////////////////////

#include <vector>
#include <string>

class Point;
class CenterOfMass;
class Dihedral;
class Atom;

class Unit {
	public:
		Unit(	std::vector<Atom*>& atoms,
				std::vector<Atom*>& hydrogens,
				Dihedral* right,
				Atom* bond_atom_left,
				Atom* bond_atom_right,
				Atom* main_atom,
				std::string resname,
				std::string description = "");

		virtual ~Unit();

		void update();

		std::vector<Atom*> getAllAtoms();

		double getMass() const;

		std::string toString();

		// Singular atoms to define the torsion around the valence bond
		Atom *right_torsion_bond_atom, *left_torsion_bond_atom;
		Point *e_right, *r_right;

		// Contents of this unit
		std::vector<Atom*> atoms, hydrogens;

		// Other useful data
		Dihedral* right_dihedral;
		CenterOfMass* com;

		// Characterization of the unit
		Atom* main_atom;
		std::string resname;
		std::string description;
};

#endif /* UNIT_H_ */
