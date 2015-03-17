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
				Dihedral* left, Dihedral* right,
				Point* e_left, Point* e_right,
				Point* r_left, Point* r_right,
				Atom* main_atom,
				std::string resname,
				std::string description = "");

		virtual ~Unit();

		Point getCenterOfMass() const;

		std::vector<Atom*> get_all_atoms();

		double getMass() const;

		Dihedral* getRightDihedral();
		Dihedral* getLeftDihedral();

		std::string toString();

		Point* e_left;
		Point* e_right;
		Point* r_left;
		Point* r_right;

		std::vector<Atom*> atoms;
		std::vector<Atom*> hydrogens;
		std::string resname;
		std::string description;

		Atom* main_atom;

	private:
		CenterOfMass* com;
		Dihedral* left_dihedral;
		Dihedral* right_dihedral;

};

#endif /* UNIT_H_ */
