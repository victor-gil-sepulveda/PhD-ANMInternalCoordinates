/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 06/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef UNITSBUILDER_H_
#define UNITSBUILDER_H_

class Dihedral;
class Unit;
class Link;
class Chain;
class Atom;
class Topology;
class CenterOfMass;

#include <vector>
#include <iosfwd>

class UnitsBuilder {
	public:
		UnitsBuilder(Chain* from_this_chain,
				bool center_coordinates_at_com = true,
				bool skip_OXT = false);

		virtual ~UnitsBuilder();

		virtual void build(std::vector<Unit*>& built_units);
		void create_CG_model(std::vector<Unit*>& units, Chain* chain);
		void merge_prolines(std::vector<Unit*>& in, std::vector<Unit*>& out);

		std::vector<Unit*> buildUnitsFromLink(Link* link);
		std::vector<Unit*> buildFirstLinkUnits(Link* link);
		std::vector<Unit*> buildMiddleLinkUnits(Link* link);
		std::vector<Unit*> buildLastLinkUnits(Link* link);

		Unit* build_N_Ca_Sidechain_Unit(Link* link);
		Unit* build_Ca_Sidechain_Unit(Link* link);
		Unit* build_C_O_O_Unit(Link* link);
		Unit* build_N_C_O_Unit(Link* link, Link* next_link);

		void get_all_hydrogens_from_N_Ca_Sidechain_Unit(Link* link, std::vector<Atom*>& atoms );
		void get_all_hydrogens_from_Ca_Sidechain_Unit(Link* link, std::vector<Atom*>& atoms );
		void get_all_hydrogens_from_C_O_O_Unit(Link* link, std::vector<Atom*>& atoms );
		void get_all_hydrogens_from_N_C_O_Unit(Link* link, Link* next_link, std::vector<Atom*>& atoms );
		CenterOfMass* centerChainCoordsAtCOM(bool skip_OXT);
		Unit* merge_units(Unit* first_unit, Unit* second_unit);

		CenterOfMass* getCOM();

	protected:
		Chain* chain;
		Topology const * chain_topology;
		CenterOfMass* global_com;

		Dihedral* get_psi_dihedral( Atom*& CA, Atom*& C, Link* link, bool is_last);
		Dihedral* get_left_phi_dihedral(  Atom*& N, Atom*& CA, Link* link);
		Dihedral* get_right_phi_dihedral( Atom*& N, Atom*& CA_after, Link* link);

		static void getAtomsFromLink(std::vector<Atom* >& atoms, std::vector<std::string>& atom_names, Link* link);
		static void getAtomsAndSidechainFromLink(std::vector<Atom* >& atoms, std::vector<std::string>& atom_names, Link* link);
		static void filterHydrogens(std::vector<Atom*>& atoms, std::vector<Atom* >& filtered_atoms);
		static void filterOXT(std::vector<Atom*>& atoms, std::vector<Atom* >& filtered_atoms);
};
#endif /* UNITSBUILDER_H_ */
