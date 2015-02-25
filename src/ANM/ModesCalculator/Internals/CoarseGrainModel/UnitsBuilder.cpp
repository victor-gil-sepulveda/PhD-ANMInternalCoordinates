#include "UnitsBuilder.h"

#include "Unit.h"
#include "../../../../Molecules/Atom.h"
#include "../../../../Tools/stringTools.h"
#include "../../../../Tools/vectorTools.h"
#include "../../../../Molecules/Topology/Topology.h"
#include "../../../../Molecules/Selection/Selector.h"
#include "../../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../../Molecules/Selection/SelectionStringsBuilder.h"
#include "../../../../Molecules/AtomNames.h"
#include "../../../../Molecules/AtomSetsTree/Links/LinkNames.h"
#include "../../../../Tools/Math/MathTools.h"
#include "../../../../Tools/Utils.h"
#include "../../../../Molecules/AtomSetsTree/Chains/Chain.h"
#include "../../../../Molecules/AtomSetsTree/Links/Link.h"
#include "../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"
#include "../../../../Tools/Math/Point.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "UnitTools.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \brief
/// The following functions can be used to generate the IC Heavy Atoms coarse grain model. Function names
/// are quite self-explanatory. Based on Jose Ramón Lopez Blanco Doctoral thesis (see
/// http://tdx.cat/handle/10803/81963).
/// When a skip_OXT parameter is needed, the function will use the OXT atom or not depending on its value. This
/// can greatly modify all results (starting from the COM centering) and makes it more difficult for comparison
/// with JR's iNMA software (imods.chaconlab.org).
///
/// build[First/Middle/Last]LinkUnits - Creates two plain units or a two special units if it is the first or last one.
/// build_[TYPE]_Unit - Builds one unit of this type.
/// get_[psi/phi]_dihedral - Returns the corresponding dihedral for one unit.
///
/// \author vgil
/// \date 1/07/2015
///////////////////////////////////////////////////////////////

UnitsBuilder::UnitsBuilder(Chain* from_this_chain, bool center_coordinates_at_com, bool skip_OXT){
	chain_topology = from_this_chain->getTopology();
	chain = from_this_chain;

	this->global_com = NULL;
	if(center_coordinates_at_com){
		this->global_com = centerChainCoordsAtCOM(skip_OXT);
	}

}

UnitsBuilder::~UnitsBuilder(){
	Utils::deleteAndSetPointerToNull(this->global_com);
}

void UnitsBuilder::build(std::vector<Unit*>& built_units ){
	vector<Unit*> tmp_units;
	create_CG_model(tmp_units, this->chain);
	merge_prolines(tmp_units , built_units);
}

void UnitsBuilder::create_CG_model(vector<Unit*>& units, Chain* chain){
	cout<<"DBG: Creating Coarse Grain Model"<<endl;
	// Pick all links
	vector<Link*> links = chain->getLinks();

	// Pick first link of the chain. As we don not trust in the
	// ordering of the links vector we use its linked list form.
	Link* first = NULL;
	for(unsigned int i = 0; i < links.size(); ++i){
		if(links[i]->isFirstInChain()){
			first = links[i];
		}
	}

	// Generate units from links (one link can generate more than one unit)
	Link* currentLink = first;
	while(currentLink != NULL){
		vector<Unit*> link_units = buildUnitsFromLink(currentLink);
		for (unsigned int j = 0; j < link_units.size(); ++j){
			units.push_back(link_units[j]);
		}
		currentLink = currentLink->getNextLink();
	}
}

void UnitsBuilder::merge_prolines(vector<Unit*>& in, vector<Unit*>& out){
	// Prolines are modeled with just one unit. We have to recheck the model.
	cout<<"DBG: Reducing model from "<<in.size()<<" units ";
	out.clear();

	for(unsigned i = 0; i < in.size(); ++i){
		// tmp_units.size() must be even (this CG model produces 2 units per residue).
		// As the second unit of each residue has atoms from residue i and i+1, when the next
		// unit is a full PRO unit, then it means that current unit owns the N atom of this residue.
		// as we want to freeze the N-CA torsion, we merge both units.
		// C-O=N and CA-SIDE
		if(i<in.size()-1 && i%2 == 1 && in[i+1]->resname == LinkNames::PRO){
			Unit* merged_unit = merge_units(in[i], in[i+1]);
			out.push_back(merged_unit);
			delete in[i];
			delete in[i+1];
			i++; // We do not want to add i+1 twice!
				 // TODO: There is surely a more elegant way to do this.
		}
		else{
			out.push_back(in[i]);
		}
	}
	cout<<" to "<< out.size() <<" units."<<endl;
}

vector<Unit*> UnitsBuilder::buildUnitsFromLink(Link* link){
	if (link->isFirstInChain()){ // Link is the first one
		return buildFirstLinkUnits(link);
	}
	else if(link->isLastInChain()){ // Link is the last one
		return buildLastLinkUnits(link);
	}
	else{ // Link is somewhere between first and last links
		return buildMiddleLinkUnits(link);
	}
}

void UnitsBuilder::getAtomsFromLink(vector<Atom* >& atoms, vector<string>& atom_names, Link* link){
	vector<Atom*> tmp_atoms;
	SelectionBuilder selection_builder;
	SelectionJsonBuilder selectionJsonBuilder;
	Selector* selector = selection_builder.create(selectionJsonBuilder.createAtomsByNameSelection(atom_names),"");
	selector->getAtoms(tmp_atoms, link);
	VectorTools::add(atoms, tmp_atoms);
	delete selector;
}

void UnitsBuilder::getAtomsAndSidechainFromLink(vector<Atom* >& atoms, vector<string>& atom_names, Link* link){
	getAtomsFromLink(atoms, atom_names, link);
	// Filter sidechain to get rid of the hydrogens
	const vector<Atom*> sidechainAtoms = link->getSideChainAtoms();
	for(unsigned int i = 0; i < sidechainAtoms.size();++i){
		if (!sidechainAtoms[i]->isHydrogen()){
			atoms.push_back(sidechainAtoms[i]);
		}
	}
}

void UnitsBuilder::filterHydrogens(vector<Atom*>& atoms, std::vector<Atom* >& filtered_atoms){
	for(unsigned int i = 0; i < atoms.size();++i){
		if (!atoms[i]->isHydrogen()){
			filtered_atoms.push_back(atoms[i]);
		}
	}
}

void UnitsBuilder::filterOXT(vector<Atom*>& atoms, std::vector<Atom* >& filtered_atoms){
	for(unsigned int i = 0; i < atoms.size();++i){
		if ( atoms[i]->name != AtomNames::OXT){
			filtered_atoms.push_back(atoms[i]);
		}
	}
}

void UnitsBuilder::get_all_hydrogens_from_N_Ca_Sidechain_Unit(Link* link, vector<Atom*>& atoms ){
	string backbone_hydrogen_names [] = {AtomNames::H1,
										 AtomNames::H2,
										 AtomNames::H3,
										 AtomNames::HA,
										 AtomNames::HA2,
										 AtomNames::HA3

	};
	vector<string> atom_names(backbone_hydrogen_names, backbone_hydrogen_names+6);

	getAtomsFromLink(atoms, atom_names, link);

	const vector<Atom*> sidechainAtoms = link->getSideChainAtoms();
	for(unsigned int i = 0; i < sidechainAtoms.size();++i){
		if (sidechainAtoms[i]->isHydrogen()){
			atoms.push_back(sidechainAtoms[i]);
		}
	}
}

void UnitsBuilder::get_all_hydrogens_from_Ca_Sidechain_Unit(Link* link, vector<Atom*>& atoms ){
	get_all_hydrogens_from_N_Ca_Sidechain_Unit(link, atoms );
}

void UnitsBuilder::get_all_hydrogens_from_C_O_O_Unit(Link* link, vector<Atom*>& atoms ){

}

void UnitsBuilder::get_all_hydrogens_from_N_C_O_Unit(Link* link, Link* next_link, vector<Atom*>& atoms ){
	vector<string> link1_hydrogen; link1_hydrogen.push_back(AtomNames::H);
	getAtomsFromLink(atoms, link1_hydrogen, next_link);
}

Unit* UnitsBuilder::build_N_Ca_Sidechain_Unit(Link* link){
	// Heavy atoms
	vector<Atom*> atoms;
	vector<string> names;
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::N));
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::CA));
	getAtomsAndSidechainFromLink(atoms,names,link);
	// Hydrogens
	vector<Atom*> hydrogens;
	get_all_hydrogens_from_Ca_Sidechain_Unit(link, hydrogens);

	// Get the other needed parameters to build the Unit
	// First unit only has one dihedral (psi)
	Atom *CA, *C;
	Dihedral* phi = NULL;
	Dihedral* psi = get_psi_dihedral(CA, C, link, link->isLastInChain());

	// Build it
	return new Unit(atoms,
			hydrogens,
			phi,// Left dihedral
			psi, // Right dihedral
			new Point(0,0,0),
			(new Point(Point::subtract(C->toPoint(),CA->toPoint())))->normalize(),
			new Point(CA->toPoint()),//Point(0,0,0),
			new Point(CA->toPoint()),
			CA,
			link->name,
			"N-CA-SIDE");
}

Unit* UnitsBuilder::build_Ca_Sidechain_Unit(Link* link){
	// Heavy atoms
	vector<Atom*> atoms;
	vector<string> names;
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::CA));
	getAtomsAndSidechainFromLink(atoms,names,link);
	// Hydrogens
	vector<Atom*> hydrogens;
	get_all_hydrogens_from_Ca_Sidechain_Unit(link, hydrogens);

	// Get the other needed parameters to build the Unit
	// Middle CA SD units have two dihedrals (phi, psi)
	Atom *CA, *C, *N;
	Dihedral* phi = get_left_phi_dihedral(N,CA, link);
	Dihedral* psi = get_psi_dihedral( CA, C, link, link->isLastInChain());
	return new Unit(atoms,
			hydrogens,
			phi,// Left dihedral
			psi, // Right dihedral
			(new Point(Point::subtract(N->toPoint(),CA->toPoint())))->normalize(),
			(new Point(Point::subtract(C->toPoint(),CA->toPoint())))->normalize(),
			new Point(CA->toPoint()),
			new Point(CA->toPoint()),
			CA,
			link->name,
			"CA-SIDE");
}

Unit* UnitsBuilder::build_C_O_O_Unit(Link* link){
	// Heavy atoms
	vector<Atom*> atoms;
	vector<string> names;
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::C));
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::O));
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::OXT));
	getAtomsFromLink(atoms,names,link);
	// Hydrogens
	vector<Atom*> hydrogens;
	get_all_hydrogens_from_C_O_O_Unit(link, hydrogens);

	// Get the other needed parameters to build the Unit
	// Middle COO units have 1 dihedral (psi)
	Atom *CA, *C;
	Dihedral* psi = get_psi_dihedral( CA, C, link, link->isLastInChain());
	Dihedral* phi = NULL;
	return new Unit(atoms,
			hydrogens,
			psi,// Left dihedral
			phi, // Right dihedral
			(new Point(Point::subtract(CA->toPoint(),C->toPoint())))->normalize(),
			new Point(0,0,0),
			new Point(C->toPoint()),
			new Point(0,0,0),
			C,
			link->name,
			"C-O=O");
}

Unit* UnitsBuilder::build_N_C_O_Unit(Link* link, Link* next_link){
	// Heavy atoms
	// Get C and O from first residue
	vector<Atom*> atoms;
	vector<string> names;
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::C));
	names.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::O));
	getAtomsFromLink(atoms,names,link);
	// Get N from the next one
	vector<string> names_next;
	names_next.push_back(StringTools::replaceSpacesByUnderscores(AtomNames::N));
	getAtomsFromLink(atoms, names_next, next_link);
	// Hydrogens
	vector<Atom*> hydrogens;
	get_all_hydrogens_from_N_C_O_Unit(link, next_link, hydrogens);


	// Get the other needed parameters to build the Unit
	// Middle N C O units have two dihedrals (psi, phi)
	Atom *CA, *C, *N_next, *CA_next;
	Dihedral* psi = get_psi_dihedral(CA, C, link, link->isLastInChain());
	Dihedral* phi = get_right_phi_dihedral(N_next, CA_next, link);
	return new Unit(atoms,
			hydrogens,
			psi,// Left dihedral
			phi, // Right dihedral
			(new Point(Point::subtract(CA->toPoint(),C->toPoint())))->normalize(),
			(new Point(Point::subtract(CA_next->toPoint(),N_next->toPoint())))->normalize(),
			new Point(C->toPoint()),
			new Point(N_next->toPoint()),
			C,
			link->name,
			"C=O-N");
}

vector<Unit*> UnitsBuilder::buildFirstLinkUnits(Link* link){
	// First link generates only one unit:
	// - N, Ca and Sidechain
	vector<Unit*> units;
	units.push_back(build_N_Ca_Sidechain_Unit(link));
	units.push_back(build_N_C_O_Unit(link, link->getNextLink()));
	return units;
}

vector<Unit*> UnitsBuilder::buildMiddleLinkUnits(Link* link){
	// Middle links generate two units:
	// - Ca, Sidechain
	// - C, O, N (from the next residue)
	vector<Unit*> units;
	units.push_back(build_Ca_Sidechain_Unit(link));
	units.push_back(build_N_C_O_Unit(link, link->getNextLink()));
	return units;
}

vector<Unit*> UnitsBuilder::buildLastLinkUnits(Link* link){
	// Last link generates two units:
	// - Ca, sidechain
	// - C, O, OXT
	vector<Unit*> units;
	units.push_back(build_Ca_Sidechain_Unit(link));
	units.push_back(build_C_O_O_Unit(link));
	return units;
}

Dihedral* UnitsBuilder::get_psi_dihedral( Atom*& CA, Atom*& C, Link* link, bool is_last){
	// In other links, psi is defined as:
	// N_before | CA C N
	Atom* N = link->getAtomsWithName(AtomNames::N)[0];
	CA = link->getAtomsWithName(AtomNames::CA)[0];
	C = link->getAtomsWithName(AtomNames::C)[0];
	Atom* last_atom;
	if (is_last){
		last_atom = link->getAtomsWithName(AtomNames::OXT)[0]; // OXT
	}
	else{
		last_atom = link->getNextLink()->getAtomsWithName(AtomNames::N)[0]; // N_next
	}
	Dihedral*  dihedral = chain_topology->locateDihedral(N,CA,C,last_atom);

	if (dihedral == 0) {
		cout << "DBG (get_psi) Null dihedral for atoms"
				<< *N
				<< *CA
				<< *C
				<< *last_atom << endl;
	}

	return dihedral;
}

Dihedral* UnitsBuilder::get_left_phi_dihedral( Atom*& N, Atom*& CA, Link* link){
	Atom* C_before = link->getPreviousLink()->getAtomsWithName(AtomNames::C)[0];
	N = link->getAtomsWithName(AtomNames::N)[0];
	CA = link->getAtomsWithName(AtomNames::CA)[0];
	Atom* C = link->getAtomsWithName(AtomNames::C)[0];

	Dihedral* dihedral = chain_topology->locateDihedral(C_before,(Atom*) N,(Atom*) CA, C);
	if (dihedral == 0) {
		cout << "DBG: (get_left_phi) Null dihedral for atoms"
				<< *C_before
				<< *N
				<< *CA
				<< *C << endl;
	}
	return dihedral;
}

Dihedral* UnitsBuilder::get_right_phi_dihedral( Atom*& N_next, Atom*& CA_next, Link* link){
	Atom* C = link->getAtomsWithName(AtomNames::C)[0];
	N_next = link->getNextLink()->getAtomsWithName(AtomNames::N)[0];
	CA_next = link->getNextLink()->getAtomsWithName(AtomNames::CA)[0];
	Atom* C_next = link->getNextLink()->getAtomsWithName(AtomNames::C)[0];

	Dihedral* dihedral = chain_topology->locateDihedral(C,(Atom*) N_next,(Atom*) CA_next, C_next);
	if (dihedral == 0) {
		cout << "DBG (get_right_phi) Null dihedral for atoms"
				<< *C
				<< *N_next
				<< *CA_next
				<< *C_next << endl;
	}
	return dihedral;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Moves the center of mass to 0,0,0 by subtracting the current center of mass to all atoms.
///
/// \param atoms [In] A list with all the heavy atoms the unit is composed of.
///
/// \author vgil
/// \date 01/07/2015
///////////////////////////////////////////////////////////////
CenterOfMass* UnitsBuilder::centerChainCoordsAtCOM(bool skip_OXT){
	cout<<"DBG: Centering molecule  [centering vector: "; // CAUTION: if you want to erase this line, erase the last one too

	vector<Atom*> non_hydrogens, non_hydrogens_no_oxt, *used_atoms;

	filterHydrogens(chain->getAllAtoms(), non_hydrogens);

	if (skip_OXT){
		filterOXT(non_hydrogens, non_hydrogens_no_oxt);
		used_atoms = &non_hydrogens_no_oxt;
		cout<<"DBG: COM IS NOT USING OXT ATOM !!"<<endl;
	}
	else{
		used_atoms = &non_hydrogens;
	}

	CenterOfMass* com = new CenterOfMass(CenterOfMass::compute(*(used_atoms)));
	vector<double> com_coords = com->getCoordinates();

	// We have to shift ALL atoms in the system (including ligands etc)
	Math::shift3D(chain->getNumberOfAllAtoms(),
			Utils::vectorToPointer<double>(chain->getAllCartesianCoordinates()),
			Utils::vectorToPointer<double>(com_coords),
			-1.);

	cout<<setprecision(16)<<com_coords[0]<<" "<<com_coords[1]<<" "<<com_coords[2]<<" ] "<<endl;

	return com;
}
///////////////////////////////////////////////////////////////
/// \remarks
/// Merges two units to create a third one. Both units must be linked or the result will be unexpected.
//  TODO: Must be moved to UnitTools (?)
///
/// \param first_unit, second_unit [In] Units to merge
///
/// \return The merged unit.
///
/// \author vgil
/// \date 1/07/2015
///////////////////////////////////////////////////////////////
Unit* UnitsBuilder::merge_units(Unit* first_unit, Unit* second_unit){

		vector<Atom*> heavy_atoms;
		VectorTools::add(heavy_atoms,first_unit->atoms);
		VectorTools::add(heavy_atoms,second_unit->atoms);

		vector<Atom*> hydrogens;
		VectorTools::add(hydrogens,first_unit->hydrogens);
		VectorTools::add(hydrogens,second_unit->hydrogens);

		Unit* merged_unit = new Unit(heavy_atoms,
				hydrogens,
				first_unit->getLeftDihedral(),
				second_unit->getRightDihedral(),
				new Point(*(first_unit->e_left)),
				new Point(*(second_unit->e_right)),
				new Point(*(first_unit->r_left)),
				new Point(*(second_unit->r_right)),
				second_unit->main_atom,
				second_unit->resname,
				first_unit->toString()+" + "+second_unit->toString());

		return merged_unit;
}


CenterOfMass* UnitsBuilder::getCOM(){
	return this->global_com;
}
