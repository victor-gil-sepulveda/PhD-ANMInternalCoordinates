/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 02/10/2014
/////////////////////////////////////////////////////////////////////////////
#include "ModesWriter.h"
#include "../../ModesCalculator/AnmEigen.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../../ModesCalculator/Internals/CoarseGrainModel/UnitTools.h"
#include "../../../Molecules/AtomNames.h"
#include "../../../Molecules/Atom.h"
#include "../../ModesCalculator/Internals/InternalModesCalculator.h"
#include "../../../Tools/Math/MathTools.h"
#include "../../AnmAtomNodeList.h"
#include "../../AnmUnitNodeList.h"
using namespace std;

ModesWriter::ModesWriter() {
	this->full_path = "";
	this-> writer_name = "";
}

ModesWriter::ModesWriter(string writer_name, string full_path) {
	this->full_path = full_path;
	this-> writer_name = writer_name;
}

ModesWriter::~ModesWriter() {
}

void ModesWriter::setPath(std::string full_path){
	this->full_path = full_path;
}

void ModesWriter::setName(std::string name){
	this->writer_name = name;
}


void ModesWriter::writeCartesianModes(AnmEigen * eigen, vector<double>& coordinates){
	ofstream file_handler;
	file_handler.open(this->full_path.c_str());
	cout<<"DBG: Writing "<<full_path<<" ..."<<endl;

	writeHeader(file_handler, eigen->type);
	writeCoordinates(file_handler, coordinates);

	for (unsigned int i = 0; i < eigen->numberOfModes;++i){
		writeOneMode(file_handler, i+1, eigen->getEigenValueOfMode(i), eigen->getEigenVectorOfMode(i));
	}

	file_handler.close();
}

void ModesWriter::writeCartesianModes(AnmEigen * eigen, std::vector<Atom*>& atoms){
	ofstream file_handler;
	cout<<"DBG: Writing "<<full_path<<" ..."<<endl;
	file_handler.open(this->full_path.c_str());

	writeHeader(file_handler, eigen->type);
	writeAtomnames(file_handler, atoms);
	writeResids(file_handler, atoms);
	writeResnames(file_handler, atoms );
	writeCoordinates(file_handler, atoms);

	for (unsigned int i = 0; i < eigen->numberOfModes;++i){
		writeOneMode(file_handler, i+1, eigen->getEigenValueOfMode(i), eigen->getEigenVectorOfMode(i));
	}

	file_handler.close();
}

void ModesWriter::writeCartesianModes(AnmEigen * eigen, const AnmNodeList* nodeList){
	vector<Atom*> atoms = dynamic_cast<const AnmAtomNodeList*>( nodeList)->getNodeList();
	writeCartesianModes(eigen, atoms);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Writes the contents of a AnmEigen object into a file. It tries to be compatible with
/// Prody and the vmd viewer.
///
/// \param filename [In] The path of the file were we will write the modes.
/// \param eigen [In] Modes to write.
/// \param units [In] The units of the CG model used to get the eigenvectors.
/// \param toCartesian [In] If true, the modes will be converted to cartesian modes prior to writing.
///
/// \author vgil
/// \date 7/01/2015
///////////////////////////////////////////////////////////////
void ModesWriter::writeInternalModes( AnmEigen * eigen, std::vector<Unit*>& units,
		bool toCartesian, ModeWritingType type){
	if (toCartesian){
		cout<<"DBG: converting to cartesian prior to writing"<<endl;
		cout<<"DBG: first eigenvalue values ("<<eigen->vectors[0].size()<< ") "<<
				eigen->vectors[0][0]<<" "<<
				eigen->vectors[0][1]<<" "<<
				eigen->vectors[0][2]<<" "<<
				eigen->vectors[0][3]<<" "<<
				eigen->vectors[0][4]<<" "<<endl;
		// First convert from internal to cartesian
		bool onlyHeavyAtoms = (not type == ALL_ATOMS);
		AnmEigen* converted = InternalModesCalculator::internalToCartesian(units, eigen, onlyHeavyAtoms);

		cout<<"DBG: first converted values ("<<converted->vectors[0].size()<< ") "<<
				converted->vectors[0][0]<<" "<<
				converted->vectors[0][1]<<" "<<
				converted->vectors[0][2]<<" "<<
				converted->vectors[0][3]<<" "<<
				converted->vectors[0][4]<<" "<<endl;
		if (type == CA_ATOMS){
			vector<int> ca_atom_indices;
			vector<Atom*> ca_atoms;
			getCAsFromUnits(units, ca_atom_indices, ca_atoms);

			// Then filter the modes
			AnmEigen filtered;
			filterEigenvectors(converted, &filtered, ca_atom_indices);

			// Write them!
			this->writeCartesianModes( &filtered, ca_atoms);
		}
		else{
			vector<Atom*> heavy_atoms;
			UnitTools::getAllAtomsFromUnits(units, heavy_atoms, true);
			this->writeCartesianModes( converted, heavy_atoms);
		}

		delete converted;
	}
	else{
		ofstream file_handler;
		file_handler.open(full_path.c_str());
		cout<<"DBG: Writing "<<full_path<<" ..."<<endl;

		writeHeader(file_handler, eigen->type);

		vector<Atom*> bb_atoms;
		getBBAtomsFromUnits(units, bb_atoms);

		writeAtomnames(file_handler, bb_atoms);
		writeResids(file_handler, bb_atoms);
		writeResnames(file_handler, bb_atoms);
		writeCoordinates(file_handler, bb_atoms);

		for (unsigned int i = 0; i < eigen->numberOfModes;++i){
			writeOneMode(file_handler, i+1, eigen->getEigenValueOfMode(i), eigen->getEigenVectorOfMode(i));
		}

		file_handler.close();
	}
}

void ModesWriter::writeInternalModes(AnmEigen * eigen, const  AnmNodeList* nodeList, bool toCartesian, ModeWritingType type){
	vector<Unit*> units = dynamic_cast<const AnmUnitNodeList*>(nodeList)->getNodeList();

	writeInternalModes( eigen, units, toCartesian, type);
}

void ModesWriter::writeHeader(std::ofstream& file_handler, ModeTypes::ModeType type){
	file_handler << "nmwiz_load " << this->full_path << endl;
	file_handler << "name " << this->writer_name << endl;
	file_handler << "type " << ModeTypes::toString(type) << endl;
}

void ModesWriter::writeResnames(ofstream& file_handler, vector<Atom*>& atoms){
	file_handler << "resnames " ;
	for (unsigned int i = 0; i< atoms.size(); ++i){
		file_handler<< atoms[i]->resName<<" ";
	}
	file_handler << endl;
}

void ModesWriter::writeAtomnames(ofstream& file_handler, vector<Atom*>& atoms){
	file_handler << "atomnames " ;
	for (unsigned int i = 0; i< atoms.size(); ++i){
		file_handler<< atoms[i]->name<<" ";
	}
	file_handler << endl;
}

void ModesWriter::writeResids(ofstream& file_handler, vector<Atom*>& atoms){
	file_handler << "resids " ;
	for (unsigned int i = 0; i< atoms.size(); ++i){
		file_handler<< atoms[i]->resSeq<<" ";
	}
	file_handler << endl;
}

void ModesWriter::writeCoordinates(ofstream& file_handler, vector<double>& coordinates){
	file_handler << "coordinates " ;
	for (unsigned int i = 0; i< coordinates.size(); ++i){
		file_handler<< coordinates[i]<<" ";
	}
	file_handler << endl;
}

void ModesWriter::writeCoordinates(ofstream& file_handler, vector<Atom*>& atoms){
	file_handler << "coordinates " ;
	for (unsigned int i = 0; i< atoms.size(); ++i){
		file_handler<< atoms[i]->getX()<<" "<<atoms[i]->getY()<<" "<<atoms[i]->getZ()<<" ";
	}
	file_handler << endl;
}

void ModesWriter::writeOneMode(ofstream& file_handler, int mode_index, double eigenvalue, vector<double>& eigenvector){
	file_handler<<"mode "<<mode_index<<" "<<eigenvalue<<"    ";
	for(unsigned int i = 0; i < eigenvector.size(); ++i){
		file_handler << setprecision(10) << eigenvector[i] << " ";
	}
	file_handler <<endl;
}

void ModesWriter::getCAsFromUnits(vector<Unit*>& units, vector<int>& ca_atom_indices, vector<Atom*>& ca_atoms){
	bool onlyHeavyAtoms = true;
	vector<Atom*> unit_atoms;
	UnitTools::getAllAtomsFromUnits(units, unit_atoms, onlyHeavyAtoms);
	for (unsigned int i =0; i< unit_atoms.size();++i){
		if(unit_atoms[i]->name == AtomNames::CA){
			ca_atom_indices.push_back(i);
			ca_atoms.push_back(unit_atoms[i]);
		}
	}
}

void ModesWriter::getBBAtomsFromUnits(vector<Unit*>& units, vector<Atom*>& bb_atoms){
	bool onlyHeavyAtoms = true;
	vector<Atom*> unit_atoms;
	UnitTools::getAllAtomsFromUnits(units, unit_atoms, onlyHeavyAtoms);
	for (unsigned int i =0; i< unit_atoms.size();++i){
		if(unit_atoms[i]->name == AtomNames::CA || unit_atoms[i]->name == AtomNames::N
				|| unit_atoms[i]->name == AtomNames::C){
			bb_atoms.push_back(unit_atoms[i]);
		}
	}
}

void ModesWriter::filterEigenvectors(AnmEigen* original,
		AnmEigen* filtered,
		vector<int>& ca_atom_indices){

	std::vector<std::vector<double> > filtered_evectors;
	vector<std::vector<double> >& original_evectors = original->vectors;

	for (unsigned int i = 0; i < original_evectors.size(); ++i){
		vector<double> one_eigen_filtered;
		for (unsigned int j = 0; j < ca_atom_indices.size(); ++j){
			int index_offset = ca_atom_indices[j] * 3;
			one_eigen_filtered.push_back(original_evectors[i][index_offset]);
			one_eigen_filtered.push_back(original_evectors[i][index_offset+1]);
			one_eigen_filtered.push_back(original_evectors[i][index_offset+2]);
		}
		filtered_evectors.push_back(one_eigen_filtered);
	}

	bool isCartesian = true;
	filtered->initialize(original->values, filtered_evectors, isCartesian);
}


///////////////////////////////////////////////////////////////
/// \remarks
/// Extracts the CA atoms of the complex.
///
/// \author vgil
/// \date 17/02/2015
///////////////////////////////////////////////////////////////
void ModesWriter::getCaAtoms(vector<Atom*>& source_atoms , vector<Atom*>& atoms){
	for (unsigned int i = 0; i < source_atoms.size(); ++i){
		Atom* atom = source_atoms[i];
		if(atom->name == AtomNames::CA){
			atoms.push_back(atom);
		}
	}
}

//////////////////////////////////////////////////////////////
/// \remarks
/// Generates an array of coordinates from the atom coordinates given as parameter.
/// \author vgil
/// \date 18/02/2015
///////////////////////////////////////////////////////////////
void ModesWriter::getCoordinatesFromAtomVector(std::vector<Atom*>& atoms, std::vector<double>& coords){
	coords.clear();
	for (unsigned int i = 0; i < atoms.size(); ++i){
		Atom* atom = atoms[i];
		coords.push_back(atom->getX());
		coords.push_back(atom->getY());
		coords.push_back(atom->getZ());
	}
}

