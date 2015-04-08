
#pragma once

#ifndef MODESWRITER_H_
#define MODESWRITER_H_
////////////////////////////////////////////////////////////////////////////
/// \brief
/// Very basic class that writes the contents of a AnmEigen object into a file (Prody compatible).
/// Some features (like writing the real coordinates of units) are missing.
/// \author vgil
/// \date 7/01/2015
/////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iosfwd>
#include "../../ModesCalculator/PreCalculated/ModeTypes.h"

class AnmEigen;
class Unit;
class Atom;
class AnmNodeList;

enum ModeWritingType{
	CA_ATOMS,
	HEAVY_ATOMS,
	ALL_ATOMS
};

class ModesWriter {
	public:
		ModesWriter();
		ModesWriter(std::string writer_name, std::string full_path);
		virtual ~ModesWriter();

		void setPath(std::string full_path);
		void setName(std::string name);

		void 	writeCartesianModes( AnmEigen * eigen, std::vector<double>& coordinates);
		void 	writeCartesianModes(AnmEigen * eigen, const AnmNodeList* nodeList);
		void 	writeCartesianModes(AnmEigen * eigen, std::vector<Atom*>& atoms);
		void 	writeInternalModes( AnmEigen * eigen, const AnmNodeList* nodeList, bool toCartesian, ModeWritingType type = CA_ATOMS);
		void 	writeInternalModes( AnmEigen * eigen, std::vector<Unit*>& units, bool toCartesian, ModeWritingType type = CA_ATOMS);

		void 	writeHeader(std::ofstream& file_handler, ModeTypes::ModeType type);
		void 	writeResnames(std::ofstream& file_handler, std::vector<Atom*>& atoms);
		void 	writeAtomnames(std::ofstream& file_handler, std::vector<Atom*>& atoms);
		void 	writeResids(std::ofstream& file_handler, std::vector<Atom*>& atoms);
		void 	writeCoordinates(std::ofstream& file_handler, std::vector<Atom*>& atoms);
		void 	writeCoordinates(std::ofstream& file_handler, std::vector<double>& coordinates);
		void 	writeOneMode(std::ofstream& file_handler, int mode_number, double eigenvalue, std::vector<double>& eigenvector);

static void 	getCaAtoms(std::vector<Atom*>& source_atoms, std::vector<Atom*>& atoms);
static void 	getCAsFromUnits(std::vector<Unit*>& units, std::vector<int>& ca_atom_indices, std::vector<Atom*>& ca_atoms);
static void 	getBBAtomsFromUnits(std::vector<Unit*>& units, std::vector<Atom*>& bb_atoms);
static void 	filterEigenvectors(AnmEigen* original, AnmEigen* filtered, std::vector<int>& ca_atom_indices);
static void 	getCoordinatesFromAtomVector(std::vector<Atom*>& atoms, std::vector<double>& coords);


	private:
		std::string full_path, writer_name;
};

#endif /* MODESWRITER_H_ */
