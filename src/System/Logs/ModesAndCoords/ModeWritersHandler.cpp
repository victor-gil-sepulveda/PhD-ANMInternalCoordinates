/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 18/02/2015
/////////////////////////////////////////////////////////////////////////////
#include "ModeWritersHandler.h"
#include "../../../ANM/Tools/Inout/ModesWriter.h"
#include <string>
#include <map>
#include <sstream>
#include "../../../Molecules/AtomSet/AtomSet.h"
#include <iostream>
#include "../../../System/SystemVars.h"
#include "../../../Tools/Math/MathTools.h"
#include "../../../System/Logs/OutputWriter.h"
#include "../../../ANM/Algorithms/AnmInternals/Movement/ANMICMovement.h"
#include "../../../ANM/Algorithms/AnmInternals/AnmInternals.h"
#include "../../../Tools/Utils.h"

using namespace std;

ModesWriterHandler::ModesWriterHandler(){
	counter = 0;
	directory = "";
}

ModesWriterHandler::~ModesWriterHandler() {
}

void ModesWriterHandler::setDirectory(string directory){
	this->directory = directory;
	cout<<"DBG: Setting mode logging directoy to: "<<directory<<endl;
}

void ModesWriterHandler::setCounter(unsigned int c){
	this->counter = c;
}

std::string ModesWriterHandler::generatePath(std::string name){
	stringstream ss;
	ss<<this->directory<<(this->directory!=""?string("/"):string(""))<<name<<"."<<this->counter<<".nmd";
	return ss.str();
}

void ModesWriterHandler::setAtomSet(AtomSet* aset){
	atomset = aset;
}

ModesWriter* ModesWriterHandler::getWriter(std::string name){
	writer.setPath(generatePath(name));
	writer.setName(name);
	return &writer;
}

void ModesWriterHandler::logStepAndVector(string logger_name, vector<double>& v){
	stringstream oss;

	oss<<this->counter<<" ";
	for (unsigned int i = 0; i< v.size(); ++i){
		oss<<v[i]<<" ";
	}

	SystemVars::getLog(logger_name,this->directory)->write(oss.str());
	SystemVars::flushLogs();
}

void ModesWriterHandler::logStepAndVector(string logger_name, vector<unsigned int>& v){
	stringstream oss;

	oss<<this->counter<<" ";
	for (unsigned int i = 0; i< v.size(); ++i){
		oss<<v[i]<<" ";
	}

	SystemVars::getLog(logger_name,this->directory)->write(oss.str());
	SystemVars::flushLogs();
}

void ModesWriterHandler::logStepAndCurrentCACoordinates(string logger_name){
	vector<Atom*> atoms;
	vector<double> ca_coords;
	ModesWriter::getCaAtoms(this->atomset->getAllAtoms(), atoms);
	ModesWriter::getCoordinatesFromAtomVector(atoms, ca_coords);
	logStepAndVector(logger_name, ca_coords);
}

void ModesWriterHandler::logStepAndCAProposal(string logger_name, vector<double>& coord_increments){
	vector<Atom*> atoms;
	vector<double> ca_coords;
	ModesWriter::getCaAtoms(this->atomset->getAllAtoms(), atoms);
	ModesWriter::getCoordinatesFromAtomVector(atoms, ca_coords);

	vector<double> proposal(ca_coords.size(), 0);
	Math::addVectors(ca_coords.size(),
				Utils::vectorToPointer(ca_coords),
				Utils::vectorToPointer(coord_increments),
				Utils::vectorToPointer(proposal));

	logStepAndVector(logger_name, proposal);
}

void ModesWriterHandler::logStepAndCurrentDihedralAngles(string logger_name, vector<Unit*>& units){
	vector<double> current_angles;
	AnmInternals::calculate_current_angles(current_angles, units);
	logStepAndVector(logger_name, current_angles);
}

void  ModesWriterHandler::logStepAndDihedralAnglesProposal(string logger_name, std::vector<double>& target_angle_increments, vector<Unit*>& units){
	vector<double> current_angles;
	AnmInternals::calculate_current_angles(current_angles, units);

	vector<double> proposal(current_angles.size(), 0);
	Math::addVectors(current_angles.size(),
			Utils::vectorToPointer(current_angles),
			Utils::vectorToPointer(target_angle_increments),
			Utils::vectorToPointer(proposal));

	logStepAndVector(logger_name, proposal);
}
void ModesWriterHandler::logStepAndDihedralToCartesianProposal(std::string logger_name,
		std::vector<double>& target_angle_increments, std::vector<Unit*>& units){
	vector<double> coords, icoords;
	atomset->saveCoordinates(coords,icoords);
	ANMICMovement::apply_rotations_to_molecule_units(units, target_angle_increments, 1.);
	logStepAndCurrentCACoordinates(logger_name);
	atomset->restoreCoordinates(coords, icoords);
}

void ModesWriterHandler::saveCoordinates(string coordset_name, vector<double>& coords){
	this->coordinates_bank[coordset_name] = coords;
}

void ModesWriterHandler::saveCACoordinates(std::string coordset_name){
	vector<Atom*> atoms;
	vector<double> ca_coords;
	ModesWriter::getCaAtoms(this->atomset->getAllAtoms(), atoms);
	ModesWriter::getCoordinatesFromAtomVector(atoms, ca_coords);
	saveCoordinates(coordset_name, ca_coords);
}

vector<double>& ModesWriterHandler::getCoordinates(string coordset_name){
	return this->coordinates_bank[coordset_name];
}

void ModesWriterHandler::removeCoordinates(std::string coordset_name){
	this->coordinates_bank.erase(coordset_name);
}
