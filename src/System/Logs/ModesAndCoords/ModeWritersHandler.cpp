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

using namespace std;

ModesWriterHandler::ModesWriterHandler() {
	counter = 0;
	directory = "";
}

ModesWriterHandler::~ModesWriterHandler() {
	for(std::map<string, ModesWriter*>::iterator iter = writers.begin(); iter != writers.end(); ++iter){
		delete iter->second;
	}
}

void ModesWriterHandler::setDirectory(string directory){
	this->directory = directory;
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

	map<string, ModesWriter*>::iterator it = writers.find(name);

	// If it's null, create it
	if(it == writers.end()){
		ModesWriter* mywriter = new ModesWriter(name, generatePath(name));
		writers[name]= mywriter;
		return mywriter;
	}
	// If not... return the already created logger.
	else{
		ModesWriter* mywriter = writers[name];
		mywriter->setPath(generatePath(name));
		return writers[name];
	}
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

