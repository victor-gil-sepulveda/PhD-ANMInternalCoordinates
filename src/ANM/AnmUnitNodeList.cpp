/*
 * AnmUnitNodeList.cpp
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */

#include "AnmUnitNodeList.h"
#include <string>
#include "../PELE/PeleTasks/Output/LogUtils.h"
#include "ModesCalculator/Internals/CoarseGrainModel/UnitsBuilder.h"
#include "../Tools/Utils.h"
#include "ModesCalculator/Internals/CoarseGrainModel/Unit.h"
using namespace std;

AnmUnitNodeList::AnmUnitNodeList(){
}

AnmUnitNodeList::~AnmUnitNodeList(){
	Utils::clearVector(this->nodeList);
}

std::vector<Unit*> AnmUnitNodeList::getNodeList() const {
	return this->nodeList;
}

void AnmUnitNodeList::centerAtCOM(){
	//TODO: IMPLEMENT!
	cout<<"***ATENTION****"<<endl<<"AnmUnitNodeList::centerAtCOM must be implemented"<<endl<<"*********"<<endl;
}

void AnmUnitNodeList::setNodeList(vector<Unit*>& units) {
	this->nodeList = units;
}

void AnmUnitNodeList::updateUnitList(){
	cout<<"DBG: Updating CG Model"<<endl;
	for(unsigned int i = 0; i< nodeList.size();++i){
		nodeList[i]->update();
	}
}

string AnmUnitNodeList::showNodeList(){
	return LogUtils::showUnitLabels(this->nodeList);
}


unsigned int AnmUnitNodeList::size(){
	return this->nodeList.size();
}
