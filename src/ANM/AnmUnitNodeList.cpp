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
	sourceAtomSet = NULL;
}

std::vector<Unit*> AnmUnitNodeList::getNodeList() const {
	return this->nodeList;
}

void AnmUnitNodeList::setNodeList(AtomSet* atomSet) {
	this->sourceAtomSet = atomSet;
}
void AnmUnitNodeList::setNodeList(vector<Unit*>& units) {
	this->nodeList = units;
}

void AnmUnitNodeList::updateUnitList(){
	Utils::clearVector<Unit>(nodeList);
	UnitsBuilder unitsBuilder((Chain*) sourceAtomSet);
	unitsBuilder.build(nodeList);
}


string AnmUnitNodeList::showNodeList(){
	return LogUtils::showUnitLabels(this->nodeList);
}


unsigned int AnmUnitNodeList::size(){
	return this->nodeList.size();
}
