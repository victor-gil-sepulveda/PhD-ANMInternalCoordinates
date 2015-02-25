/*
 * AnmUnitNodeList.cpp
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */

#include "AnmUnitNodeList.h"
#include <string>
#include "../PELE/PeleTasks/Output/LogUtils.h"
using namespace std;

std::vector<Unit*> AnmUnitNodeList::getNodeList() const {
	return this->nodeList;
}

void AnmUnitNodeList::setNodeList(const std::vector<Unit*> &nodeList) {
	this->nodeList = nodeList;
}


string AnmUnitNodeList::showNodeList(){
	return LogUtils::showUnitLabels(this->nodeList);
}


unsigned int AnmUnitNodeList::size(){
	return this->nodeList.size();
}
