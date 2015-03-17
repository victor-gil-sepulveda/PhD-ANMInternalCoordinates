/*
 * AnmAtomNodeList.cpp
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */


#include "AnmAtomNodeList.h"
#include <string>
#include "../PELE/PeleTasks/Output/LogUtils.h"
using namespace std;

std::vector<Atom*> AnmAtomNodeList::getNodeList() const {
	return this->nodeList;
}

void AnmAtomNodeList::setNodeList(const std::vector<Atom*> &nodeList) {
	this->nodeList = nodeList;
}

string AnmAtomNodeList::showNodeList(){
	return LogUtils::showAtomsLabels(this->nodeList);
}

unsigned int AnmAtomNodeList::size(){
	return this->nodeList.size();
}
