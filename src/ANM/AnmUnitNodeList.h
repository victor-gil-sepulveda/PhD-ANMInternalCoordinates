/*
 * AnmUnitNodeList.h
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */

#ifndef ANMUNITNODELIST_H_
#define ANMUNITNODELIST_H_

#include <vector>

#include "AnmNodeList.h"

class Unit;

class AnmUnitNodeList : public AnmNodeList {
	public:
		std::vector<Unit*> getNodeList() const;
		void setNodeList(const std::vector<Unit*> &nodeList);
		std::string showNodeList();
		unsigned int size();

	private:
		std::vector<Unit*> nodeList;
};




#endif /* ANMUNITNODELIST_H_ */
