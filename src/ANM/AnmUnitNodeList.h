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
class AtomSet;
class AnmUnitNodeList : public AnmNodeList {
	public:
		AnmUnitNodeList();
		~AnmUnitNodeList();

		std::vector<Unit*> getNodeList() const;
		void setNodeList(std::vector<Unit*>& units);

		std::string showNodeList();
		unsigned int size();

		void updateUnitList();

	private:
		std::vector<Unit*> nodeList;
		bool updated;
};




#endif /* ANMUNITNODELIST_H_ */
