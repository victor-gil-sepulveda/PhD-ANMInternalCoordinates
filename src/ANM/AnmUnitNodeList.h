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
		std::vector<Unit*> getNodeList() const;
		void setNodeList(AtomSet* atomSet);
		std::string showNodeList();
		unsigned int size();

		void updateUnitList();

		// TODO: must handle Unit deletion when destructed?

	private:
		std::vector<Unit*> nodeList;
		AtomSet* sourceAtomSet;
		bool updated;
};




#endif /* ANMUNITNODELIST_H_ */
