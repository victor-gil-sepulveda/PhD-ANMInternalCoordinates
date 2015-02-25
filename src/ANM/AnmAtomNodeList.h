/*
 * AnmAtomNodeList.h
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */

#ifndef ANMATOMNODELIST_H_
#define ANMATOMNODELIST_H_


#include <vector>

#include "AnmNodeList.h"
#include <string>

class Atom;

class AnmAtomNodeList : public AnmNodeList {
	public:
		std::vector<Atom*> getNodeList() const;
		void setNodeList(const std::vector<Atom*> &nodeList);
		std::string showNodeList();
		unsigned int size();

	private:
		std::vector<Atom*> nodeList;
};



#endif /* ANMATOMNODELIST_H_ */
