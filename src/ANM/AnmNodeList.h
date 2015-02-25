/*
 * AnmNodeList.h
 *
 *  Created on: 28/07/2014
 *      Author: jestrada
 */

#ifndef ANMNODELIST_H_
#define ANMNODELIST_H_

#include <vector>
#include <string>

// abstract parent class of all ANM list of nodes' types
// (such as vector of atoms or vector of units).
// It serves as a way to do generic algorithms that
// let the specifics of dealing with the node list to
// internal algorithms, which must cast the AnmNodeList
// to the specific class (AnmAtomNodeList or AnmUnitNodeList)
// before operating.
class AnmNodeList {
	public:
		virtual ~AnmNodeList() {};

		virtual std::string showNodeList() = 0;

		virtual unsigned int size() = 0;
};



#endif /* ANMNODELIST_H_ */
