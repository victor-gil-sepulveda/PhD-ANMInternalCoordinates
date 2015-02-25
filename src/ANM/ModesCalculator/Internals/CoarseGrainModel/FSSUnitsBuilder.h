/////////////////////////////////////////////////////////////////////////////
/// Frozen Secondary Structure Units Builder
///
/// \author myName
/// \date 03/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef FSSUNITSBUILDER_H_
#define FSSUNITSBUILDER_H_

#include "UnitsBuilder.h"
#include <string>


// Frozen Secondary Structure Units Builder

class FSSUnitsBuilder: public UnitsBuilder {
	public:
		FSSUnitsBuilder(Chain* from_this_chain,std::string secondary_structure_description);
		virtual ~FSSUnitsBuilder();

		void build(std::vector<Unit*>& built_units);

	private:
		std::string ssd; //secondary structure description

};

#endif /* FSSUNITSBUILDER_H_ */
