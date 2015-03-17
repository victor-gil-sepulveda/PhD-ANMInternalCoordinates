/////////////////////////////////////////////////////////////////////////////
/// Frozen Secondary Structure Units Builder
///
/// \author vgil
/// \date 03/02/2015
/////////////////////////////////////////////////////////////////////////////

#include "FSSUnitsBuilder.h"
#include <string>
using namespace std;

FSSUnitsBuilder::FSSUnitsBuilder(Chain* from_this_chain, string secondary_structure_description): UnitsBuilder(from_this_chain, true, false) {
	this->ssd = secondary_structure_description;
}

FSSUnitsBuilder::~FSSUnitsBuilder() {
}

void FSSUnitsBuilder::build(std::vector<Unit*>& built_units){
	UnitsBuilder::build(built_units);
}


