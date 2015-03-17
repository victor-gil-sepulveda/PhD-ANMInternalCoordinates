/////////////////////////////////////////////////////////////////////////////
/// \file AnmAlgorithmTypes.h
///
/// \brief Different types of ANM algorithm and ANM algorithm constants
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 01/11/2011
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMALGORITHMTYPES_H_
#define ANMALGORITHMTYPES_H_

#include <string>

/////////////////////////////////////////////////////////////////////////////
/// \brief Different ANM algorithm types
///
/// \author mrivero
/// \date 01/11/2011
/////////////////////////////////////////////////////////////////////////////
enum AnmAlgorithmType
{
	ANM_CARTESIAN_ALPHACARBONS = 1,
	ANM_CARTESIAN,
	ANM_INTERNAL
};

/////////////////////////////////////////////////////////////////////////////
/// \brief Some ANM predefined node selections
///
/// \author mrivero
/// \date 01/11/2011
/////////////////////////////////////////////////////////////////////////////
namespace AnmAlgorithmPredefinedSelections
{
	const std::string ALPHA_CARBONS_SELECTION_STRING = "{\"atoms\" : { \"names\":[\"_CA_\"]} }";
	const std::string EVERYTHING_BUT_LIGAND = "{ \"chains\" : \"everythingButLigands\" }";
}

#endif /* ANMALGORITHMTYPES_H_ */
