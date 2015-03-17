/////////////////////////////////////////////////////////////////////////////
/// AnmTargetUpdater.cpp
///
/// Implementation of AnmTargetUpdater class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 13/09/2012
/////////////////////////////////////////////////////////////////////////////

#include "AnmTargetUpdater.h"

#include "../Tools/AnmMagnitudeCalculator.h"
#include "../Tools/AnmMoveVectorCalculator.h"
#include "../../ANM/Parameters/AnmParameters.h"
#include "../../Molecules/Atom.h"
#include "../../Tools/Utils.h"
#include <algorithm>
#include <sstream>
#include <ostream>
#include "../../Molecules/AtomSet/AtomSet.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// Constructor
///
/// \param magnitudeCalculator [In] Collaborator that computes the magnitude of the displacement
/// \param moveVectorCalculator [In] Collaborator that computes the movement vector
///
/// \author mrivero
/// \date 13/09/2012
///////////////////////////////////////////////////////////////
AnmTargetUpdater::AnmTargetUpdater(AnmMagnitudeCalculator * magnitudeCalculator,
		 	 	 	 	 	 	 	 AnmMoveVectorCalculator * moveVectorCalculator)
{
	this->magnitudeCalculator = magnitudeCalculator;
	this->moveVectorCalculator = moveVectorCalculator;
}

AnmTargetUpdater::~AnmTargetUpdater()
{
	Utils::deleteAndSetPointerToNull(magnitudeCalculator);
	Utils::deleteAndSetPointerToNull(moveVectorCalculator);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It generates the target coordinate updater report
///
/// \param anmParameters [In] ANM parameters
///
/// \return The target coordinate updater report
///
/// \author mrivero
/// \date 23/04/2013
///////////////////////////////////////////////////////////////
std::string AnmTargetUpdater::generateReport(AnmParameters * anmParameters) const
{
	ostringstream oss(ostringstream::out);

	oss<<moveVectorCalculator->generateReport(anmParameters)<<endl;
	oss<<magnitudeCalculator->generateReport();

	return oss.str();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It computes the displacement magnitude (and sets it in the ANM parameters)
///
/// \param eigen [In] Eigen vectors and eigen values
/// \param chosenModeIndex [In] Index of chosen mode
///
/// \param anmParameters [In/Out] ANM parameters
///
/// \author mrivero
/// \date 16/12/2013
///////////////////////////////////////////////////////////////
void AnmTargetUpdater::computeDisplacementMagnitude(AnmParameters * anmParameters, AnmEigen * eigen,
															unsigned int chosenModeIndex)
{
	double magnitude = magnitudeCalculator->computeDisplacementMagnitude(anmParameters, eigen, chosenModeIndex);
	anmParameters->setDisplacementMagnitude(magnitude);
}
