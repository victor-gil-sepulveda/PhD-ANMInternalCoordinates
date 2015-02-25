/////////////////////////////////////////////////////////////////////////////
/// \file AnmTargetUpdater.h
///
/// \brief It updates the target coordinates of the nodes
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 13/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMTARGETUPDATER_H_
#define ANMTARGETUPDATER_H_

#include <vector>
#include <string>

class AnmParameters;
class AnmEigen;
class AnmMagnitudeCalculator;
class AnmMoveVectorCalculator;
class Atom;
class AtomSet;
class AnmNodeList;

/////////////////////////////////////////////////////////////////////////////
/// \brief It updates the target coordinates of the nodes
///
/// \author mrivero
/// \date 13/09/2012
/////////////////////////////////////////////////////////////////////////////
class AnmTargetUpdater{
	public:
		AnmTargetUpdater(AnmMagnitudeCalculator * magnitudeCalculator,
	 	 	 	 	 	 AnmMoveVectorCalculator * moveVectorCalculator);

		virtual ~AnmTargetUpdater();

		virtual void updateTargetCoordinates(std::vector<double> & targetCoordinates, AnmParameters *anmParameters,
									AnmEigen * eigen, const AnmNodeList & nodes,
									unsigned int chosenModeIndex) = 0;


		std::string generateReport(AnmParameters * anmParameters) const;

	protected:
		AnmMagnitudeCalculator * magnitudeCalculator;
		AnmMoveVectorCalculator * moveVectorCalculator;

		void computeDisplacementMagnitude(AnmParameters * anmParameters, AnmEigen * eigen,
												unsigned int chosenModeIndex);
};

#endif /* ANMTARGETUPDATER_H_ */
