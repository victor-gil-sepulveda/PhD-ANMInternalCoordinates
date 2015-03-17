/*
 * AnmInternalsTargetUpdater.h
 *
 *  Created on: 25/07/2014
 *      Author: jestrada
 */

#ifndef ANMINTERNALSTARGETUPDATER_H_
#define ANMINTERNALSTARGETUPDATER_H_

#include <vector>
#include <string>

#include "AnmTargetUpdater.h"

class AnmParameters;
class AnmInternalsEigen;
class AnmMagnitudeCalculator;
class AnmInternalsMoveVectorCalculator;
class Atom;
class AtomSet;
class Unit;

class AnmInternalsTargetUpdater : public AnmTargetUpdater
{
	public:
		AnmInternalsTargetUpdater(AnmMagnitudeCalculator * magnitudeCalculator,
	 	 	 	 	 	 AnmMoveVectorCalculator * moveVectorCalculator);

		virtual ~AnmInternalsTargetUpdater();

		virtual void updateTargetCoordinates(std::vector<double> & targetCoordinates, AnmParameters *anmParameters,
									AnmEigen * eigen, const AnmNodeList & nodes,
									unsigned int chosenModeIndex);

	private:

		void updateTargetCoordinates(std::vector<double> & targetCoordinates, const std::vector<double> & moveVector,
									 const std::vector<Unit*> & nodes);



};




#endif /* ANMINTERNALSTARGETUPDATER_H_ */
