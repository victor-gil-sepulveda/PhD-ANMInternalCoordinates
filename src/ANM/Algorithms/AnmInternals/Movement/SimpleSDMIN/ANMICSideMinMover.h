/////////////////////////////////////////////////////////////////////////////
/// \author vgil
/// \date 13/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMICSIDEMINMOVER_H_
#define ANMICSIDEMINMOVER_H_

#include "../ANMMover.h"
#include <vector>

class Unit;
class AtomSet;
class EnergyCalculator;

class ANMICSideMinMover: public ANMICMover {
	public:
		ANMICSideMinMover(EnergyCalculator * enerCalc, AnmParameters * anmParameters, AtomSet* atomset);
		virtual ~ANMICSideMinMover();

		void perform_one_movement_cycle(std::vector<Unit*>& units, std::vector<double>& targetCoords);
		void perform_side_minim(AtomSet* atomset, EnergyCalculator * energyCalculator);

	protected:
		AtomSet* atomset;
};

#endif /* ANMICSIDEMINMOVER_H_ */
