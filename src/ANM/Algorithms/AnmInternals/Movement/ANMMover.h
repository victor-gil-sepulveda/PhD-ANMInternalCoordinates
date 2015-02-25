/////////////////////////////////////////////////////////////////////////////
/// \author vgil
/// \date 04/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANM_IC_MOVER_H_
#define ANM_IC_MOVER_H_

#include <vector>

class EnergyCalculator;
class AnmParameters;
class Unit;

class ANMICMover {
	public:
		 ANMICMover(EnergyCalculator * enerCalc,
					AnmParameters * anmParameters){
			this->enerCalc = enerCalc;
			this->anmParameters = anmParameters;
			this->targetCoords = targetCoords;
		 }

		virtual ~ANMICMover(){}

		virtual void perform_one_movement_cycle(std::vector<Unit*> & units, std::vector<double>& targetCoords) = 0;

	protected:
		EnergyCalculator * enerCalc;
		AnmParameters * anmParameters;
		std::vector<double> targetCoords;

};
#endif /* ANMMOVER_H_ */
