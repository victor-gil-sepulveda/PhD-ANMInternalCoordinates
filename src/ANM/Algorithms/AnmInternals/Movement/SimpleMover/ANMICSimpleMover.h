/////////////////////////////////////////////////////////////////////////////
///
/// \author vgil
/// \date 04/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ICSIMPLEMOVER_H_
#define ICSIMPLEMOVER_H_

#include "../ANMMover.h"
class ANMICSimpleMover : public ANMICMover{

	public:
		ANMICSimpleMover(EnergyCalculator * enerCalc,
				AnmParameters * anmParameters);

		virtual ~ANMICSimpleMover();

		void perform_one_movement_cycle(std::vector<Unit*>& units, std::vector<double>& targetCoords);


};

#endif /* ICSIMPLEMOVER_H_ */
