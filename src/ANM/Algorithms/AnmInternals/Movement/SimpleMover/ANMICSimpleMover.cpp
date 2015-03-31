/////////////////////////////////////////////////////////////////////////////
///
/// \author vgil
/// \date 04/02/2015
/////////////////////////////////////////////////////////////////////////////
#include "ANMICSimpleMover.h"
#include "../../../../../Tools/Math/RandomGenerator/RandomGenerator.h"
#include "../../../../../MonteCarlo/MetropolisAcceptanceCriterion.h"
#include <iostream>
#include "../../../../../Energy/EnergyCalculator.h"
#include "../../../../Parameters/AnmParameters.h"
#include "../ANMICMovement.h"
#include "../../../../../System/SystemVars.h"
#include "../../../../../System/Logs/ModesAndCoords/ModeWritersHandler.h"
using namespace std;

ANMICSimpleMover::ANMICSimpleMover(EnergyCalculator * enerCalc, AnmParameters * anmParameters) :
		ANMICMover(enerCalc, anmParameters){
}

ANMICSimpleMover::~ANMICSimpleMover(){
}

void ANMICSimpleMover::perform_one_movement_cycle(std::vector<Unit*>& units,
		std::vector<double>& targetCoords){
	RandomGenerator random_gen(time(NULL));
	MetropolisAcceptanceCriterion m_criterion(&random_gen);

	double initial_energy = this->enerCalc->calcEnergyWithOnlyShortRangeNonBonding();
	double new_energy;

	bool accepted = true;
	unsigned int m;
	for( m = 0; m < this->anmParameters->getICMaximumStepsPerCycle() && accepted; m++){

		// Apply the rotations
		ANMICMovement::apply_rotations_to_molecule_units( units,
				targetCoords,
				1./anmParameters->getICMaximumStepsPerCycle());

		// Do a metropolis test (only LJ and electrostatic)
		new_energy = enerCalc->calcEnergyWithOnlyShortRangeNonBonding();
		accepted = m_criterion.acceptsThisStep(new_energy,
									initial_energy,
									anmParameters->getICCycleCheckTemp());

		cout<<"DBG: ANM IC Metropolis "<<m_criterion.generateReport()<<endl;
	}

	cout<<"DBG: "<<m<<" of "<<this->anmParameters->getICMaximumStepsPerCycle() <<" iterations performed."<<endl;
	vector<double > tmp(1, anmParameters->getDisplacement()*(m/anmParameters->getICMaximumStepsPerCycle()));
	SystemVars::getModesWriterHandler()->logStepAndVector("max_displacement", tmp);
}
