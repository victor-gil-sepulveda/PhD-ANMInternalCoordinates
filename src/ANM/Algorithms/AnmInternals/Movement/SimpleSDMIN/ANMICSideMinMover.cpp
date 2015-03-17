/////////////////////////////////////////////////////////////////////////////
///
/// \author vgil
/// \date 13/02/2015
/////////////////////////////////////////////////////////////////////////////
#include "ANMICSideMinMover.h"
#include "../ANMICMovement.h"
#include "../../../../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../../../../../Molecules/Selection/AtomsListSelection.h"
#include "../../../../../System/StateUpdater.h"
#include "../../../../../Molecules/Solvent/Solvent.h"
#include "../../../../../Minimizers/MinimizerBuilder.h"
#include "../../../../../Minimizers/Minimizer.h"
#include "../../../../../Minimizers/Parameters/MinParams.h"
#include "../../../../../Molecules/Freezers/Freezer.h"
#include "../../../../../Molecules/AtomSet/AtomSet.h"
#include "../../../../../Tools/Math/RandomGenerator/RandomGenerator.h"
#include "../../../../../MonteCarlo/MetropolisAcceptanceCriterion.h"
#include "../../../../../Energy/EnergyCalculator.h"
#include "../../../../Parameters/AnmParameters.h"
#include "../../../../../Minimizers/Parameters/TNParams.h"

#include <iostream>

using namespace std;

ANMICSideMinMover::ANMICSideMinMover(EnergyCalculator * enerCalc, AnmParameters * anmParameters, AtomSet* atomset)  :
				ANMICMover(enerCalc, anmParameters){
	this->atomset = atomset;
}

ANMICSideMinMover::~ANMICSideMinMover() {
}


void ANMICSideMinMover::perform_one_movement_cycle(vector<Unit*>& units, vector<double>& targetCoords){
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

		// Try to fix sidechains only
		perform_side_minim(this->atomset, this->enerCalc);

		// Do a metropolis test (only LJ and electrostatic)
		new_energy = enerCalc->calcEnergyWithOnlyShortRangeNonBonding();
		accepted = m_criterion.acceptsThisStep(new_energy,
									initial_energy,
									anmParameters->getICCycleCheckTemp());

		cout<<"DBG: ANM IC Metropolis "<<m_criterion.generateReport()<<endl;
	}
	cout<<"DBG: "<<m<<" iterations performed."<<endl;

	// Small free minimization
	MinimizerBuilder minimizerBuilder;;
	Minimizer * min = minimizerBuilder.createMinimizer(TN_MIN);
	MinParams *min_params = new TNParams(1.0, 0.001, 1, 20);
	min->minimize(min_params, this->atomset, this->enerCalc);

	// Free mem
	delete min;
	delete min_params;
}

void ANMICSideMinMover::perform_side_minim(AtomSet* atom_set, EnergyCalculator * energyCalculator) {

	Solvent* solvent = energyCalculator->getSolvent();

	// Frozen
	vector<Atom*> allButExtendedBackBoneAtoms;
	AtomsListSelection::getAllButBackBone(allButExtendedBackBoneAtoms, atom_set);
	atom_set->freeze();
	Freezer::unfreezeAtoms(allButExtendedBackBoneAtoms);
	StateUpdater::updateBecauseOfFrozenChanges(atom_set, solvent, energyCalculator);

	// Minimize
	MinimizerBuilder minimizerBuilder;
	Minimizer * min = minimizerBuilder.createMinimizer(TN_MIN);
	MinParams * min_params = new TNParams(1.0, 0.001, 1, 100, false);
	min->minimize(min_params, atom_set, energyCalculator);

	// Unfrozen
	StateUpdater::fullUpdate(atom_set, solvent, energyCalculator);
	delete min;
	delete min_params;

//	// Small free minimization
//	min = minimizerBuilder.createMinimizer(TN_MIN);
//	min_params = new TNParams(1.0, 0.001, 1, 20);
//	min->minimize(min_params, atom_set, energyCalculator);
//
//	// Free mem
//	delete min;
//	delete min_params;
}

