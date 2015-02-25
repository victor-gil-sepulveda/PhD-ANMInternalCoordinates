/////////////////////////////////////////////////////////////////////////////
/// AnmInternalsTargetUpdater.cpp
///
/// Implementation of AnmInternalsTargetUpdater class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author vgil
/// \author jestrada
/// \date 28/07/2014
/////////////////////////////////////////////////////////////////////////////

#include "AnmInternalsTargetUpdater.h"

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
#include "../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../../ANM/AnmUnitNodeList.h"
#include "../../Molecules/ForceFieldTerms/Dihedral.h"
#include "../../Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/HarmonicDihedralConstraintFunctions.h"
#include <cmath>
#include <iostream>
#include "../../Tools/Math/MathTools.h"
#include "../../System/SystemVars.h"
#include "../../System/Logs/OutputWriter.h"

using namespace std;

AnmInternalsTargetUpdater::AnmInternalsTargetUpdater(AnmMagnitudeCalculator * magnitudeCalculator,
		 	 	 	 	 	 	 	 AnmMoveVectorCalculator * moveVectorCalculator)
  : AnmTargetUpdater(magnitudeCalculator, moveVectorCalculator)
{}

AnmInternalsTargetUpdater::~AnmInternalsTargetUpdater(){}



void AnmInternalsTargetUpdater::updateTargetCoordinates(std::vector<double> & targetCoordinates,
													AnmParameters *anmParameters,
													AnmEigen * eigen, const AnmNodeList & nodes,
													unsigned int chosenModeIndex)
{

	const AnmUnitNodeList& unitNodeList = dynamic_cast<const AnmUnitNodeList &>(nodes);

	// TODO: See what to do with this guy. It overwrites the "displacementFactor" attribute
	// to create the displacementMagnitude ANMParams attribute.
	// In our case is better to have
	anmParameters->setDisplacementMagnitude(anmParameters->getDisplacement());

	vector<double> moveVector;

	moveVectorCalculator->computeMoveVector(moveVector, eigen, chosenModeIndex, anmParameters);

	updateTargetCoordinates(targetCoordinates, moveVector, unitNodeList.getNodeList());

}

void AnmInternalsTargetUpdater::updateTargetCoordinates(std::vector<double> & targetCoordinates,
								const std::vector<double> & moveVector,
								const std::vector<Unit*> & units) {
	targetCoordinates.clear();

	cout << "DBG: Units size: " << units.size() << endl;

	ostringstream oss;
	oss <<"Move vector: ";
	for (unsigned int i = 0; i< moveVector.size(); ++i){
		oss <<moveVector[i]<<" ";
	}
	SystemVars::getLog("move_vector")->write(oss.str());

	ostringstream oss2;
	oss2 <<"Current angles: ";
	for(unsigned int i=0; i<units.size()-1; ++i)
	{
		Dihedral* right = units[i]->getRightDihedral();
		std::vector<Atom*> the_atoms_of_the_dihedral = right->getAtoms();
		double current_angle = HarmonicDihedralConstraintFunctions::calculateDihedralAngleWithArcTanFunction(
				the_atoms_of_the_dihedral[0]->getX(), the_atoms_of_the_dihedral[0]->getY(), the_atoms_of_the_dihedral[0]->getZ(),
				the_atoms_of_the_dihedral[1]->getX(), the_atoms_of_the_dihedral[1]->getY(), the_atoms_of_the_dihedral[1]->getZ(),
				the_atoms_of_the_dihedral[2]->getX(), the_atoms_of_the_dihedral[2]->getY(), the_atoms_of_the_dihedral[2]->getZ(),
				the_atoms_of_the_dihedral[3]->getX(), the_atoms_of_the_dihedral[3]->getY(), the_atoms_of_the_dihedral[3]->getZ());

		oss2 <<Math::degToRad(current_angle)<<" ";

/*************************************************
		// Energies need eq. angles in radians. calc dihed. angle returns degrees (may be good to change to rads).
		double target_angle = moveVector[i] + Math::degToRad(current_angle);

		target_angle = HarmonicDihedralConstraintFunctions::put_in_pi_minus_pi_range(target_angle);

		targetCoordinates.push_back(target_angle);
**************/
	}
	//Copiamos el vector en este caso
	for(unsigned int i=0; i<moveVector.size(); ++i)
		targetCoordinates.push_back(moveVector[i]);

	//--------------------------------
	SystemVars::getLog("move_vector")->write(oss2.str());

	ostringstream oss3;
	oss3 <<"Target Coordinates: ";
	for (unsigned int i = 0; i< targetCoordinates.size(); ++i){
		oss3 <<targetCoordinates[i]<<" ";
	}
	SystemVars::getLog("move_vector")->write(oss3.str());
}
