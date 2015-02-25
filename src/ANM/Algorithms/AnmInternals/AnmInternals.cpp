/////////////////////////////////////////////////////////////////////////////
/// AnmInternals.cpp
///
/// Implementation of AnmInternals class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \author vgil
/// \date 31/07/2013
/////////////////////////////////////////////////////////////////////////////

#include "AnmInternals.h"

#include <sstream>
#include "../../../Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintTerm.h"
#include "../../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../../../Energy/EnergyCalculator.h"
#include "../../../Molecules/ForceFieldTerms/Dihedral.h"
#include "../../AnmUnitNodeList.h"
#include "../../ModesCalculator/AnmEigen.h"
#include "../../TargetUpdater/AnmInternalsTargetUpdater.h"
#include "../../ModesCalculator/Internals/InternalModesCalculator.h"
#include "../../Parameters/AnmParameters.h"
#include "../../Tools/Inout/ModesWriter.h"
#include "../../../Molecules/Atom.h"
#include "DualQuaternion/point3.hpp"
#include "DualQuaternion/dual_quat_cu.hpp"
#include "DualQuaternion/vec3.hpp"
#include "DualQuaternion/quat_cu.hpp"
#include "../../../Tools/Math/Point.h"
#include "../../../Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/HarmonicDihedralConstraintFunctions/HarmonicDihedralConstraintFunctions.h"
#include "../../ModesCalculator/ModesCalculator.h"
#include "../../../Energy/EnergyValues/EnergyValues.h"
#include "../../../MonteCarlo/MetropolisAcceptanceCriterion.h"
#include "../../../Tools/Math/RandomGenerator/RandomGenerator.h"
#include <ctime>
#include "../../../Tools/Math/MathTools.h"
#include "Movement/ANMICMovement.h"
#include "Movement/ANMMover.h"
#include "Movement/SimpleMover/ANMICSimpleMover.h"
#include "Movement/SimpleSDMIN/ANMICSideMinMover.h"
#include "../../../System/SystemVars.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// Constructor
///
/// \param randomGenerator [In] Collaborator that generates random numbers
/// \param modesCalculator [In] Collaborator that computes the eigen values and vectors
/// \param anmTargetUpdater [In] Collaborator that computes the target coordinates
/// \param picker [In] Collaborator that picks a mode
///
/// \author mrivero
/// \date 31/07/2013
///////////////////////////////////////////////////////////////
AnmInternals::AnmInternals(RandomGenerator * randomGenerator, ModesCalculator * modesCalculator,
								AnmInternalsTargetUpdater * anmTargetUpdater, Picker * picker)
 : AnmAlgorithm(randomGenerator, modesCalculator, anmTargetUpdater, picker)
{
}

AnmInternals::~AnmInternals()
{
}

/////////////////////////////////////////////////////////////////////////////////////
/// \remarks
/// It constrains the current dihedrals
///
/// \param complex [In] AtomSet object
/// \param enerCalc [In] EnergyCalculator object
/// \param nodesList [In] List of nodes
/// \param anmParameters [In] ANM parameters
///
/// \author vgil
/// \date 31/07/2013
///////////////////////////////////////////////////////////////////////////////////
void AnmInternals::constrainCurrentAnmNodes(EnergyCalculator * enerCalc,
											AnmParameters * anmParameters,
											const AnmNodeList & nodeList)
{

	// We have to constraint the dihedral angles
	double springConstant = anmParameters->getSteeringForce();
	const AnmUnitNodeList& unitNodeList = dynamic_cast<const AnmUnitNodeList&>(nodeList);
	std::vector<Unit*> units = unitNodeList.getNodeList();

	for(unsigned int i=0; i < units.size()-1; ++i){
		Dihedral* right = units[i]->getRightDihedral();
		std::vector<Atom*> the_atoms_of_the_dihedral = right->getAtoms();
		double current_angle = HarmonicDihedralConstraintFunctions::calculateDihedralAngleWithArcTanFunction(
				the_atoms_of_the_dihedral[0]->getX(), the_atoms_of_the_dihedral[0]->getY(), the_atoms_of_the_dihedral[0]->getZ(),
				the_atoms_of_the_dihedral[1]->getX(), the_atoms_of_the_dihedral[1]->getY(), the_atoms_of_the_dihedral[1]->getZ(),
				the_atoms_of_the_dihedral[2]->getX(), the_atoms_of_the_dihedral[2]->getY(), the_atoms_of_the_dihedral[2]->getZ(),
				the_atoms_of_the_dihedral[3]->getX(), the_atoms_of_the_dihedral[3]->getY(), the_atoms_of_the_dihedral[3]->getZ());

		// Adding a constraint for each dihedral (torsion)
		current_angle = HarmonicDihedralConstraintFunctions::put_in_pi_minus_pi_range(Math::degToRad(current_angle));
		ConstraintTerm * term = new HarmonicDihedralConstraintTerm(
				springConstant,
				current_angle,  // Current angle is the eq. angle, as we want to 'fix' it but let it oscillate
								// a bit if needed
				the_atoms_of_the_dihedral);

		enerCalc->addAnmConstraintTerm(term);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
/// \remarks
/// Does the actual movement of the protein using quaternions to move the torsion angles.
///
/// \param enerCalc [In] EnergyCalculator object
/// \param minParameters [In] Parameters for the minimization
/// \param minimizer [In] The minimizer
///
/// \author vgil
/// \date 15/01/2015
///////////////////////////////////////////////////////////////////////////////////
void AnmInternals::performMovement(EnergyCalculator * enerCalc, 
									AnmParameters * anmParameters,
    								AtomSet * workingAtomSet, 
									const AnmNodeList & nodeList,
									std::vector<double> targetCoords)
{
	const AnmUnitNodeList& unitNodeList = dynamic_cast<const AnmUnitNodeList&>(nodeList);
	std::vector<Unit*> units = unitNodeList.getNodeList();

	// Create a ANM Mover (TODO: This can be moved to a factory when it gets more logic)
	ANMICMover* mover = NULL;
	cout<<"DBG Mover type: "<<anmParameters->getICMoverType()<<endl;
	switch(anmParameters->getICMoverType()){
		case SIMPLE:
			mover = new ANMICSimpleMover(enerCalc, anmParameters);
			break;
		case SIMPLE_PLUS_SIDECHAIN_MIN:
			mover = new ANMICSideMinMover(enerCalc, anmParameters, workingAtomSet);
			break;
		case SIMPLE_PLUS_SIDECHAIN_MD:
			break;
	}

	cout <<"DBG: Mover "<<mover<<endl;

	//Perform the movement
	mover->perform_one_movement_cycle(units, targetCoords);
	cout <<"DBG: Movement cycle performed"<<endl;

//	unsigned int NUMBER_OF_STEPS = 10;
//	temp = 2000
	delete mover;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It generates the ANM internals algorithm report
///
/// \param anmParameters [In] ANM parameters
///
/// \return The ANM internals algorithm report
///
/// \author mrivero
/// \date 31/07/2013
///////////////////////////////////////////////////////////////
std::string AnmInternals::generateReport(AnmParameters * anmParameters) const
{
	ostringstream oss(ostringstream::out);

	oss<<"ANM Internals algorithm:"<<endl;
	oss<<AnmAlgorithm::generateReport(anmParameters)<<endl;

	return oss.str();
}

void AnmInternals::calculateModes(AnmParameters * anmParameters, const AnmNodeList & nodeList, AnmEigen * eigen)
{

	modesCalculator->calculateEigenValuesAndVectors(anmParameters, nodeList, eigen);
	SystemVars::getModesWriterHandler()->getWriter("original_modes_cc")->
			writeInternalModes(eigen, &nodeList, true);
	SystemVars::getModesWriterHandler()->getWriter("original_modes_ic")->
				writeInternalModes(eigen, &nodeList, false);

	normalizeEigenVectors(eigen);
	SystemVars::getModesWriterHandler()->getWriter("normalized_modes_cc")->
					writeInternalModes(eigen, &nodeList, true);
	SystemVars::getModesWriterHandler()->getWriter("normalized_modes_ic")->
						writeInternalModes(eigen, &nodeList, false);

}

void AnmInternals::normalizeEigenVectors(AnmEigen * eigen)
{
	eigen->normalizeByLargestValue();
}
