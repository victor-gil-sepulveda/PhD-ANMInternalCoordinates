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
#include <stdlib.h>
#include <cmath>
#include "../../../Tools/stringTools.h"
#include "../../../System/Logs/OutputWriter.h"
#include "../../Tools/AnmNormalizer.h"
#include "../../../Tools/Utils.h"
#include <algorithm>
#include "../../../Molecules/AtomSet/AtomSet.h"

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
 : AnmAlgorithm(randomGenerator, modesCalculator, anmTargetUpdater, picker){
}

AnmInternals::~AnmInternals(){
}

void AnmInternals::calculate_current_angles(vector<double>& current_angles, vector<Unit*>& units){

	for(unsigned int i=0; i < units.size()-1; ++i){
		Dihedral* right = units[i]->right_dihedral;
		std::vector<Atom*> the_atoms_of_the_dihedral = right->getAtoms();
		double current_angle = HarmonicDihedralConstraintFunctions::calculateDihedralAngleWithArcTanFunction(
				the_atoms_of_the_dihedral[0]->getX(), the_atoms_of_the_dihedral[0]->getY(), the_atoms_of_the_dihedral[0]->getZ(),
				the_atoms_of_the_dihedral[1]->getX(), the_atoms_of_the_dihedral[1]->getY(), the_atoms_of_the_dihedral[1]->getZ(),
				the_atoms_of_the_dihedral[2]->getX(), the_atoms_of_the_dihedral[2]->getY(), the_atoms_of_the_dihedral[2]->getZ(),
				the_atoms_of_the_dihedral[3]->getX(), the_atoms_of_the_dihedral[3]->getY(), the_atoms_of_the_dihedral[3]->getZ());
		current_angle = HarmonicDihedralConstraintFunctions::put_in_pi_minus_pi_range(Math::degToRad(current_angle));
		current_angles.push_back(current_angle);
	}
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

	vector<double> current_angles;
	calculate_current_angles(current_angles, units);

	for(unsigned int i=0; i < units.size()-1; ++i){
		Dihedral* right = units[i]->right_dihedral;
		std::vector<Atom*> the_atoms_of_the_dihedral = right->getAtoms();
		// Adding a constraint for each dihedral (torsion)
		ConstraintTerm * term = new HarmonicDihedralConstraintTerm(
				springConstant,
				current_angles[i],  // Current angle is the eq. angle, as we want to 'fix' it but let it oscillate
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
	cout<<"DBG: Anm IC Mover type: "<<anmParameters->getICMoverType()<<endl;
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

	// Perform the movement
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
	// Calculates frequencies and sets frequency counters to a random offset;
	modesCalculator->calculateEigenValuesAndVectors(anmParameters, nodeList, eigen);

	if(anmParameters->getOverridePickingAndMixing() != DO_NOT_OVERRIDE)
		resetFrequencies(anmParameters, eigen);

	// Logging
	SystemVars::getModesWriterHandler()->getWriter("original_modes_cc")->
			writeInternalModes(eigen, &nodeList, true);
	SystemVars::getModesWriterHandler()->getWriter("original_modes_ic")->
				writeInternalModes(eigen, &nodeList, false);

	// Normalization
	normalizeEigenVectors(eigen);

	// Logging
	SystemVars::getModesWriterHandler()->getWriter("normalized_modes_cc")->
					writeInternalModes(eigen, &nodeList, true);
	SystemVars::getModesWriterHandler()->getWriter("normalized_modes_ic")->
						writeInternalModes(eigen, &nodeList, false);

}

void AnmInternals::resetFrequencies(AnmParameters * anmParameters, AnmEigen * eigen){

	cout<<"DBG: Resetting mode counters "<<endl;
	unsigned int max_cycle = anmParameters->getMajorModePeriod()*2;

	// Prepare cycle frequencies
	mode_frequency.clear();

	// Calculate real frequencies
	vector<double> real_mode_frequencies;
	for (unsigned int i =0; i < eigen->values.size(); ++i){
		if(anmParameters->getMainModeWeightForMixModes() < 1){
			real_mode_frequencies.push_back(1./eigen->values[i]); // TODO: Use real formula
		}
		else{
			real_mode_frequencies.push_back(eigen->values[i]);
		}
		//cout<<"KK "<<eigen->values[i]<<endl;
	}
	// Divide by the maximum value, as eigenvalues are ordered, the maximum is the first one :)
	double max_real_cycle =  real_mode_frequencies[0];
	for (unsigned int i =0; i < eigen->values.size(); ++i){
		//cout<<"LL "<<real_mode_frequencies[i]<<endl;
		real_mode_frequencies[i] = real_mode_frequencies[i] /max_real_cycle;
	}
	// Assign the cycle frequencies
	for (unsigned int i =0; i < eigen->values.size(); ++i){
		mode_frequency.push_back(max((int)(floor(max_cycle*real_mode_frequencies[i])),1)); //Half of the cycle in one direction, half in the other
		//cout<<"XX "<<max_cycle<<" "<<max((int)(floor(max_cycle*real_mode_frequencies[i])),2)<<endl;
	}

	// Randomize initial counters
	mode_frequency_counter.clear();
	srand((unsigned)time(NULL));
	for (unsigned int i =0; i < mode_frequency.size(); ++i){
		mode_frequency_counter.push_back(((double)rand()/RAND_MAX)*mode_frequency[i]);
	}

	// Log the periods
	SystemVars::getModesWriterHandler()->logStepAndVector("mode_counters", mode_frequency);
}

void AnmInternals::calculateTargetCoords(AnmParameters * anmParameters, AnmEigen * eigen,
						const AnmNodeList & nodeList, std::vector<double> & targetCoords,
						unsigned int chosenMode){

	vector<Unit*> units = dynamic_cast<AnmUnitNodeList&>(const_cast<AnmNodeList&>(nodeList)).getNodeList();

	// Log current coordinates and proposal
	cout<<"DBG: Logging initial coordinates and conformation proposal. "<<endl;
	SystemVars::getModesWriterHandler()->logStepAndCurrentCACoordinates("initial_cc");
	SystemVars::getModesWriterHandler()->logStepAndCurrentDihedralAngles("initial_ic", units);

	// If the mixing case is the eigenvalues one, we intercept the whole calculation protocol
	if(anmParameters->getOverridePickingAndMixing() == DO_NOT_OVERRIDE){
		cout<<"DBG: Calculating target coords using Regular Control File Parameters"<<endl;
		anmTargetUpdater->updateTargetCoordinates(targetCoords, anmParameters, eigen, nodeList,
				chosenMode);
	}
	else{
		if(anmParameters->getOverridePickingAndMixing() == EIGEN_MIXING_EIGEN_WEIGHT){
			cout<<"DBG: Calculating target coords using EIGEN_MIXING_EIGEN_WEIGHT"<<endl;
			unsigned int eigensize = eigen->vectors[0].size();

			// Initialize target coordinates
			targetCoords.clear();
			targetCoords.resize(eigensize,0);
			for (unsigned int j = 0; j < eigensize; ++j){
				targetCoords[j] = 0;
			}

			// Ensure modes are normalized (norm of the mode must be one)
			for(unsigned int i = 0; i < eigen->vectors.size(); ++i){
				vector<double> & eigenvector = eigen->vectors[i];
				Math::normalizeVector(Utils::vectorToPointer(eigenvector), eigensize);
			}

			// Calculate target angle increments
			for (unsigned int i = 0; i < eigen->vectors.size(); ++i){
				double sense = mode_frequency_counter[i] > mode_frequency[i]/2? 1:-1;
				double atenuation = eigen->values[i] / eigen->values[0];
				for (unsigned int j = 0; j < eigensize; ++j){
					targetCoords[j] += sense*atenuation*eigen->vectors[i][j];
				}
			}

			// Put all angles in the correct range
			for (unsigned int i = 0; i < targetCoords.size(); ++i){
				targetCoords[i] = HarmonicDihedralConstraintFunctions::put_in_pi_minus_pi_range(targetCoords[i]);
			}

			// And normalize the result!
			AnmNormalizer::normalizeByLargestValue(targetCoords);

			// Then multiply the maximum displacement
			Math::multiplyVectorByScalar(targetCoords, anmParameters->getDisplacement());

			// Log frequency signs
			vector<double> signs;
			for (unsigned int i =0; i < mode_frequency.size(); ++i){
				signs.push_back((mode_frequency_counter[i] > mode_frequency[i]/2? 1:-1));
			}
			SystemVars::getModesWriterHandler()->logStepAndVector("mode_counters", mode_frequency_counter);
			SystemVars::getModesWriterHandler()->logStepAndVector("mode_signs", signs);

			// Update frequency counters
			for (unsigned int i =0; i < mode_frequency.size(); ++i){
				mode_frequency_counter[i] = (mode_frequency_counter[i]+1)%mode_frequency[i];
			}

		}
		else{
			cout<<"DBG:  ERROR, Picking and mixing override type has not been defined."<<endl;
			exit(-1);
		}
	}

	SystemVars::getModesWriterHandler()->logStepAndVector("proposal_inc_ic",targetCoords);
	SystemVars::getModesWriterHandler()->logStepAndDihedralAnglesProposal("proposal_ic",targetCoords,units);
	SystemVars::getModesWriterHandler()->logStepAndDihedralToCartesianProposal("proposal_cc", targetCoords, units);
}

void AnmInternals::normalizeEigenVectors(AnmEigen * eigen){
	eigen->normalizeByLargestValue();
}

void AnmInternals::logAfterANMCoords(AnmNodeList* nodeList){
	vector<Unit*> units = dynamic_cast<AnmUnitNodeList*>(const_cast<AnmNodeList*>(nodeList))->getNodeList();
	cout << "DBG: Logging coordinates after anm."<<endl;
	SystemVars::getModesWriterHandler()->logStepAndCurrentCACoordinates("after_anm_cc");
	SystemVars::getModesWriterHandler()->logStepAndCurrentDihedralAngles("after_anm_ic",units);
}

void AnmInternals::logAfterMinimizationCoords(AnmNodeList* nodeList){
	vector<Unit*> units = dynamic_cast<AnmUnitNodeList*>(const_cast<AnmNodeList*>(nodeList))->getNodeList();
	cout << "DBG: Logging final step coordinates."<<endl;
	SystemVars::getModesWriterHandler()->logStepAndCurrentCACoordinates("after_minimization_cc");
	SystemVars::getModesWriterHandler()->logStepAndCurrentDihedralAngles("after_minimization_ic",units);
}
