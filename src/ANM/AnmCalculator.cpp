///////////////////////////////////////////////////////////
/// AnmCalculator.cpp
///
/// Implementation of AnmCalculator class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author atarraco
/// \author vgil
/// \author mrivero
/// \author icabeza
/// \date  28/10/2010
///////////////////////////////////////////////////////////

#include "AnmCalculator.h"

#include "../Molecules/AtomSet/AtomSet.h"
#include "../Energy/EnergyCalculator.h"
#include "../Minimizers/Minimizer.h"
#include "Parameters/AnmParameters.h"
#include "ModesCalculator/AnmEigen.h"
#include "../Tools/Utils.h"
#include "AnmAlgorithm.h"
#include "Options/AnmOptions.h"
#include <sstream>
#include <ostream>
#include <stdexcept>

#include "../System/SystemVars.h"
#include "../PELE/PeleTasks/Output/Loggers/EventsLoggers/Types/PeleSimulationLogger.h"
#include "../PELE/PeleTasks/Output/LogUtils.h"
#include "Tools/DirectionSelection/DirectionSelector.h"
#include "AnmNodeList.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// Class constructor.
///
/// \param anmParameters [In] ANM parameters
/// \param algorithm [In] ANM algorithm
/// \param atomSet [In] AtomSet object
/// \param anmOptions [In] ANM options
/// \param nodeList [In] List of nodes
/// \param directionSelector [In] DirectionSelector collaborator
///
/// \author mrivero
/// \author atarraco
/// \author xoro
/// \date 03/09/2012
///////////////////////////////////////////////////////////////
AnmCalculator::AnmCalculator(AnmAlgorithm * algorithm, AnmParameters * anmParameters, AnmOptions * anmOptions,
								AtomSet * atomSet, AnmNodeList *nodeList, DirectionSelector * directionSelector)
{
	this->anmParameters = anmParameters;

	this->anmOptions = anmOptions;

	this->algorithm = algorithm;

	this->eigen = new AnmEigen;

	this->step = 0;

	this->workingAtomSet = atomSet;

	chosenModeIndex = 0;

	this->nodeList = nodeList;

	this->directionSelector = directionSelector;
}

AnmCalculator::~AnmCalculator()
{
	Utils::deleteAndSetPointerToNull<AnmAlgorithm>(algorithm);

	Utils::deleteAndSetPointerToNull<AnmEigen>(eigen);

	Utils::deleteAndSetPointerToNull<AnmOptions>(anmOptions);

	Utils::deleteAndSetPointerToNull<DirectionSelector>(directionSelector);

	Utils::deleteAndSetPointerToNull<AnmNodeList>(nodeList);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Computes the target coordinates of the Nodes
///
/// \author atarraco
/// \author xoro
/// \author mrivero
/// \date 03/09/2012
///////////////////////////////////////////////////////////////
void AnmCalculator::computeTargetCoordinatesForAnmNodes()
{
	cout<<"DBG: Entering AnmCalculator::computeTargetCoordinatesForAnmNodes()"<<endl;

	if(isTimeToUpdateModes()) {
		SystemVars::SIMULATION_LOGGER->logLine("Computing Eigenvalues and Eigenvectors...");
		cout << "Calculate Eigenvalues and EigenVectors..."<<endl;
		algorithm->calculateModes(anmParameters, *nodeList, eigen);
	}

	if(isTimeToPickNewMode()) {
		cout << "Pick new mode..."<<endl;
		SystemVars::SIMULATION_LOGGER->logLine("Pick new mode...");
		chosenModeIndex = algorithm->pickNewMode(eigen, anmParameters);
		directionSelector->updateDirection();
	}

	cout<<"Calculate move vector and new target coords..."<<endl;

	algorithm->calculateTargetCoords(anmParameters, eigen, *nodeList, targetCoords, chosenModeIndex);

	this->step++;
}


/////////////////////////////////////////////////////////////////////////////////////
/// \remarks
/// Does the actual movement of the protein, in a way which depends on the selected algorithm
///
/// \param enerCalc [In] EnergyCalculator object
/// \param minParameters [In] Parameters for the minimization
/// \param minimizer [In] The minimizer
///
/// \author priera
/// \date 15/01/2015
///////////////////////////////////////////////////////////////////////////////////
void AnmCalculator::performMovement(EnergyCalculator * enerCalc){
	cout<<"ANM movement..."<<endl;

	algorithm->performMovement(enerCalc, anmParameters, workingAtomSet, *nodeList, targetCoords);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function updates the constraint list
///
/// \param enerCalc [In] Energy Calculator
///
/// \author mrivero
/// \date 04/09/2012
///////////////////////////////////////////////////////////////
void AnmCalculator::constrainCurrentAnmNodes(EnergyCalculator * enerCalc)
{
	algorithm->constrainCurrentAnmNodes(enerCalc, anmParameters, *nodeList);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It returns a string showing the list of eigenvalues and
/// eigenVectors
///
/// \return A string showing the list of eigenvalues and eigenVectors
///
/// \author mrivero
/// \date 18/09/2012
///////////////////////////////////////////////////////////////
std::string AnmCalculator::listEigenValuesAndEigenVectors()
{
	return eigen->toString();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function checks the modes updating frequency to
/// find out if it's time to update them or not
///
/// \return True it it's, false otherwise
///
/// \author mrivero
/// \date 18/09/2012
///////////////////////////////////////////////////////////////
bool AnmCalculator::isTimeToUpdateModes()
{
	return step % anmParameters->getEigenUpdateFrequency() == 0;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function checks the pick-new-mode frequency to
/// find out if it's time to pick a new one or not
///
/// \return True it it's, false otherwise
///
/// \author mrivero
/// \date 18/09/2012
///////////////////////////////////////////////////////////////
bool AnmCalculator::isTimeToPickNewMode()
{
	return step % anmParameters->getModesChangeFrequency() == 0;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It generates the ANM report
///
/// \return The ANM report
///
/// \author mrivero
/// \date 23/04/2013
///////////////////////////////////////////////////////////////
std::string AnmCalculator::generateReport() const
{
	ostringstream oss(ostringstream::out);
	oss<<"** ANM:"<<endl;
	oss<<"ANM options:"<<endl;
	oss<<this->algorithm->generateReport(anmParameters)<<endl;
	oss<<"ANM parameters:"<<endl;
	oss<<anmParameters->showNumberOfModes()<<endl;
	oss<<anmParameters->showUseThermalScalingInPicking()<<endl;
	oss<<anmParameters->showCutoff()<<endl;
	oss<<anmParameters->showEigenUpdateFrequency()<<endl;
	oss<<anmParameters->showModesChangeFrequency()<<endl;
	oss<<showNodeList()<<endl;

	return oss.str();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It shows the node list
///
/// \return The list of nodes
///
/// \author dlecina
/// \date 02/05/2013
///////////////////////////////////////////////////////////////
std::string AnmCalculator::showNodeList() const {
	ostringstream oss (ostringstream::out);
	oss << "Number of nodes: " << this->nodeList->size() << endl;
	oss << "Node list: " << endl;
	oss << this->nodeList->showNodeList();
	return oss.str();
}

AnmAlgorithm * AnmCalculator::getAlgorithm(){
	return this->algorithm;
}
