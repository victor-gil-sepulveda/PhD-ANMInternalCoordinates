/////////////////////////////////////////////////////////////////////////////
/// AnmAlgorithm.cpp
///
/// Implementation of AnmAlgorithm class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author atarraco
/// \author mrivero
/// \author vgil
/// \date 31/08/2012
/////////////////////////////////////////////////////////////////////////////

#include "AnmAlgorithm.h"
#include "ModesCalculator/ModesCalculator.h"
#include "TargetUpdater/AnmTargetUpdater.h"
#include "Pickers/Picker.h"
#include "../Tools/Utils.h"
#include <sstream>
#include <ostream>
#include "../Molecules/AtomSet/AtomSet.h"
#include "Parameters/AnmParameters.h"
#include "ModesCalculator/AnmEigen.h"
#include "Tools/Inout/ModesWriter.h"
#include "../System/SystemVars.h"

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
/// \author atarraco
/// \date 03/09/2012
///////////////////////////////////////////////////////////////
AnmAlgorithm::AnmAlgorithm(RandomGenerator * randomGenerator, ModesCalculator * modesCalculator,
							AnmTargetUpdater * anmTargetUpdater, Picker * picker)
{
	this->randomGenerator = randomGenerator;
	this->modesCalculator = modesCalculator;
	this->anmTargetUpdater = anmTargetUpdater;
	this->picker = picker;
}

AnmAlgorithm::~AnmAlgorithm()
{
	Utils::deleteAndSetPointerToNull<ModesCalculator>(modesCalculator);
	Utils::deleteAndSetPointerToNull<AnmTargetUpdater>(anmTargetUpdater);
	Utils::deleteAndSetPointerToNull<Picker>(picker);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function picks a new mode
///
/// \param eigen [In] Eigen vectors and eigen values
/// \param anmParameters [In] ANM parameters
///
/// \return The index of the chosen mode
///
/// \author mrivero
/// \date 02/10/2012
///////////////////////////////////////////////////////////////
unsigned int AnmAlgorithm::pickNewMode(AnmEigen * eigen, AnmParameters * anmParameters)
{
	return picker->selectMode(eigen, anmParameters);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It generates the Anm algorithm report
///
/// \param anmParameters [In] ANM parameters
///
/// \return The Anm algorithm report
///
/// \author mrivero
/// \date 23/04/2013
///////////////////////////////////////////////////////////////
std::string AnmAlgorithm::generateReport(AnmParameters * anmParameters) const
{
	ostringstream oss(ostringstream::out);

	oss<<anmTargetUpdater->generateReport(anmParameters)<<endl;
	oss<<"Normal modes: \n"<< modesCalculator->generateReport() <<endl;
	oss<<"Modes picking strategy:"<<endl;
	oss<<picker->generateReport();

	return oss.str();
}



///////////////////////////////////////////////////////////////
/// \remarks
/// Normalize the eigenvectors
///
/// \param useThermalScaling [In] If true thermal scaling is used
/// \param constantForHessian [In] Constant used in normalization with thermal scaling
///
/// \param eigen [In/Out] Eigenvalues and eigenVectors
///
/// \author atarraco
/// \author mrivero
/// \date 05/07/2012
///////////////////////////////////////////////////////////////
void AnmAlgorithm::normalizeEigenVectors(AnmEigen * eigen,
															bool useThermalScaling,
															double constantForHessian)
{
	if(useThermalScaling) {
		eigen->normalizeByInverseNormWithThermal(constantForHessian);
	}
	else {
		eigen->normalizeByInverseLargestNorm();
	}
}

