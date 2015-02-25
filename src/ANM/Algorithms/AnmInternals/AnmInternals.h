/////////////////////////////////////////////////////////////////////////////
/// \file AnmInternals.h
///
/// \brief ANM algorithm using internal coordinates
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \author vgil
/// \date 31/07/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMINTERNALS_H_
#define ANMINTERNALS_H_

#include "../../AnmAlgorithm.h"

class ModesCalculator;
class AnmInternalsTargetUpdater;

/////////////////////////////////////////////////////////////////////////////
/// \brief ANM algorithm using internal coordinates
///
/// For information about the algorithms in this class, check the following references:
///Anisotropy of Fluctuation Dynamics of Proteins with an Elastic Network Model, Atilgan et al 2001, Biophys J. 2001 Jan; 80(1)
///
///Normal Mode Analysis (Mathematical and Computational Biology) (Qiang Cui (Editor), Ivet Bahar)
///
///
/// \author mrivero
/// \author vgil
/// \date 31/07/2013
/////////////////////////////////////////////////////////////////////////////
class AnmInternals : public AnmAlgorithm
{
	public:
		AnmInternals(RandomGenerator * randomGenerator, ModesCalculator * modesCalculator,
						AnmInternalsTargetUpdater * anmTargetUpdater, Picker * picker);

		virtual ~AnmInternals();

		void constrainCurrentAnmNodes(EnergyCalculator * enerCalc,
													AnmParameters * anmParameters,
													const AnmNodeList & nodeList);

		void performMovement(EnergyCalculator * enerCalc,
				AnmParameters * anmParameters,
				AtomSet * workingAtomSet,
				const AnmNodeList & nodeList,
				std::vector<double> targetCoords);

		std::string generateReport(AnmParameters * anmParameters) const;

		virtual void calculateModes(AnmParameters * anmParameters,
									 const AnmNodeList & node_list, AnmEigen * eigen);

		virtual void logVectorAsMode(std::string name, std::vector<double>& vector_mode,const AnmNodeList* nodeList){};

	private:
		// Atributes
		AnmInternals(const AnmInternals&);
		AnmInternals& operator=(const AnmInternals&);

		void normalizeEigenVectors(AnmEigen * eigen);


};

#endif /* ANMINTERNALS_H_ */
