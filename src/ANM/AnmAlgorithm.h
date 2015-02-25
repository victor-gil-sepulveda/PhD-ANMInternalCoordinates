///////////////////////////////////////////////////////////
/// \file AnmAlgorithm.h
///
/// \brief Common interface for all ANM algorithm classes.
/// It's an abstract class. It includes the methods that
/// are common to all ANM algorithm.
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author atarraco
/// \author mrivero
/// \author vgil
/// \date  06/05/2012
///////////////////////////////////////////////////////////

#pragma once

#ifndef ANMALGORITHM_H_
#define ANMALGORITHM_H_

#include <vector>
#include <string>


class Atom;
class AtomSet;
class AnmParameters;
class RandomGenerator;
class MinParams;
class EnergyCalculator;
class Minimizer;
class AnmEigen;
class ModesCalculator;
class AnmTargetUpdater;
class Picker;
class AnmNodeList;
class ModesWriter;

///////////////////////////////////////////////////////////
/// \brief Common interface for all ANM algorithm classes.
/// It's an abstract class. It includes the methods that
/// are common to all ANM algorithm.
///
/// \author atarraco
/// \author mrivero
/// \author vgil
/// \date  06/05/2012
///////////////////////////////////////////////////////////
class AnmAlgorithm
{
	public:
		AnmAlgorithm(RandomGenerator * randomGenerator, ModesCalculator * modesCalculator,
						AnmTargetUpdater * anmTargetUpdater, Picker * picker);
		virtual ~AnmAlgorithm();

		virtual void calculateModes(AnmParameters * anmParameters,
									 const  AnmNodeList & node_list, AnmEigen * eigen) = 0;

		virtual void calculateTargetCoords(AnmParameters * anmParameters, AnmEigen * eigen,
											const AnmNodeList & nodeList, std::vector<double> & targetCoords,
											unsigned int chosenMode);

		virtual unsigned int pickNewMode(AnmEigen * eigen, AnmParameters * anmParameters);

		virtual void constrainCurrentAnmNodes(EnergyCalculator * enerCalc,
													AnmParameters * anmParameters,
													const AnmNodeList & nodeList) = 0;

		virtual void performMovement( EnergyCalculator * enerCalc,
				AnmParameters * anmParameters,
				AtomSet * workingAtomSet,
				const AnmNodeList & nodeList,
				std::vector<double> targetCoords) = 0;

		virtual std::string generateReport(AnmParameters * anmParameters) const;

		// Mode logging
		virtual void logVectorAsMode(std::string name, std::vector<double>& vector_mode,const AnmNodeList* nodeList) = 0;

	protected:
		// Atributes
		RandomGenerator * randomGenerator;
		ModesCalculator * modesCalculator;
		AnmTargetUpdater * anmTargetUpdater;
		Picker * picker;

		void normalizeEigenVectors(AnmEigen * eigen, bool useThermalScaling, double constantForHessian);

	private:
		AnmAlgorithm(const AnmAlgorithm&);
		AnmAlgorithm& operator=(const AnmAlgorithm&);

		friend class TestAnm;
		friend class TestAnmTargetUpdater;
};
#endif /* ANMALGORITHM_H_ */
