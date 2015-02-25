///////////////////////////////////////////////////////////
/// \file AnmCalculator.h
///
/// \brief Computes different ANM algorithms.
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author atarraco
/// \author vgil
/// \author mrivero
/// \author icabeza
/// \date 28/10/2010
///////////////////////////////////////////////////////////

#pragma once

#ifndef ANMCALCULATOR_H_
#define ANMCALCULATOR_H_

#include "Algorithms/AnmAlgorithmTypes.h"
#include<vector>

class AnmParameters;
class AnmAlgorithm;
class AtomSet;
class Atom;
class MinParams;
class Minimizer;
class EnergyCalculator;
class AnmEigen;
class AnmOptions;
class DirectionSelector;
class AnmNodeList;

///////////////////////////////////////////////////////////
/// \brief Computes different ANM algorithms.
///
/// \author atarraco
/// \author vgil
/// \author mrivero
/// \author icabeza
/// \date 28/10/2010
///////////////////////////////////////////////////////////
class AnmCalculator
{
	public:
		AnmCalculator(AnmAlgorithm * algorithm, AnmParameters * anmParameters, AnmOptions * anmOptions,
						AtomSet * atomSet, AnmNodeList *nodeList, DirectionSelector * directionSelector);
		virtual ~AnmCalculator();

		// Methods
		void computeTargetCoordinatesForAnmNodes();

		void performMovement(EnergyCalculator * enerCalc);

		void constrainCurrentAnmNodes(EnergyCalculator * enerCalc);

		std::string listEigenValuesAndEigenVectors();

		std::string generateReport() const;
		std::string showNodeList() const;

		AnmAlgorithm * getAlgorithm();

	private:
		// Attributes
		AnmAlgorithm * algorithm;
		AnmParameters * anmParameters;
		AnmOptions * anmOptions;
		DirectionSelector * directionSelector;

		AnmNodeList *nodeList;
		std::vector<double> targetCoords;
		AnmEigen * eigen;
		unsigned int step;
		AtomSet * workingAtomSet;
		unsigned int chosenModeIndex;

		bool anm_calculated_at_least_once;

		// Methods
		bool isTimeToUpdateModes();
		bool isTimeToPickNewMode();

		AnmCalculator(const AnmCalculator&);
		AnmCalculator& operator=(const AnmCalculator&);

	friend class TestAnm;
	friend class TestPeleComponentsBuilder;
};

#endif /* ANMCALCULATOR_H_ */
