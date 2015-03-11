/////////////////////////////////////////////////////////////////////////////
/// \file InternalModesCalculator.h
///
/// \brief It computes modes in internal coordinates
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef INTERNALMODESCALCULATOR_H_
#define INTERNALMODESCALCULATOR_H_
#include "../ModesCalculator.h"

#include "MatrixCalculationFunctions/TriangularMatrices.h"
#include <vector>
#include "CoarseGrainModel/Unit.h"

/////////////////////////////////////////////////////////////////////////////
/// \brief It computes modes in internal coordinates
///
/// \author vgil
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////
class InternalModesCalculator : public ModesCalculator
{
	public:
		InternalModesCalculator();
		~InternalModesCalculator();

		void calculateEigenValuesAndVectors(AnmParameters*, const AnmNodeList& , AnmEigen*);
		
		static AnmEigen* internalToCartesian(std::vector<Unit*>& units, AnmEigen* in, bool onlyHeavy);
		static AnmEigen* cartesianToInternal(std::vector<Unit*>& units, AnmEigen* in);
		static AnmEigen* cartesianToInternalGeometrical(std::vector<Unit*>& units, AnmEigen* in);

		std::string generateReport() const;

	private:
		static void calculate_modes(AnmParameters * anmParameters, TriangularMatrix* H,
				 TriangularMatrix* K, AnmEigen* eigen);


	friend class TestICModesCalculator;
};

#endif /* INTERNALMODESCALCULATOR_H_ */
