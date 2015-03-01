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
/// \author mrivero
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////
class InternalModesCalculator : public ModesCalculator
{
	public:
		InternalModesCalculator();
		virtual ~InternalModesCalculator();

		void calculateEigenValuesAndVectors(AnmParameters * anmParameters,
											const AnmNodeList & node_list, AnmEigen * eigen);
		
		std::string generateReport() const;
		
		static AnmEigen* internalToCartesian(std::vector<Unit*>& units, AnmEigen* in, bool onlyHeavy);
		AnmEigen* cartesianToInternal(std::vector<Unit*>& units, AnmEigen* in);

	private:
		static void calculate_modes(AnmParameters * anmParameters, TriangularMatrix* H,
				 TriangularMatrix* K, AnmEigen* eigen);

	friend class TestICModesCalculator;
};

#endif /* INTERNALMODESCALCULATOR_H_ */
