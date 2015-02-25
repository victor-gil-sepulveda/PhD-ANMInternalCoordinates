/////////////////////////////////////////////////////////////////////////////
/// \file ModesCalculator.h
///
/// \brief Interface for different algorithms to compute ANM modes
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef MODESCALCULATOR_H_
#define MODESCALCULATOR_H_

#include <vector>
#include <iosfwd>

class AnmParameters;
class AnmEigen;
class Atom;
class AnmNodeList;

/////////////////////////////////////////////////////////////////////////////
/// \brief Interface for different algorithms to compute ANM modes
///
/// \author mrivero
/// \date 05/09/2012
/////////////////////////////////////////////////////////////////////////////
class ModesCalculator
{
	public:
		ModesCalculator(){};
		virtual ~ModesCalculator(){};

		virtual void calculateEigenValuesAndVectors(AnmParameters * anmParameters,
														const AnmNodeList & node_list,
														AnmEigen * eigen) = 0;

		virtual std::string generateReport() const = 0;
};

#endif /* MODESCALCULATOR_H_ */
