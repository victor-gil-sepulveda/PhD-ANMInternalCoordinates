/////////////////////////////////////////////////////////////////////////////
/// \file TestICModesCalculator.h
///
/// \brief Class to test the ICModesCalculator class
///
/// \author vgil
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTICMODESCALCULATOR_H_
#define TESTICMODESCALCULATOR_H_
#include "../../../../Tools/TestTools.h"

#include <functional>
#include <algorithm>
#include <iterator>


/////////////////////////////////////////////////////////////////////////////
/// \brief Class to test the ICModesCalculator class
///
/// \author vgil
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////
class TestICModesCalculator : public Test
{
    public:
        TestICModesCalculator(std::string name);
        virtual ~TestICModesCalculator();
        void run();
        void init();

    private:
        void finish();

        bool testCompleteEigencalculation(const char* pdb_path,
        		const char* eigenvalues_path, const char* eigenvectors_path,
        		double, double);

        bool testEigencalculation(const char* h_path, const char* k_path,
        		const char* eigenvalues_path, const char* eigenvectors_path,
        		double, double);

        bool testInternalToCartesian(const char* prot_path,
        		const char* initial_ic_path,
        		const char* final_cc_path,
        		double tolerance);

        bool testCartesianToInternal(const char* prot_path,
        		double tolerance);

        bool testGeometricCartesianToInternal(const char* prot_path,
                		const char* initial_cc_path,
                		const char* final_ic_path,
                		double tolerance);
};

#endif /* TESTICMODESCALCULATOR_H_ */
