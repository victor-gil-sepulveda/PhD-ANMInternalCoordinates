/////////////////////////////////////////////////////////////////////////////
/// \file TestDihedralConstraintFunctions.h
///
/// \brief Class to test the Template class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \\author arincon
/// \date 25/05/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTDIHEDRALCONSTRANTFUNCTION_H_
#define TESTDIHEDRALCONSTRANTFUNCTION_H_

#include <iosfwd>
#include "../../../../../../Tools/TestTools.h"

/////////////////////////////////////////////////////////////////////////////
/// \brief Class to test the Template class
///
/// \\author arincon
/// \date 25/05/2012
/////////////////////////////////////////////////////////////////////////////
class TestDihedralConstraintFunctions : public Test
{
    public:
		TestDihedralConstraintFunctions(std::string name);
        virtual ~TestDihedralConstraintFunctions();
        void run();
        void init();

    private:
        void finish();

        // Test methods
        bool testAngleFixing(double angle, double expected);
        bool testPIto2PIRange(double angle, double expected);
        bool testRadSubtraction(double a, double b, double expected_increment);
        bool testDihedralAngleWithArcTanFunction(double precision);
        bool testEnergyWithArcTanFunction(double precision);
        bool testDerivatives(double precision);
        bool testHessian(double precision);
};

#endif /* TESTDIHEDRALCONSTRANTFUNCTION_H_ */
