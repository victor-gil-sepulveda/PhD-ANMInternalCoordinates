/////////////////////////////////////////////////////////////////////////////
/// \file TestANMFunctions.h
///
/// \brief Class to test the ANMFunctions class
///
/// \author arincon
/// \date 18/11/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTANMICMATH_H_
#define TESTANMICMATH_H_
#include "../../../../../Tools/TestTools.h"

class TestANMICMath : public Test
{
    public:
        TestANMICMath(std::string name);
        virtual ~TestANMICMath();
        void run();
        void init();

    private:
        void finish();

        // matrix operations
        bool testMultiplyColumnByRow();
        bool testMultiplyMatrixByScalar();
        bool testMultiplyMatrixByScalarGeneral();
        bool testAddMatrices();
		bool testAddMatricesGeneral();
		bool testAddMatricesInPlace();
		bool testAddMatricesInPlaceGeneral();
		bool testSubtractMatrices();
        bool testSubtractMatricesGeneral();
        bool testMultitplyRowByMatrix();
        bool testMultiplyRowByColumn();
        bool testMultiplyMatrixByMatrix();
        bool testMultiplyMatrixByMatrixGeneral();
        bool testMultiplyIMatrixByEVector();
        bool testInvertMatrix();
        bool testInvertMatrixGeneral();
        bool testTransposeMatrix();
        bool testTransposeMatrixGeneral();
        bool testGetABMatrix();
        bool testAddMatrixToU();
        bool testConvertMatrixToArray();
};

#endif /* TESTANMICMATH_H_ */
