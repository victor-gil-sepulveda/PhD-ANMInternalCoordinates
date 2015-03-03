/////////////////////////////////////////////////////////////////////////////
/// \file TestANMKineticMatrixFunctions.h
///
/// \brief Class to test the ANMKineticMatrixFunctions class
///
/// \author vgil
/// \author arincon
/// \date 12/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTANMKINETICMATRIXFUNCTIONS_H_
#define TESTANMKINETICMATRIXFUNCTIONS_H_
#include "../../../../../Tools/TestTools.h"
#include "../TensorCalcTypes.h"

class TestANMKineticMatrixFunctions : public Test
{
    public:
        TestANMKineticMatrixFunctions(std::string name);
        virtual ~TestANMKineticMatrixFunctions();
        void run();
        void init();

    private:
        void finish();

        bool testCalculateI(const char* ,const char* ,ICalcType ,bool ,double);
        bool testCalculateK(const char*,const char*, ICalcType, double);
        bool testJacobi(const char*,	double );

};

#endif /* TESTANMKINETICMATRIXFUNCTIONS_H_ */
