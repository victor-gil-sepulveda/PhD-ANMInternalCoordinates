/////////////////////////////////////////////////////////////////////////////
/// \file TestANMFunctions.h
///
/// \brief Class to test the ANMFunctions class
///
/// \author arincon
/// \date 18/11/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTANMFUNCTIONS_H_
#define TESTANMFUNCTIONS_H_
#include "../../../../../Tools/TestTools.h"
class Unit;
class Chain;
class Complex;


/////////////////////////////////////////////////////////////////////////////
/// \brief Class to test the ANMFunctions class
///
/// \author arincon
/// \date 18/11/2012
/////////////////////////////////////////////////////////////////////////////
class TestANMInternalHessianFunctions : public Test
{
    public:
        TestANMInternalHessianFunctions(std::string name);
        virtual ~TestANMInternalHessianFunctions();
        void run();
        void init();

    private:
        void finish();

        bool testCalculateD(const char* golden_file);
        bool testCalculateT(const char* complex_pdb, const char* goldenT, bool skip_OXT, double test_precision);
        bool testCalculateU(const char* complex_pdb, const char* goldenU, bool skip_OXT, double test_precision);
        bool testCalculateR(const char* complex_pdb, const char* goldenR, bool skip_OXT, double test_precision);
        bool testCalculateH(const char* complex_pdb, const char* goldenH, bool skip_OXT, double test_precision);
        bool testModifyHessianWithExtraTorsion(const char* initialH, const char* goldenHmodif);
};

#endif /* TESTANMFUNCTIONS_H_ */
