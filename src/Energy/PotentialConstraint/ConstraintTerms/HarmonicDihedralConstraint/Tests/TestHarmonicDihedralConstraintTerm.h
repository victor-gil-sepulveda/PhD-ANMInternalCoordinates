/////////////////////////////////////////////////////////////////////////////
/// \file TestHarmonicDihedralConstraintTerm.h
///
/// \brief Class to test the HarmonicDihedralConstraintTerm class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author arincon
/// \date 20/09/2012
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTHARMONICDIHEDRALCONSTRAINTTERM_H_
#define TESTHARMONICDIHEDRALCONSTRAINTTERM_H_
#include "../../../../../Tools/TestTools.h"

/////////////////////////////////////////////////////////////////////////////
/// \brief Class to test the HarmonicDihedralConstraintTerm class
///
/// \author arincon
/// \date 20/09/2012
/////////////////////////////////////////////////////////////////////////////
class TestHarmonicDihedralConstraintTerm : public Test
{
    public:
        TestHarmonicDihedralConstraintTerm(std::string name);
        virtual ~TestHarmonicDihedralConstraintTerm();
        void run();
        void init();

    private:
        void finish();

        bool testGenerateFullHessianMatrixFromTriangularMatrix();

        bool testCreation();

        bool testMinimizationWithDihedralConstraint(double angle,
        		double force_const,
        		const char* pdb,
        		const char* mol_template_path,
        		double ang_range_start,
        		double ang_range_end);


        Atom** atoms;
};

#endif /* TESTHARMONICDIHEDRALCONSTRAINTTERM_H_ */
