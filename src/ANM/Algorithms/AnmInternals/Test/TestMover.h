/////////////////////////////////////////////////////////////////////////////
/// \file TestICModesCalculator.h
///
/// \brief Class to test the ICModesCalculator class
///
/// \author vgil
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef __TESTICMOVER_H_
#define __TESTICMOVER_H_
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
class TestMover : public Test
{
    public:
					TestMover(std::string name);
        virtual 	~TestMover();
        void run();
        void init();

    private:
        void finish();
        bool testApplyRotations();
};

#endif /* TESTICMODESCALCULATOR_H_ */
