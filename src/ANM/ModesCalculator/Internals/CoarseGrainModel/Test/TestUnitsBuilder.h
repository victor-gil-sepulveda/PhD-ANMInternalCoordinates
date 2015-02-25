/////////////////////////////////////////////////////////////////////////////
/// \file TestUnitsBuilder.h
///
/// \brief Class to test the UnitsBuilder class
///
/// \author vgil
/// \date 06/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TESTUNITSBUILDER_H_
#define TESTUNITSBUILDER_H_
#include "../../../../../Tools/TestTools.h"


/////////////////////////////////////////////////////////////////////////////
/// \brief Class to test Units and UnitsBuilder class
///
/// \author vgil
/// \date 06/08/2013
/////////////////////////////////////////////////////////////////////////////
class TestUnitsBuilder : public Test
{
    public:
        TestUnitsBuilder(std::string name);
        virtual ~TestUnitsBuilder();
        void run();
        void init();

        bool testUnitaryBondVectors(const char* complex_path, const char* vectors_file_path);

        bool testNunitsCenterAndMass(const char* complex_path, const char* expected_data_file);

    private:
        void finish();
};

#endif /* TESTUNITSBUILDER_H_ */
