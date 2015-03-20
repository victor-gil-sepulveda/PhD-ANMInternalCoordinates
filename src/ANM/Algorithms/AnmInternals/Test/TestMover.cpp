/////////////////////////////////////////////////////////////////////////////
/// TestMover.cpp
///
/// Implementation of TestMover class
///
/// \author vgil
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "TestMover.h"
#include "../../../../Tools/Utils.h"
#include "../../../ModesCalculator/Internals/MatrixCalculationFunctions/Tests/TestANMICTools.h"
#include "../../../../Molecules/ForceFieldTerms/Dihedral.h"
#include "../../../../Molecules/Atom.h"
#include "../../../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../../../../System/System.h"
#include "../../../../Molecules/Complex.h"
#include "../AnmInternals.h"
#include "../../../../Tools/Math/MathTools.h"
#include "../Movement/ANMICMovement.h"
#include "../../../../Tools/Assertion.h"

using namespace std;

TestMover::TestMover(string name)
{
    test_name = name;
}

TestMover::~TestMover()
{}

void TestMover::init()
{}

void TestMover::run()
{
    Test::run();

    TEST_FUNCTION(testApplyRotations)

    finish();
}

void TestMover::finish()
{
    Test::finish();
}

bool TestMover::testApplyRotations(){

	const char* complex_pdb = "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala5/ala5.fixed.pdb";
	cout << "Testing rotations with "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;

	bool skipOXT = true;
	TestANMICTools::createUnitsFromFile(complex_pdb,
					units,
					complex,
					not skipOXT);

//	for (unsigned int i = 0; i< units.size(); ++i){
//		cout<<"get_dihedral ";
//		Dihedral* d = units[i]->right_dihedral;
//		vector<Atom*> d_atoms = d->getAtoms();
//		for (unsigned int j = 0; j< d_atoms.size(); ++j){
//			cout<<"resi "<<d_atoms[j]->resSeq<<" and name "<<d_atoms[j]->name<<", ";
//		}
//		cout<<endl;
//	}

	// Get the initial angles
	vector<double> expected_initial_angles_deg, expected_initial_angles_rad, calculated_initial_angles_rad;
	TestTools::load_vector( expected_initial_angles_deg, "src/ANM/Algorithms/AnmInternals/Test/data/ala5/ala5_initial_rot.txt");
	// Convert to rads
	for(unsigned int i =0; i< expected_initial_angles_deg.size(); ++i){
		expected_initial_angles_rad.push_back(Math::degToRad(expected_initial_angles_deg[i]));
	}
	// Get the ones we calculate
	AnmInternals::calculate_current_angles(calculated_initial_angles_rad, units);
//	for(unsigned int i =0; i< expected_initial_angles_deg.size(); ++i){
//		cout<<calculated_initial_angles_rad[i]<<" "<<expected_initial_angles_rad[i]<<endl;
//	}
	bool we_can_calculate_angles = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_initial_angles_rad,calculated_initial_angles_rad,1e-4);

	// move once
	vector<double> rotation_increments_deg, rotation_increments_rad, expected_final_angles;
	TestTools::load_vector(rotation_increments_deg, "src/ANM/Algorithms/AnmInternals/Test/data/ala5/ala5_rot_inc.txt");
	for(unsigned int i =0; i< rotation_increments_deg.size(); ++i){
		rotation_increments_rad.push_back(Math::degToRad(rotation_increments_deg[i]));
		expected_final_angles.push_back(rotation_increments_rad[i]+calculated_initial_angles_rad[i]);
	}
	ANMICMovement::apply_rotations_to_molecule_units(units, rotation_increments_rad, 1);
	vector<double> calculated_final_angles_rad;
	AnmInternals::calculate_current_angles(calculated_final_angles_rad, units);
	for(unsigned int i =0; i< calculated_initial_angles_rad.size(); ++i){
		cout<<"initial: "<<calculated_initial_angles_rad[i]
		    <<" increment: "<<rotation_increments_rad[i]
		    <<" expected: "<<expected_final_angles[i]
		    <<" calculated: "<<calculated_final_angles_rad[i]<<endl;
	}

	// move twice
	//	ANMICMovement()

	bool we_can_move_them = false;

	delete complex;
	Utils::clearVector<Unit>(units);

	return we_can_calculate_angles && we_can_move_them;
}
