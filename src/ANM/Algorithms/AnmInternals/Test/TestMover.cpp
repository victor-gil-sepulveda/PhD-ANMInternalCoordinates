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
#include "../../../../Inout/Structures/PDBWriter.h"
#include "../../../ModesCalculator/AnmEigen.h"
#include "../../../ModesCalculator/Internals/InternalModesCalculator.h"
#include "../../../ModesCalculator/PreCalculated/SimpleModesLoader.h"
#include "../../../Tools/AnmNormalizer.h"
#include "../../../AnmNodeList.h"
#include "../../../AnmUnitNodeList.h"
#include "../../../../System/Logs/ModesAndCoords/ModeWritersHandler.h"
#include "../../../../Tools/Math/Point.h"
#include "../../../ModesCalculator/PreCalculated/ModeTypes.h"
#include "../../../ModesCalculator/Internals/CoarseGrainModel/UnitTools.h"

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

    //TEST_FUNCTION(testApplyRotations)

    //TEST_FUNCTION(testIterativeAngleApplication)

    TEST_FUNCTION(testIterativeAngleApplicationWithConversionUpdate)

    //TEST_FUNCTION(testIterativeAngleApplicationWithConversionUpdateToTarget)

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


	// Get the initial angles
	vector<double> expected_initial_angles_deg, expected_initial_angles_rad, calculated_initial_angles_rad;
	TestTools::load_vector( expected_initial_angles_deg, "src/ANM/Algorithms/AnmInternals/Test/data/ala5/ala5_initial_rot.txt");
	// Convert to rads
	for(unsigned int i =0; i< expected_initial_angles_deg.size(); ++i){
		expected_initial_angles_rad.push_back(Math::degToRad(expected_initial_angles_deg[i]));
	}

	// Get the ones we calculate
	AnmInternals::calculate_current_angles(calculated_initial_angles_rad, units);
	bool we_can_calculate_angles = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_initial_angles_rad,calculated_initial_angles_rad,1e-4);

	// move once
	vector<double> rotation_increments_deg, rotation_increments_rad, expected_final_angles;
	TestTools::load_vector(rotation_increments_deg, "src/ANM/Algorithms/AnmInternals/Test/data/ala5/ala5_rot_inc.txt");
	for(unsigned int i =0; i< rotation_increments_deg.size(); ++i){
		rotation_increments_rad.push_back(Math::degToRad(rotation_increments_deg[i]));
		expected_final_angles.push_back(rotation_increments_rad[i]+calculated_initial_angles_rad[i]);
	}
	ANMICMovement::apply_rotations_to_molecule_units(units, rotation_increments_rad, 1);
//	PDBWriter writer;
//	writer.write("test.pdb", complex);

	vector<double> calculated_final_angles_rad;
	AnmInternals::calculate_current_angles(calculated_final_angles_rad, units);
	for(unsigned int i =0; i< calculated_initial_angles_rad.size(); ++i){
		cout<<"initial: "<<calculated_initial_angles_rad[i]
		    <<" increment: "<<rotation_increments_rad[i]
		    <<" expected: "<<expected_final_angles[i]
		    <<" calculated: "<<calculated_final_angles_rad[i]<<endl;
	}

	bool we_can_move_them = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_final_angles,calculated_final_angles_rad,1e-4);

	// Chain two rotations to move them to the original position
	vector<double> backward_increments;
	for(unsigned int i =0; i< rotation_increments_rad.size(); ++i){
		backward_increments.push_back(-rotation_increments_rad[i]);
	}

	ANMICMovement::apply_rotations_to_molecule_units(units, backward_increments, 0.5);
	ANMICMovement::apply_rotations_to_molecule_units(units, backward_increments, 0.5);
	AnmInternals::calculate_current_angles(calculated_final_angles_rad, units);
	for(unsigned int i =0; i< calculated_initial_angles_rad.size(); ++i){
		cout<<" calculated: "<<calculated_final_angles_rad[i]<<endl;
	}

	bool we_can_undo_the_move = Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_initial_angles_rad,calculated_final_angles_rad,1e-4);


	delete complex;
	Utils::clearVector<Unit>(units);

	return we_can_calculate_angles && we_can_move_them && we_can_undo_the_move;
}

bool apply_rots(PDBWriter& writer,
		vector<double>& angular_increments,
		vector<double>& current_angles,
		vector<Unit*>& units,
		Complex* complex,
		const char* pdb_file,
		unsigned int iteration){

	vector<double> expected_angles(angular_increments.size(), 0);

	AnmInternals::calculate_current_angles(current_angles, units);
	cout<<current_angles.size()<<" "<<angular_increments.size()<<endl;
	for(unsigned int i=0; i< current_angles.size(); ++i){
		expected_angles[i] = current_angles[i] + angular_increments[i];
	}
	cout<<"Applying rotations ("<<iteration<<")"<<endl;
	ANMICMovement::apply_rotations_to_molecule_units(units, angular_increments, 1);
	AnmInternals::calculate_current_angles(current_angles, units);
	writer.writeTrajectory(pdb_file, complex);

	return Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_angles,
			current_angles,
			1e-5);
}

bool TestMover::testIterativeAngleApplication(){
	const char* complex_pdb = "src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/9WVG.minim.pdb";
	cout << "Testing Iterative Angle application with  "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;

	bool skipOXT = true;
	TestANMICTools::createUnitsFromFile(complex_pdb,
					units,
					complex,
					not skipOXT);

	vector<double> angular_increments, current_angles;
	TestTools::load_vector( angular_increments, "src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/increments.txt");

	PDBWriter writer;
	for(unsigned int i = 0; i < 20; ++i){
		apply_rots(writer,	angular_increments,	current_angles,
				units, complex, "test1/test.pdb", i);
	}

	delete complex;
	Utils::clearVector<Unit>(units);

	return false;
}

AnmEigen* angles_from_eigen_conversion(const char* nmd_file,
		vector<Unit*>& units,
		vector<double>& fake_target_increments){

	// This nmd file has only one mode
	AnmEigen* eigen = SimpleModesLoader::load(nmd_file, true);
	AnmEigen* pca_ic_eigen =  InternalModesCalculator::cartesianToInternal(units, eigen);

	fake_target_increments.clear();
	for(unsigned int i = 0; i < pca_ic_eigen->vectors[0].size(); ++i){
		fake_target_increments.push_back(pca_ic_eigen->vectors[0][i]);
	}

	AnmNormalizer::normalizeByLargestValue(fake_target_increments);
	Math::multiplyVectorByScalar(fake_target_increments, 0.2);

	delete eigen;

	return pca_ic_eigen;
}

bool TestMover::testIterativeAngleApplicationWithConversionUpdate(){
	const char* complex_pdb = "src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/9WVG.minim.pdb";
	cout << "Testing Iterative Angle application with  "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;

	bool skipOXT = true;
	TestANMICTools::createUnitsFromFile(complex_pdb,
					units,
					complex,
					not skipOXT);

	vector<double> angular_increments, current_angles;

	PDBWriter writer;
	ModesWriterHandler whandler;
	whandler.setDirectory("/home/user/workspace/ANMIC_AZ/test2");

	for(unsigned int i = 0; i < 20; ++i){
		AnmEigen* pca_ic_eigen = angles_from_eigen_conversion("src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/cc_pca_aa.nmd",
				units,
				angular_increments);

		whandler.logStepAndVector("angular_increments", angular_increments);

		AnmNodeList* node_list = new AnmUnitNodeList;
		dynamic_cast<AnmUnitNodeList*>(node_list)->setNodeList(units);
		whandler.getWriter("iteration")->writeInternalModes(pca_ic_eigen,
				node_list,
				true);

		apply_rots(writer,	angular_increments,	current_angles,
				units, complex, "/home/user/workspace/ANMIC_AZ/test2/test.pdb", i);

		for(unsigned int j = 0; j < units.size(); ++j){
			units[j]->update();
		}

		whandler.setCounter(i);

		delete pca_ic_eigen;
	}

	delete complex;
	Utils::clearVector<Unit>(units);

	return false;
}

AnmEigen* angles_from_eigen_conversion(AnmEigen* eigen,
		vector<Unit*>& units,
		vector<double>& fake_target_increments){

	// This nmd file has only one mode
	AnmEigen* pca_ic_eigen =  InternalModesCalculator::cartesianToInternal(units, eigen);

	cout<<"Converted"<<endl;
	fake_target_increments.clear();
	for(unsigned int i = 0; i < pca_ic_eigen->vectors[0].size(); ++i){
		fake_target_increments.push_back(pca_ic_eigen->vectors[0][i]);
	}

	AnmNormalizer::normalizeByLargestValue(fake_target_increments);
	Math::multiplyVectorByScalar(fake_target_increments, 0.2);

	return pca_ic_eigen;
}

bool TestMover::testIterativeAngleApplicationWithConversionUpdateToTarget(){
	const char* complex_pdb = "src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/9WVG.minim.pdb";
	cout << "Testing Biased Iterative Angle application with  "<<complex_pdb<<endl;

	System sys;
	vector<Unit*> units;
	Complex* complex;

	bool skipOXT = true;
	TestANMICTools::createUnitsFromFile(complex_pdb,
					units,
					complex,
					not skipOXT);

	vector<double> angular_increments, current_angles;

	// Get pca mode (2)
	AnmEigen* eigen = SimpleModesLoader::load("src/ANM/Algorithms/AnmInternals/Test/data/it_ang_app/cc_pca_aa.nmd", true);

	// Calculate target points
	vector<Point> targetCCPoints;
	for (unsigned int i =0; i< eigen->vectors[0].size()/3; ++i){
		targetCCPoints.push_back(Point(&(Utils::vectorToPointer(eigen->vectors[0])[i*3])));
	}

	PDBWriter writer;
	ModesWriterHandler whandler;
	whandler.setDirectory("/home/user/workspace/ANMIC_AZ/test3");

	for(unsigned int i = 0; i < 20; ++i){

		// Calculate translations -> to AnmEigen structure
		vector<double> mode;
		vector<Atom*> all_unit_atoms;
		UnitTools::getAllAtomsFromUnits(units,all_unit_atoms, true);
		for (unsigned int j =0; j< all_unit_atoms.size();++j){
			Point translation = Point::subtract(targetCCPoints[j], all_unit_atoms[j]->toPoint());
			mode.push_back(translation.getX());
			mode.push_back(translation.getY());
			mode.push_back(translation.getZ());
		}


		vector<vector<double> > eigvecs; eigvecs.push_back(mode);
		vector<double> eigvals(1,0);
		AnmEigen * new_trans_eigen = new AnmEigen;
		new_trans_eigen->type = ModeTypes::PCA_CC;
		new_trans_eigen->initialize(eigvals, eigvecs, true);

		AnmEigen* pca_ic_eigen = angles_from_eigen_conversion(new_trans_eigen,
				units,
				angular_increments);


		vector<double> expected_angles(angular_increments.size(), 0);

		apply_rots(writer, angular_increments, current_angles,
				units, complex, "/home/user/workspace/ANMIC_AZ/test3/test.pdb", i);

		for(unsigned int j = 0; j < units.size(); ++j)
			units[j]->update();

		cout<<"units updated"<<endl;
		whandler.setCounter(i);


		AnmNodeList* node_list = new AnmUnitNodeList;
		dynamic_cast<AnmUnitNodeList*>(node_list)->setNodeList(units);

		whandler.getWriter("iteration")->writeInternalModes(pca_ic_eigen,
				node_list,
				true);

		delete pca_ic_eigen;
		delete new_trans_eigen;
	}

	delete complex;
	Utils::clearVector<Unit>(units);

	return false;
}
