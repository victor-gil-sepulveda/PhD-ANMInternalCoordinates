/////////////////////////////////////////////////////////////////////////////
/// TestHarmonicDihedralConstraintTerm.cpp
///
/// Implementation of TestHarmonicDihedralConstraintTerm class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author arincon
/// \author vgil
/// \date 20/09/2012
/////////////////////////////////////////////////////////////////////////////

#include "TestHarmonicDihedralConstraintTerm.h"

#include "../HarmonicDihedralConstraintTerm.h"
#include "../../../../../Tools/Assertion.h"
#include "../../../../../System/System.h"
#include "../../../../../Molecules/ComplexBuilder.h"
#include "../../../../../Molecules/Complex.h"
#include "../../../../../Minimizers/Parameters/TNParams.h"
#include "../../../../../Minimizers/MinimizerBuilder.h"
#include "../../../../../Minimizers/Minimizer.h"
#include "../../../../EnergyCalculatorBuilder.h"
#include "../../../../EnergyCalculator.h"
#include "../../../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../../../Molecules/Selection/Selector.h"
#include "../../../ConstrainedTerms.h"
#include "../../../../../Tools/Math/Point.h"
#include "../../../../Tests/ToolsForEnergyTests.h"
#include "../../../../../Tools/Math/MathTools.h"
#include "../../../../../Minimizers/TruncatedNewton/TruncatedNewton.h"
#include "../../../../../Inout/Structures/PDBWriter.h"
#include "../../../../../Molecules/Templates/ImpactTemplateReader.h"
#include "../../../../../Energy/ForceField/ForceFieldBuilder.h"
#include "../../../../../Energy/ForceField/ForceField.h"
#include "../../../../../Molecules/Solvent/SolventGenerator.h"
#include "../../../../../Molecules/Solvent/Solvent.h"
#include "../../../../../Molecules/AtomNames.h"

using namespace std;




/**
 * This test uses the structure created ex profeso called test_dihedral.pdb
 * which defines a residue like this:
 *
 *  OH
 *  |
 *  |
 *  |
 *  |
 *  C-------------CZ
 *                |
 *                |
 *                |
 *                |
 *                |
 *                O
 *
 * The names of the atom have been chosen to be inside the 'AtomNames' namespace.
 */

TestHarmonicDihedralConstraintTerm::TestHarmonicDihedralConstraintTerm(string name)
{
    test_name = name;
    this->atoms = NULL;
}

TestHarmonicDihedralConstraintTerm::~TestHarmonicDihedralConstraintTerm(){}

void TestHarmonicDihedralConstraintTerm::init() {}

void TestHarmonicDihedralConstraintTerm::run()
{
    Test::run();

	TEST_FUNCTION(testCreation)

	TEST_FUNCTION(testGenerateFullHessianMatrixFromTriangularMatrix)

	// Move close to 90 degrees with an extra constraint and a mild force (it doesn't move too much
	// from the initial angle which is 0 / 180)
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			90,
			2,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_180deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez_with_extra_phi",
			168.,169.)

	// Move close to 90 degrees with an extra constraint and a very strong force (it moves)
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			90,
			20,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_180deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez_with_extra_phi",
			90.,92.)


	// Move from 180º to the third quadrant (overflow)
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			210,
			20,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_180deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez",
			-151,-149) // angles always in -pi,pi range

    // Move from 180º to the second quadrant (overflow)
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			-210,
			20,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_180deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez",
			149.,151.)// angles always in -pi,pi range

	// Move from 180º to the 4th quadrant
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			-45,
			20,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_180deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez",
			-46.,-44)// angles always in -pi,pi range

	// Move from 0º to the 4th quadrant
	TEST_REGRESSION_FUNCTION(testMinimizationWithDihedralConstraint,
			-45,
			20,
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral_0deg.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez",
			-46.,-44)// angles always in -pi,pi range

	finish();
}

void TestHarmonicDihedralConstraintTerm::finish()
{
    Test::finish();
}

bool TestHarmonicDihedralConstraintTerm::testGenerateFullHessianMatrixFromTriangularMatrix() {
	double triangular_hessian_array[] = {24, 35, 68, 67, 16, 25, 1, 70, 27, 42, 6, 34, 18, 15, 19,
			11, 12, 52, 37, 28, 59, 63, 75, 32, 21, 14, 43, 48, 50, 10, 60, 44, 57, 62, 54,
			4, 55, 9, 3, 17, 64, 0, 76, 31, 69, 58, 74, 23, 8, 53, 22, 51, 56, 29, 66, 65,
			49, 72, 46, 39, 7, 30, 13, 61, 45, 33, 26, 36, 71, 40, 5, 41, 73, 2, 77, 38, 20, 47};
	vector<double> triangular_hessian (triangular_hessian_array, triangular_hessian_array + sizeof(triangular_hessian_array) / sizeof(double));

	double full_matrix[][12] =
	     {{24, 35, 68, 67, 16, 25,  1, 70, 27, 42,  6, 34},
		 {0, 18, 15, 19, 11, 12, 52, 37, 28, 59, 63, 75},
		 {0,  0, 32, 21, 14, 43, 48, 50, 10, 60, 44, 57},
		 {0,  0,  0, 62, 54,  4, 55,  9,  3, 17, 64,  0},
		 {0,  0,  0,  0, 76, 31, 69, 58, 74, 23,  8, 53},
		 {0,  0,  0,  0,  0, 22, 51, 56, 29, 66, 65, 49},
		 {0,  0,  0,  0,  0,  0, 72, 46, 39,  7, 30, 13},
		 {0,  0,  0,  0,  0,  0,  0, 61, 45, 33, 26, 36},
		 {0,  0,  0,  0,  0,  0,  0,  0, 71, 40,  5, 41},
		 {0,  0,  0,  0,  0,  0,  0,  0,  0, 73,  2, 77},
		 {0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 38, 20},
		 {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 47}};

	double returned_full_matrix[12][12];
	HarmonicDihedralConstraintTerm::generateFullHessianMatrixFromTriangularMatrix(triangular_hessian, returned_full_matrix);

	bool ok = true;
	for(int i=0; i<12; ++i) {
		ok = Assertion::expectedVectorEqualsCalculatedWithinPrecision(returned_full_matrix[i], full_matrix[i], 12, 1e-13)
				and ok;
	}

	return ok;
}

Complex* load_complex(const char* pdb, const char* mol_template_path){
	TestTools::addOneImpactProteinTemplate(mol_template_path, OPLS2005, "DIEZ");
	ComplexBuilder complexBuilder(OPLS2005);
	Complex * complex = complexBuilder.build(pdb);
	TestTools::deleteOneImpactProteinTemplate(OPLS2005, "DIEZ");
	return complex;
}

void generate_constraint_nodes(vector<Atom*>& nodes, vector<string>& nodeNames, AtomSet* complex) {
	nodes.clear();
	nodeNames.clear();
	nodeNames.push_back(AtomNames::OH);
	nodeNames.push_back(AtomNames::C);
	nodeNames.push_back(AtomNames::CZ);
	nodeNames.push_back(AtomNames::O);

	for(unsigned int i=0; i<nodeNames.size(); i++) {
		SelectionBuilder selectionBuilder;
		vector<Atom*> selectedAtoms;
		//cout <<"{\"atoms\" : { \"names\":[\""+nodeNames[i]+"\"]} }"<< endl;
		Selector* selector = selectionBuilder.create("{\"atoms\" : { \"names\":[\""+nodeNames[i]+"\"]} }");
		selector->getAtoms(selectedAtoms, complex);
		delete selector;

		if(selectedAtoms.size() != 1) {
			cout << "Selected atoms length different than 1" << endl;
		}
		nodes.push_back(selectedAtoms[0]);
	}
}

bool TestHarmonicDihedralConstraintTerm::testCreation(){
	// Load complex and create nodes
	System sys;
	Complex* complex = load_complex(
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/test_dihedral.pdb",
			"src/Energy/PotentialConstraint/ConstraintTerms/HarmonicDihedralConstraint/Tests/data/diez");
	vector<Atom*> nodes;
	vector<string> nodeNames;
	generate_constraint_nodes(nodes, nodeNames, complex);

	// Generate the dihedral constraint
	HarmonicDihedralConstraintTerm constraintTerm(
			2,
			Math::degToRad(90),
			nodes);


	// Check
	// Get the names actually stored into the constraint
	vector<string> stored_names;
	for(unsigned int i = 0; i < constraintTerm.atoms.size(); i++){
		stored_names.push_back(constraintTerm.atoms[i]->name);
	}
	delete complex;

	// Do the assertions
	return Assertion::expectedEqualsCalculated(constraintTerm.springConstant, 2.) and
	Assertion::expectedEqualsCalculated(constraintTerm.equilibriumAngle, Math::degToRad(90)) and
	Assertion::expectedVectorEqualsCalculated(nodeNames, stored_names);
}

bool TestHarmonicDihedralConstraintTerm::testMinimizationWithDihedralConstraint(
		double angle,
		double force_const,
		const char* pdb,
		const char* mol_template_path,
		double ang_range_start,
		double ang_range_end) {

	// Load complex and create nodes
	System sys;
	Complex* complex = load_complex(pdb, mol_template_path);
	vector<Atom*> nodes;
	vector<string> nodeNames;
	generate_constraint_nodes(nodes, nodeNames, complex);

	// Create energy calculator for use in minimization
	SolventGenerator solventGenerator;
	Solvent * vacuum = solventGenerator.createVacuum();
	EnergyCalculatorBuilder energyCalculatorBuilder;
	EnergyCalculator* enerCalc = ToolsForEnergyTests::createSimpleEnergyCalculator( complex, vacuum, 12.);

	// With this torsional constraint
	HarmonicDihedralConstraintTerm* dihedral_constraint = new HarmonicDihedralConstraintTerm(
			force_const,
			Math::degToRad(angle),
			nodes);
	enerCalc->addAnmConstraintTerm(dihedral_constraint);

	// Write the initial angle
	double initial_angle = dihedral_constraint->getCurrentAngle(&complex->getAllCartesianCoordinates()[0]);
	cout<<"Initial angle: "<<initial_angle<<endl;

	// Create minimizer and use it
	MinParams* min_params = new TNParams;
	Minimizer minimizer(new TruncatedNewton);
	minimizer.minimize(min_params,
			complex,
			enerCalc);

	// Calculate the final angle
	double current_angle = dihedral_constraint->getCurrentAngle(&complex->getAllCartesianCoordinates()[0]);
	cout<<"Final angle: "<<current_angle<<endl;

//	PDBWriter writer;
//	std::ostringstream out;
//	writer.writeAllChainsToPdb(complex, out);
//	cout<<out.str()<<endl;

	delete enerCalc;
	delete min_params;
	delete complex;
	delete vacuum;

	return Assertion::isInRange(current_angle,ang_range_start,ang_range_end);
}
