/////////////////////////////////////////////////////////////////////////////
/// TestICModesCalculator.cpp
///
/// Implementation of TestICModesCalculator class
///
/// \author vgil
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "TestICModesCalculator.h"
#include "../../../../Tools/TestTools.h"
#include "../MatrixCalculationFunctions/TriangularMatrices.h"
#include "../../../../Tools/Assertion.h"
#include "../InternalModesCalculator.h"
#include "../../../ModesCalculator/AnmEigen.h"
#include "../../../Parameters/AnmParameters.h"
#include "../CoarseGrainModel/Unit.h"
#include "../../../../Molecules/Atom.h"
#include "../../../../Tools/Math/Point.h"
#include "../MatrixCalculationFunctions/ANMICMath.h"
#include "../MatrixCalculationFunctions/Tests/TestANMICTools.h"
#include "../../../../System/System.h"
#include "../../../../Molecules/Complex.h"
#include "../../../../Tools/Math/MathTools.h"

#include <algorithm>
#include "../../../../Tools/Utils.h"
#include "../CoarseGrainModel/UnitTools.h"

#include<iomanip>
#include "../../../AnmUnitNodeList.h"
#include "../../../AnmNodeList.h"
#include "../../../Tools/Inout/ModesWriter.h"
#include "../../PreCalculated/SimpleModesLoader.h"

using namespace std;

TestICModesCalculator::TestICModesCalculator(string name)
{
    test_name = name;
}

TestICModesCalculator::~TestICModesCalculator()
{}

void TestICModesCalculator::init()
{}

void TestICModesCalculator::run()
{
    Test::run();

    TEST_FUNCTION(testCompleteEigencalculation,
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/1AKE.pdb",
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/eigenvalues.txt",
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/6_eigenvectors.txt",
			1e-4,1e-7);

    TEST_FUNCTION(testCompleteEigencalculation,
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/eigenvalues.txt",
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/6_eigenvectors.txt",
			1e-4,1e-7);

    TEST_FUNCTION(testEigencalculation,
    		"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc1/h.txt",
    		"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc1/k.txt",
    		"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc1/eigenvalues.txt",
    		"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc1/eigenvectors.txt",
    		1e-4,1e-6);

    // Our calculated evectors can reproduce iNMA eigenstuff
    TEST_FUNCTION(testEigencalculation,
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/H.txt",
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/1AKE/K.txt",
			"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc2/eigenvalues.txt",
			"src/ANM/ModesCalculator/Internals/Tests/data/eigencalc2/eigenvectors.txt",
			1e-4,1e-6);

    TEST_FUNCTION(testInternalToCartesian,
    		"src/ANM/ModesCalculator/Internals/Tests/data/conversion/4AKE/4AKE.pdb",
    		"src/ANM/ModesCalculator/Internals/Tests/data/conversion/4AKE/4akeicevec.txt",
    		"src/ANM/ModesCalculator/Internals/Tests/data/conversion/4AKE/4akeccevec.txt",
    		0.005);

    TEST_FUNCTION(testInternalToCartesian,
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1DDT_s1/1DDT_s1.pdb",
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1DDT_s1/ic_evec.txt",
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1DDT_s1/cc_evec.txt",
			0.05);

    TEST_FUNCTION(testInternalToCartesian,
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1SU4/1su4.pdb",
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1SU4/ic_evec.txt",
			"src/ANM/ModesCalculator/Internals/Tests/data/conversion/1SU4/cc_evec.txt",
			0.08);

    FAILING_TEST_FUNCTION(testCartesianToInternal,
    			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
    			1e-12);


    FAILING_TEST_FUNCTION(testConvert9WVGStuff)


    FAILING_TEST_FUNCTION(testGeometricCartesianToInternal,
    			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ala3.pdb",
    			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/cc_evec.txt",
    			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/ala3/ic_evec.txt",
    			1e-12);

    finish();
}

void TestICModesCalculator::finish()
{
    Test::finish();
}

bool check_eigen_results(AnmEigen& eigen,
		unsigned int number_of_modes,
		const char* eigenvalues_path,
		const char* eigenvectors_path,
		double eig_value_precision,
		double eig_vector_precision){

	// Load expected golden data
	vector<double> expected_eigenvalues;
	vector< vector<double> > expected_eigenvectors;
	TestTools::load_vector(expected_eigenvalues, eigenvalues_path);
	TestTools::load_vector_of_vectors(expected_eigenvectors, eigenvectors_path);

	// Check the results
	bool ok = true;
	cout<<"Checking eigenvalues..."<<endl;
	for(unsigned int i = 0; i< number_of_modes;++i){
		ok = ok && Assertion::expectedEqualsCalculatedWithinPrecision(expected_eigenvalues[i],
				eigen.getEigenValueOfMode(i), eig_value_precision);
	}

	cout<<"Checking eigenvectors..."<<endl;
	for(unsigned int i = 0; i< number_of_modes;++i){
		ok = ok && Assertion::expectedVectorEqualsCalculatedWithinPrecision(expected_eigenvectors[i],
				eigen.getEigenVectorOfMode(i), eig_vector_precision);
	}

	return ok;
}

bool TestICModesCalculator::testCompleteEigencalculation(const char* pdb_path,
		const char* eigenvalues_path,
		const char* eigenvectors_path,
		double eig_value_precision,
		double eig_vector_precision){

	// Prepare params : cutoff= 10.0, k= 1.000000, x0= 3.8
	AnmParameters anmParameters;
	anmParameters.setNumberOfModes(1);
	anmParameters.setCutoff(10.0);
	anmParameters.setConstantForHessian(1.0);

	// Load structure
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model"<<endl;
	TestANMICTools::createUnitsFromFile(pdb_path, units, complex, false);

	//	anmParameters.getNumberOfModes()
	AnmEigen eigen;
	AnmUnitNodeList* node_list = new AnmUnitNodeList;
	node_list->setNodeList(complex);
	node_list->setNodeList(units); // OJOOOOOOOOOOOOOO
	InternalModesCalculator calctor;
	calctor.calculateEigenValuesAndVectors(&anmParameters,
											*((AnmNodeList*)node_list),
											&eigen);
	delete node_list;

	return check_eigen_results(eigen,
			anmParameters.getNumberOfModes(),
			eigenvalues_path,
			eigenvectors_path,
			eig_value_precision, eig_vector_precision);

}

bool TestICModesCalculator::testEigencalculation(const char* h_path,
		const char* k_path,
		const char* eigenvalues_path,
		const char* eigenvectors_path,
		double eig_value_precision,
		double eig_vector_precision){

	vector<vector<double> > Hv, Kv;

	TestTools::load_vector_of_vectors(Hv, h_path);
	TriangularMatrix H(Hv.size(),Hv.size());
	for(unsigned int i = 0; i< Hv.size();++i){
		for(unsigned int j = i; j < Hv.size();++j){
			H(i,j) = Hv[i][j];
		}
	}

	TestTools::load_vector_of_vectors(Kv, k_path);
	TriangularMatrix K(Kv.size(),Kv.size());
	for(unsigned int i = 0; i< Kv.size();++i){
		for(unsigned int j = i; j < Kv.size();++j){
			K(i,j) = Kv[i][j];
		}
	}

	AnmParameters anmParameters;
	anmParameters.setNumberOfModes(5);

	AnmEigen eigen;
	InternalModesCalculator::calculate_modes(&anmParameters, &H, &K, &eigen);

	// 1e-4 1e-6
	return check_eigen_results(eigen,
			anmParameters.getNumberOfModes(),
			eigenvalues_path,
			eigenvectors_path,
			eig_value_precision, eig_vector_precision);
}

bool TestICModesCalculator::testInternalToCartesian(const char* prot_path,
		const char* initial_ic_path,
		const char* final_cc_path,
   		double tolerance){
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model"<<endl;
	TestANMICTools::createUnitsFromFile(prot_path, units, complex, false);

	//UnitTools::printUnits(units);

	AnmEigen eigen_ic;
	double eigenvalue = 1; // Dummy value for eigenvalue
	vector<double> ic_eigenvector;
	cout << "Loading vector (ic)"<<endl;
	TestTools::load_vector(ic_eigenvector, initial_ic_path);
	eigen_ic.initialize( &eigenvalue, 	// values
			&(ic_eigenvector[0]), 		// vectors
			1, 							// num. modes
			ic_eigenvector.size(), 		// num. nodes
			false);						// is cartesian?

	cout << "Loading vector (cc golden)"<<endl;
	vector<double> cc_eigenvector;
	TestTools::load_vector(cc_eigenvector, final_cc_path);


	cout << "Calculating cc from ic"<<endl;
	AnmEigen* eigen_cc = InternalModesCalculator::internalToCartesian(units, &eigen_ic, true);

	vector<double> calculated_cc = eigen_cc->getEigenVectorOfMode(0);
	cout<<"Mode size: "<<calculated_cc.size()<<endl;
//	for(unsigned int i = 0; i < calculated_cc.size(); i++){
//		cout<<calculated_cc[i]<<endl;
//	}

	Math::normalizeVector(&(calculated_cc[0]),calculated_cc.size());
	Math::normalizeVector(&(cc_eigenvector[0]), cc_eigenvector.size());

	// Atom ordering from INMA is different from ours. As our golden data file for testing comes from
	// INMA a direct comparison is not possible. A workaround would be to order both vectors and compare.
	// The test however is weaker than the real one-to-one comparison.

	std::stable_sort (cc_eigenvector.begin(), cc_eigenvector.end());
	std::stable_sort (calculated_cc.begin(), calculated_cc.end());

	bool ok = Assertion::expectedVectorEqualsCalculatedWithTolerance(
							cc_eigenvector,
							calculated_cc,
							tolerance);

	delete complex;
	Utils::clearVector<Unit>(units);
	return ok;
}

bool TestICModesCalculator::testCartesianToInternal(const char* prot_path,
       		double tolerance){

	// Prepare params : cutoff= 10.0, k= 1.000000, x0= 3.8
	AnmParameters anmParameters;
	anmParameters.setNumberOfModes(4);
	anmParameters.setCutoff(10.0);
	anmParameters.setConstantForHessian(1.0);

	// Load structure
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model."<<endl;
	TestANMICTools::createUnitsFromFile(prot_path, units, complex, false);

	//	Calculate the modes
	cout << "IC modes calculation."<<endl;
	AnmEigen eigen;
	AnmUnitNodeList* node_list = new AnmUnitNodeList;
	node_list->setNodeList(complex);
	node_list->setNodeList(units);
	InternalModesCalculator calctor;
	calctor.calculateEigenValuesAndVectors(&anmParameters,
											*((AnmNodeList*)node_list),
											&eigen);

	// Calculate the cartesian modes
	cout << "IC -> CC conversion."<<endl;
	AnmEigen* eigen_cc = InternalModesCalculator::internalToCartesian(units, &eigen, true);

	// Calculate the internal modes from cartesian modes
	cout << "CC -> IC conversion."<<endl;
	AnmEigen* output_eigen =  InternalModesCalculator::cartesianToInternal(units, eigen_cc);

	for(unsigned int i = 0; i < eigen.vectors.size();++i){
		cout<<"---------------"<<endl;
		cout<<"Calculated: ";ANMICMath::printVector(Utils::vectorToPointer(output_eigen->vectors[i]), output_eigen->vectors[i].size());
		cout<<"Expected:   ";ANMICMath::printVector(Utils::vectorToPointer(eigen.vectors[i]), eigen.vectors[i].size());
	}

	delete node_list;
	delete eigen_cc;
	delete output_eigen;
	delete complex;

	return false;
}

bool TestICModesCalculator::testGeometricCartesianToInternal(const char* prot_path,
       		const char* initial_cc_path,
       		const char* final_ic_path,
       		double tolerance){

	// Load structure
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model"<<endl;
	TestANMICTools::createUnitsFromFile(prot_path, units, complex, false);

	// Load input and golden
	vector<vector<double> > input_cartesian, output_internal, expected_internal;
	TestTools::load_vector_of_vectors(input_cartesian, initial_cc_path);
	TestTools::load_vector_of_vectors(expected_internal, final_ic_path);

	// Create a AnmEigen object to hold input data
	vector<double> eigenvalues(input_cartesian.size(), 1);
	AnmEigen input_cc;
	bool usingCC = true;
	input_cc.initialize(eigenvalues, input_cartesian, usingCC);

//	do i = 1, nomega
//	         base = iomega(1,i)
//	         partner = iomega(2,i)
//	         call rotlist (base,partner)
//	         xdist = x(base) - x(partner)
//	         ydist = y(base) - y(partner)
//	         zdist = z(base) - z(partner)
//	         norm = sqrt(xdist**2 + ydist**2 + zdist**2)
//	         xdist = xdist / norm
//	         ydist = ydist / norm
//	         zdist = zdist / norm
	// calculate ea

//	         do j = 1, nrot
//	            k = rot(j)
//	            xatom = x(k) - x(base)
//	            yatom = y(k) - y(base)
//	            zatom = z(k) - z(base)
	// calculate r-ra
//	            xterm = ydist*zatom - zdist*yatom
//	            yterm = zdist*xatom - xdist*zatom
//	            zterm = xdist*yatom - ydist*xatom
	// calculate (r-ra) x ea
//	            teb(i) = teb(i) + deb(1,k)*xterm + deb(2,k)*yterm
//	     &                              + deb(3,k)*zterm
//	         end do
//	      end do
//
//	      do i = 1, nomega
//	         tesum(i) = teb(i)
//	         derivs(i) = tesum(i)
//	      end do
//
//	nomega   number of dihedral angles allowed to rotate
//	iomega   numbers of two atoms defining rotation axis
//	rotlist     creates a list of all atoms affectd by a rotation
//	nrot        total number of atoms moving when bond rotates
//	rot         atom numbers of atoms moving when bond rotates

	vector<double>& orig_mode = input_cartesian[0];
	unsigned int number_of_torsions = units.size()-1;
	vector<double> conversion(number_of_torsions, 0);
	for (unsigned int alpha=0; alpha < number_of_torsions; ++alpha){
		vector<Atom*> right_atoms;
		UnitTools::getAllAtomsFromUnitRange(units, right_atoms, 0, alpha,true);
		Point* ea = units[alpha]->e_right;
		Point* ra = units[alpha]->r_right;
		for(unsigned int i = 0; i < right_atoms.size(); ++i){
			unsigned int offset = i*3;
			Point grad(orig_mode[offset], orig_mode[offset+1], orig_mode[offset+2]);
			Point r = right_atoms[i]->toPoint();
			Point r_ra = Point::subtract(r,*ra);
			Point eaxr_ra = ANMICMath::crossProduct(*ea, r_ra);
			conversion[alpha] += Point::dotProduct(grad, eaxr_ra);
		}
	}

	ANMICMath::printVector(Utils::vectorToPointer(conversion), conversion.size());

	return false;
}


bool TestICModesCalculator::testConvert9WVGStuff(){
		//	cout<<"K loading..."<<endl;
		//	vector< vector<double> > Kv;
		//	TestTools::load_vector_of_vectors(Kv, "src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/K.txt");
		//	// recover symmetry and printout matrix
		//	cout<<"********************************"<<endl;
		//	for(unsigned int i = 0; i< Kv.size();++i){
		//		for (unsigned int j = 0; j< Kv.size();++j){
		//			if(i<=j){
		//				cout<<Kv[i][j]<<" ";
		//			}
		//			else{
		//				cout<<Kv[j][i]<<" ";
		//			}
		//		}
		//		cout<<endl;
		//	}
		//
		//	cout<<"********************************"<<endl;


	// Load structure
	System sys;
	vector<Unit*> units;
	Complex* complex;
	cout << "Loading model"<<endl;
	TestANMICTools::createUnitsFromFile(
			"src/ANM/ModesCalculator/Internals/MatrixCalculationFunctions/Tests/data/9WVG/9WVG.pdb",
			units, complex, false);

//	// Write the atoms
//	vector<Atom*> atoms;
//	UnitTools::getAllAtomsFromUnits(units, atoms, true);
//	for(unsigned int i = 0; i< atoms.size();++i){
//		cout<<atoms[i]->name<<" "<<atoms[i]->resSeq<<" "<<atoms[i]->getX()<<" "<<atoms[i]->getY()<<" "<<atoms[i]->getZ()<<" "<<endl;
//	}

	// Get pca modes
	vector<double> values;
	vector<vector<double> > vectors;
	TestTools::load_vector(values, "/home/user/Desktop/ANM_analysis/IC/ClusterExperiments/cc_pca/BB_nozeros/cc_pca_aa.values");
	TestTools::load_vector_of_vectors(vectors, "/home/user/Desktop/ANM_analysis/IC/ClusterExperiments/cc_pca/BB_nozeros/cc_pca_aa.vectors");
	AnmEigen pca_modes;
	pca_modes.initialize(values,vectors,true);

	// Perform the conversion
	AnmEigen* pca_ic_eigen =  InternalModesCalculator::cartesianToInternal(units, &pca_modes);

	// Write them!
	ModesWriter writer;

	AnmUnitNodeList nodeList;
	nodeList.setNodeList(complex);
	nodeList.setNodeList(units);

	writer.setPath("/home/user/Desktop/ANM_analysis/IC/ClusterExperiments/cc_pca/BB_nozeros/9WVG_PCA_CCtoIC.nmd");
	writer.setName("9WVG_PCA_CCtoIC");
	writer.writeInternalModes( pca_ic_eigen, &nodeList, false);

	writer.setPath("/home/user/Desktop/ANM_analysis/IC/ClusterExperiments/cc_pca/BB_nozeros/9WVG_PCA_CCtoICtoCC.nmd");
	writer.setName("9WVG_PCA_CCtoICtoCC");
	writer.writeInternalModes( pca_ic_eigen, &nodeList, true);

	delete complex;
	Utils::clearVector(units);
	delete pca_ic_eigen;

	return false;
}
