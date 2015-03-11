#include "TestANMICTools.h"
#include "../../../../../Tools/TestTools.h"
#include "../../../../../Molecules/Atom.h"
#include "../../../../../Molecules/AtomElements.h"
#include "../../../../../Tools/Utils.h"
#include "../../../../../Molecules/ComplexBuilder.h"
#include "../../../../../Molecules/Selection/Selector.h"
#include "../../../../../Tools/Math/Point.h"
#include "../../../../../Molecules/Selection/SelectionBuilder.h"
#include "../../../../../Molecules/AtomSetsTree/Chains/Chain.h"
#include "../../CoarseGrainModel/UnitsBuilder.h"

#include <vector>
#include <cmath>
#include "../../../../../Molecules/Complex.h"
#include "../../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
/// bla,bla.cpp
///
/// Implementation of bla,bla class
///
/// \author myName
/// \date 13/08/2013
/////////////////////////////////////////////////////////////////////////////

Chain* TestANMICTools::createChain(const char* file_path, Complex*& r_complex){
	ComplexBuilder complexBuilder(OPLS2005);
	Complex* complex = complexBuilder.build(file_path);
	vector<Chain*> chains;
	SelectionBuilder selection_builder;
	Selector* selector = selection_builder.create("{\"chains\": \"all\"}");
	selector->getChains(chains, (AtomSet*) complex);
	r_complex = complex;
	delete selector;
	return chains[0];
}

Chain* TestANMICTools::createUnitsFromFile(const char* file_path,
												vector<Unit*>& units,
												Complex*& r_complex,
												bool skip_OXT){
	Chain * chain = createChain(file_path, r_complex);
	changeAtomMasses(chain->getAllAtoms());
	UnitsBuilder builder(chain, true, skip_OXT);
	builder.build(units);
	return chain;
}

Complex* TestANMICTools::createUnitsFromFilePickingCOM(const char* file_path,
		vector<Unit*>& units, CenterOfMass*& com, bool injectINMAMasses){

	ComplexBuilder complexBuilder(OPLS2005);
	Complex* complex = complexBuilder.build(file_path);

	if (injectINMAMasses){
		TestANMICTools::changeAtomMasses(complex->getAllAtoms());
	}

	vector<Chain*> chains;
	SelectionBuilder selection_builder;
	Selector* selector = selection_builder.create("{\"chains\": \"all\"}");
	selector->getChains(chains, complex);
	UnitsBuilder builder(chains[0]);
	CenterOfMass* tmp_cmp = builder.getCOM();
	com = new CenterOfMass(*tmp_cmp);

	builder.build(units);
	delete selector;
	return complex;
}

void TestANMICTools::createUnitsFromFileChangingCoordinates(const char* file_path,
															const char* coordinates_file_name,
															std::vector<Unit*>& units,
															Complex*& r_complex){
	Chain * chain = createChain(file_path, r_complex);
	changeAtomMasses(chain->getAllAtoms());
	UnitsBuilder builder(chain);
	changeCoordsByThoseInTheFile(chain, coordinates_file_name);
	builder.build(units);
}

bool TestANMICTools::pointsLookSimilar(double*a, double*b){
	if( abs(a[0]-b[0]) < 1e-6 &&
		abs(a[1]-b[1]) < 1e-6 &&
		abs(a[2]-b[2]) < 1e-6 ){
		return true;
	}
	else{
		return false;
	}
}

void TestANMICTools::changeCoordsByThoseInTheFile(Chain* chain, const char* coordinates_file_name){
	vector<Atom*>& atoms = chain->getAllAtoms();
	unsigned int atoms_size = atoms.size();
	vector<vector<double> > positions;
	TestTools::load_vector_of_vectors(positions, coordinates_file_name);
	vector<bool> used(atoms_size,false);
	int changed = 0;
	double* coordinates = atoms[0]->getCoordinates();
	cout<<"Changing coordinates of atoms to look like the ones in the file: "<<coordinates_file_name<<endl;
	cout<<"The file contains "<<positions.size()<<" coordinate triplets."<<endl;
	for(unsigned int i = 0; i < positions.size();++i){
		double* i_coords = Utils::vectorToPointer(positions[i]);
		for(unsigned int j = 0; j < atoms_size; ++j){
			vector<double> coords = atoms[j]->toPoint().getCoordinates();
			double* j_coords =  Utils::vectorToPointer<double>(coords);
			if(pointsLookSimilar(i_coords,j_coords) &&	!used[j]){
				coordinates[atoms[j]->ix] =i_coords[0];
				coordinates[atoms[j]->iy] =i_coords[1];
				coordinates[atoms[j]->iz] =i_coords[2];
				used[j] = true;
				changed++;
				break;
			}
		}
	}
	cout<<"Done. Changed "<<changed<<" atoms of "<<atoms_size<<endl;
}

void TestANMICTools::changeAtomMasses(std::vector<Atom*>& atoms){
	for(unsigned int i = 0; i< atoms.size(); ++i){
		if(atoms[i]->element == AtomElements::C){ // 12.010999679565
			atoms[i]->setMass(12.011);
		}
		else{
			if(atoms[i]->element == AtomElements::N){ // 14.006739616394
				atoms[i]->setMass(14.00674);
			}
			else{
				if(atoms[i]->element == AtomElements::O){ // 15.999400138855
					atoms[i]->setMass(15.9994);
				}
				else{
					if (atoms[i]->element == AtomElements::S){ // 32.066001892090
						atoms[i]->setMass(32.066);
					}
					else{
						if(atoms[i]->element == AtomElements::H){ // 1.00797
							atoms[i]->setMass(1.00797);
						}
					}
				}
			}
		}
	}
}


