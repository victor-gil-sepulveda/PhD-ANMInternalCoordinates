/////////////////////////////////////////////////////////////////////////////
/// \author vgil
/// \date 17/03/2015
/////////////////////////////////////////////////////////////////////////////

#include "TriangularMatrices.h"

//SQUARE vector of vectors to triangular SQUARE matrix
TriangularMatrix* TriangularMatrixTools::vectorMatrixToTriangularMatrix(std::vector<std::vector<double> >& vm){
	TriangularMatrix* tm = new TriangularMatrix(vm.size(), vm.size());

	for(unsigned int i = 0; i < vm.size();++i){
		for(unsigned int j = i; j < vm[i].size(); ++j){
			(*tm)(i,j) = vm[i][j];
		}
	}

	return tm;
}

//SQUARE triangular SQUARE matrix to vector of vectors
void TriangularMatrixTools::triangularMatrixToVectorMatrix(TriangularMatrix* tm,
		std::vector<std::vector<double> >& vm){
	vm.clear();
	vm.resize(tm->size1());
	for(unsigned int i = 0; i < vm.size();++i){
		vm[i].resize(tm->size1());
	}

	for(unsigned int i = 0; i < vm.size();++i){
		for(unsigned int j = i; j < vm[i].size(); ++j){
			vm[i][j] = (*tm)(i,j);
			vm[j][i] = (*tm)(i,j);
		}
	}
}
