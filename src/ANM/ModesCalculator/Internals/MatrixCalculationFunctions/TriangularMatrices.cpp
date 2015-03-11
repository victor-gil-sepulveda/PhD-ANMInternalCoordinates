#include "TriangularMatrices.h"

void TriangularMatrixTools::triangularMatrixToVectorMatrix(TriangularMatrix* t, std::vector<std::vector<double> >& v){
	// both t and v must be squared
	v.clear();
	v.resize(t->size1());
	for (unsigned int i = 0; i < t->size1(); ++i){
		v[i].resize(t->size1());
	}

	for (unsigned int i = 0; i < t->size1(); ++i){
		TriangularMatrixRow trow (*t,i);
		for (unsigned int j = i; j < t->size2(); ++j){
			v[i][j] = trow(j);
			v[j][i] = trow(j);
		}
	}
}

TriangularMatrix* TriangularMatrixTools::vectorMatrixToTriangularMatrix(std::vector<std::vector<double> >& v){
	// both t and v must be squared
	TriangularMatrix* t = new TriangularMatrix(v.size(), v[0].size());

	for (unsigned int i = 0; i < v.size(); ++i){
			for (unsigned int j = i; j < v[0].size(); ++j){
				(*t)(i, j) = v[i][j];
			}
	}
	return t;
}
