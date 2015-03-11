/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 15/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef TRIANGULARMATRICES_H_
#define TRIANGULARMATRICES_H_

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace ublas = boost::numeric::ublas;

// Column major triangular matrix allows us to directly interface with lapack eigen functions.
typedef ublas::triangular_matrix<double, ublas::upper, ublas::column_major> TriangularMatrix;

// Shortcut to access to rows.
typedef ublas::matrix_row<TriangularMatrix> TriangularMatrixRow;

namespace TriangularMatrixTools{

	void triangularMatrixToVectorMatrix(TriangularMatrix* t, std::vector<std::vector<double> >& v);

	TriangularMatrix* vectorMatrixToTriangularMatrix(std::vector<std::vector<double> >& v);

};

#endif /* TRIANGULARMATRICES_H_ */
