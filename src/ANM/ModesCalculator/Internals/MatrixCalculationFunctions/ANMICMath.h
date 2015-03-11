/////////////////////////////////////////////////////////////////////////////
/// \file ANMICMath.h
///
/// \brief Mathematical functions used in the computation of the hessian
///        and kinetic matrix
///
/// \author arincon
/// \author vgil
/// \date 10/08/2013
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMICMATH_H_
#define ANMICMATH_H_

#include <vector>
#include "../../../../Tools/Math/Point.h"
#include "../../../../Molecules/Atom.h"
#include "../CoarseGrainModel/Unit.h"
#include "TriangularMatrices.h"


namespace ANMICMath {

	void multiplyColumnByRow(double a[6], double b[6], double result[6][6]);

	void multiplyColumnByRow(double *a, double *b, double *result, unsigned int size);

	void multiplyMatrixByScalar(double scalar, double result[3][3]);

	void multiplyMatrixByScalar(double scalar, double result[6][6]);

	void multiplyMatrixByScalar(double scalar, double *result, unsigned int size);

	void addMatrices(double a[6][6], double b[6][6], double result[6][6]);

	void addMatrices(double *a, double *b, double *result, unsigned int size);

	void addMatrices(double const matrixA[3][3], double const matrixB[3][3], double result[3][3]);

	void addMatricesInPlace(double a[3][3], double b[3][3]);

	void addMatricesInPlace(double *a, double *b, unsigned int size);

	void subtractMatrices(double const matrixA[3][3], double const matrixB[3][3], double result[3][3]);

	void subtractMatrices(double const *matrixA, double const *matrixB, double *result, unsigned int size);

	void multitplyRowByMatrix(double *row, double *matrix, double *result, unsigned int size);

	Point multitplyPointByMatrix(Point& row, const double matrix[3][3]);

	double multiplyRowByColumn(double *a, double *b, unsigned int size);

	void multiplyMatrixByMatrix(double first[3][3], double second[3][3], double result[3][3]);

	void multiplyMatrixByMatrix(double *first, double *second, double *result, unsigned int size);

	void multiplyMatrixByMatrix(std::vector<std::vector<double> >& first,
			std::vector<std::vector<double> >& second,
			std::vector<std::vector<double> >& out);

	Point multiplyIMatrixByEVector(double const I[3][3], double e[3]);

	Point multiplyIMatrixByEVector(double const I[3][3], Point& e);

	void invertIMatrix(double matrix[3][3], double result[3][3]);

	void invertMatrix(double *matrix, double *result, unsigned int size);

	void transposeMatrix(double matrix[3][3], double result[3][3]);

	void transposeMatrix(double *matrix, double *result, unsigned int size);

	void getABMatrix(std::vector< std::vector<double> > &UMatrix, double *result, int a, int b);

	void addMatrixToU(std::vector< std::vector<double> > &UMatrix, double *Tab, int a, int b);

	void subsMatrixFromU(std::vector< std::vector<double> > &UMatrix, double *Tab, int a, int b);

	Point crossProduct(Point a, Point b);

	double dotProduct(Point a, Point b);

	void printMatrix(double matrix[6][6]);

	void printMatrix(double matrix[3][3]);

	void printMatrix(std::vector<std::vector<double> >& matrix);

	void printMatrix(double *matrix, int length);

	void printVector(double *vector, int length);

	void convertMatrixToArray(double **matrix, double *array, int rows, int columns);

	Unit* getDummyUnit(double x, double y, double z);

	void transpose(std::vector<std::vector<double> >& in,
			std::vector<std::vector<double> >& out);

	void invertMatrix(std::vector<std::vector<double> >& in);

	void invertMatrix(TriangularMatrix* in_out);

};

#endif /* ANMICMATH_H_ */
