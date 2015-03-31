/////////////////////////////////////////////////////////////////////////////
/// ANMICMath.cpp
///
/// Implementation of ANMICMath functions
///
/// \author arincon
/// \author vgil
/// \date 10/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "ANMICMath.h"
#include <vector>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include "../../../../Tools/Utils.h"
using namespace std;

extern "C" {
	// Symmetric positive definite inversion
	void dpptrf_(char* uplo, int* n, double* ap, int* info);
	void dpptri_(char* uplo, int* n, double* ap, int* info);

	// Symmetric indefinite inversion
	void dsptrf_( char* uplo, int* n, double* ap, int* ipv, int* info);
	void dsptri_( char* uplo, int* n, double* ap, int* ipv, double* work, int* info);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the matrix resulting of the
/// multiplication of a column matrix by a row matrix
///
/// \param a [In] Coordinates of the point a
/// \param b [In] Coordinates of the point a
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyColumnByRow(double a[6], double b[6], double result[6][6]) {
	for(unsigned int i=0; i<6; i++) {
		for(unsigned int j=0; j<6; j++) {
			result[i][j] = a[i] * b[j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the matrix resulting of the
/// multiplication of a column matrix by a row matrix
///
/// \param a [In]
/// \param b [In]
/// \param size [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyColumnByRow(double *a, double *b, double *result, unsigned int size) {
	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			result[i*size+j] = a[i] * b[j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of a 6x6
/// matrix by a scalar
///
/// \param scalar [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByScalar(double scalar, double result[6][6]) {
	for(unsigned int i=0; i<6; i++) {
		for(unsigned int j=0; j<6; j++) {
			result[i][j] *= scalar;
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of a 3x3
/// matrix by a scalar
///
/// \param scalar [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByScalar(double scalar, double result[3][3]) {
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			result[i][j] *= scalar;
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of a matrix
/// by a scalar
///
/// \param scalar [In]
/// \param size [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByScalar(double scalar, double *result, unsigned int size) {
	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			result[i*size+j] *= scalar;
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the sum of two 6x6 matrices
///
/// \param a [In] matrix a
/// \param b [In] matrix b
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 02/12/2012
///////////////////////////////////////////////////////////////
void ANMICMath::addMatrices(double a[6][6], double b[6][6], double result[6][6]) {
	for(unsigned int i=0; i<6; i++) {
		for(unsigned int j=0; j<6; j++) {
			result[i][j] = a[i][j] + b[i][j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the sum of two matrices
///
/// \param a [In] matrix a
/// \param b [In] matrix a
/// \param size [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
void ANMICMath::addMatrices(double *a, double *b, double *result, unsigned int size) {
	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			result[i*size+j] = a[i*size+j] + b[i*size+j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the sum of two 3x3 matrices
///
/// \param matrixA [In]
/// \param matrixB [In]
///
/// \param result [Out] Result matrix
///
/// \author vgil
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
void ANMICMath::addMatrices(double const matrixA[3][3], double const matrixB[3][3], double result[3][3]) {
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			result[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the addition of two 3x3 matrices,
/// adds "b" to "a" and returns the result in "a" matrix
///
/// \param b [In] matrix b
///
/// \param a [In/Out] Result matrix
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::addMatricesInPlace(double a[3][3], double b[3][3]) {
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			a[i][j] = a[i][j] + b[i][j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the addition of two matrices,
/// adds "b" to "a" and returns the result in "a" matrix
///
/// \param b [In] matrix b
/// \param size [In]
///
/// \param a [In/Out] Result matrix
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::addMatricesInPlace(double *a, double *b, unsigned int size) {
	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			a[i*size+j] = a[i*size+j] + b[i*size+j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the subtraction of two 3x3
/// matrices
///
/// \param matrixA [In]
/// \param matrixB [In]
///
/// \param result [Out] Result matrix
///
/// \author vgil
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
void ANMICMath::subtractMatrices(double const matrixA[3][3], double const matrixB[3][3], double result[3][3]) {
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			result[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the subtraction of two matrices
///
/// \param matrixA [In]
/// \param matrixB [In]
/// \param size [In]
///
/// \param result [Out] Result matrix
///
/// \author vgil
/// \author arincon
/// \date 01/09/2013
///////////////////////////////////////////////////////////////
void ANMICMath::subtractMatrices(double const *matrixA, double const *matrixB, double *result, unsigned int size) {
	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			result[i*size+j] = matrixA[i*size+j] - matrixB[i*size+j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of a one row matrix by a square matrix
///
/// \param row [In] Row matrix
/// \param matrix [In] Square matrix
/// \param size [In]
///
/// \param result [Out] Result matrix
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multitplyRowByMatrix(double *row, double *matrix, double *result, unsigned int size) {
	double value = 0;

	for(unsigned int i=0; i<size; i++) {
		for(unsigned int j=0; j<size; j++) {
			value += row[j] * matrix[i+j*size];
		}
		result[i] = value;
		value = 0;
	}
}

Point ANMICMath::multitplyPointByMatrix(Point& row, const double matrix[3][3]) {

	double result[3];
	vector<double> row_coords = row.getCoordinates();

	for(unsigned int i=0; i < 3; i++) {
		result[i] = 0;
		for(unsigned int j=0; j < 3; j++) {
			result[i] += row_coords[j] * matrix[i][j];
		}
	}

	return Point(result[0], result[1], result[2]);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of a one row matrix by one column matrix
///
/// \param a [In] Row matrix
/// \param b [In] Column matrix
/// \param size [In]
///
/// \return A double with the result of the multiplication
///
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
double ANMICMath::multiplyRowByColumn(double *a, double *b, unsigned int size) {
	double result = 0;

	for(unsigned int i=0; i<size; i++) {
		result += a[i] * b[i];
	}

	return result;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of two 3x3 matrices
///
/// \param first [In] matrix
/// \param second [In] matrix
///
/// \param result [Out]
///
/// \author vgil
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByMatrix(double first[3][3], double second[3][3], double result[3][3]) {
	int temp = 0;
	int a, b, c;
	for(a = 0; a < 3; a++) {
	   for(b = 0; b < 3; b++) {
		   for(c = 0; c < 3; c++) {
			   temp += first[b][c] * second[c][a];
		   }
		   result[b][a] = temp;
		   temp = 0;
	   }
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of two
/// square matrices
///
/// \param first [In] matrix
/// \param second [In] matrix
/// \param size [In]
///
/// \param result [Out]
///
/// \author vgil
/// \author arincon
/// \date 18/12/2012
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByMatrix(double *first, double *second, double *result, unsigned int size) {
	int temp = 0;
	unsigned int a, b, c;
	for(a = 0; a < size; a++) {
	   for(b = 0; b < size; b++) {
		   for(c = 0; c < size; c++) {
			   temp += first[b*size+c] * second[c*size+a];
		   }
		   result[b*size+a] = temp;
		   temp = 0;
	   }
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Multiplies two matrices of arbitrary size and compatible dimensions.
///
/// \param first [In] matrix as a vector of vectors (inner vector is the row)
/// \param second [In] matrix as a vector of vectors (inner vector is the row)
///
/// \param result [Out]
///
/// \author vgil
/// \date 03/03/2015
///////////////////////////////////////////////////////////////
void ANMICMath::multiplyMatrixByMatrix(std::vector<std::vector<double> >& first,
		std::vector<std::vector<double> >& second,
		std::vector<std::vector<double> >& out){

	out.clear();

	// check dimensions (first's row size must be equal to second's colum size)
	if (first[0].size() != second.size()){
		cout<<"ERROR matrix dimensions are not compatible (row x column) ("<<first[0].size()<<"x"<<first.size()
				<<" and "<<second[0].size()<<"x"<<second.size()<<")"<<endl;
		exit(-1);
	}

	for (unsigned int col_index=0; col_index < first.size(); ++col_index){
		// The row is first[col_index]
		vector<double> new_row;
		for (unsigned int col_2_index=0; col_2_index < second[0].size(); ++col_2_index){
			double dot = 0;
			for (unsigned int i=0; i < first[0].size(); ++i){ // row size
				dot += first[col_index][i]*second[i][col_2_index];
			}
			new_row.push_back(dot);
		}
		out.push_back(new_row);
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the multiplication of the
/// "I" matrix by the "e" vector
///
/// \param I [In] matrix
/// \param e [In] vector
///
/// \param result [Out]
///
/// \author vgil
/// \date 10/08/2012
///////////////////////////////////////////////////////////////
Point ANMICMath::multiplyIMatrixByEVector(double const I[3][3], double e[3]) {
	double value = 0;
	double result[3];

	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			value += I[i][j] * e[j];
		}
		result[i] = value;
		value = 0;
	}

	return Point(result[0], result[1], result[2]);
}

Point ANMICMath::multiplyIMatrixByEVector(double const I[3][3], Point& e){
	vector<double> coordinates = e.getCoordinates();
	return multiplyIMatrixByEVector(I, Utils::vectorToPointer(coordinates));
}


///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the inverse of a 3x3 matrix
///
/// \param a [In] matrix
///
/// \param result [Out]
///
/// \author arincon
/// \date 10/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::invertIMatrix(double a[3][3], double result[3][3]) {
	Point x0(a[0][0], a[1][0], a[2][0]);
	Point x1(a[0][1], a[1][1], a[2][1]);
	Point x2(a[0][2], a[1][2], a[2][2]);

	Point cross_x1_x2 = ANMICMath::crossProduct(x1, x2);
	double detA = ANMICMath::dotProduct(x0, cross_x1_x2);

	Point x1_x2 = ANMICMath::crossProduct(x1, x2);
	Point x2_x0 = ANMICMath::crossProduct(x2, x0);
	Point x0_x1 = ANMICMath::crossProduct(x0, x1);

	x1_x2.multiplyByScalar(1/detA);
	x2_x0.multiplyByScalar(1/detA);
	x0_x1.multiplyByScalar(1/detA);

	result[0][0] = x1_x2.getX();
	result[0][1] = x1_x2.getY();
	result[0][2] = x1_x2.getZ();

	result[1][0] = x2_x0.getX();
	result[1][1] = x2_x0.getY();
	result[1][2] = x2_x0.getZ();

	result[2][0] = x0_x1.getX();
	result[2][1] = x0_x1.getY();
	result[2][2] = x0_x1.getZ();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the inverse of a matrix
///
/// \param a [In] matrix
/// \param size [In]
///
/// \param result [Out]
///
/// \author arincon
/// \date 10/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::invertMatrix(double *a, double *result, unsigned int size) {
	Point x0(a[0*3+0], a[1*size+0], a[2*size+0]);
	Point x1(a[0*3+1], a[1*size+1], a[2*size+1]);
	Point x2(a[0*3+2], a[1*size+2], a[2*size+2]);

	Point cross_x1_x2 = ANMICMath::crossProduct(x1, x2);
	double detA = ANMICMath::dotProduct(x0, cross_x1_x2);

	Point x1_x2 = ANMICMath::crossProduct(x1, x2);
	Point x2_x0 = ANMICMath::crossProduct(x2, x0);
	Point x0_x1 = ANMICMath::crossProduct(x0, x1);

	x1_x2.multiplyByScalar(1/detA);
	x2_x0.multiplyByScalar(1/detA);
	x0_x1.multiplyByScalar(1/detA);

	result[0*3+0] = x1_x2.getX();
	result[0*3+1] = x1_x2.getY();
	result[0*3+2] = x1_x2.getZ();

	result[1*3+0] = x2_x0.getX();
	result[1*3+1] = x2_x0.getY();
	result[1*3+2] = x2_x0.getZ();

	result[2*3+0] = x0_x1.getX();
	result[2*3+1] = x0_x1.getY();
	result[2*3+2] = x0_x1.getZ();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the inverse of a matrix with arbitrary size
///
/// \param a [In] matrix
///
/// \param result [Out]
///
/// \author vgil
/// \date 17/03/2015
///////////////////////////////////////////////////////////////
void ANMICMath::invertMatrix(std::vector<std::vector<double> >& out){

}

void ANMICMath::invertMatrix(TriangularMatrix* inout){

	char uplo = 'U';
	int info = 0;
	int N = inout->size1();

	// Temporary copy
	TriangularMatrix tmp (*inout);

	dpptrf_(&uplo, &N, inout->data().begin(), &info);
	cout<< "DBG: INFO dpptrf_: "<<info<<endl;

	if(info == 0){
		dpptri_(&uplo, &N, inout->data().begin(), &info);
		cout<< "DBG: INFO dpotri_: "<<info<<endl;
	}
	else{
		info = 0;
		inout->swap(tmp);
		int ipv[N];
		double work[N];
		dsptrf_( &uplo, &N, inout->data().begin(), ipv, &info);
		cout<< "DBG: INFO dsptrf_: "<<info<<endl;
		dsptri_( &uplo, &N, inout->data().begin(), ipv, work, &info);
		cout<< "DBG: INFO dsptri_: "<<info<<endl;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the transpose of a 3x3 matrix
///
/// \param matrix [In]
///
/// \param result [Out]
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::transposeMatrix(double matrix[3][3], double result[3][3]) {
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			result[j][i] = matrix[i][j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function calculates the transpose of a matrix
///
/// \param matrix [In]
/// \param size [In]
///
/// \param result [Out]
///
/// \author arincon
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
void ANMICMath::transposeMatrix(double *matrix, double *result, unsigned int size) {
	for (unsigned int i = 0; i < size; ++i) {
		for (unsigned int j = 0; j < size; ++j) {
			result[j*size+i] = matrix[i*size+j];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function gets the 6x6 matrix Tab from the U matrix
/// in the positions indicated by (a, b)
///
/// \param UMatrix [In]
/// \param a [In]
/// \param b [In]
///
/// \param result [Out]
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::getABMatrix(std::vector< std::vector<double> > &UMatrix, double *result, int a, int b) {
	unsigned int i_ini = a * 6;
	unsigned int i_fi = i_ini + 6;
	unsigned int j_ini = b * 6;
	unsigned int j_fi = j_ini + 6;

	unsigned int n = 0;
	unsigned int m = 0;
	for(unsigned int i=i_ini; i<i_fi; i++) {
		for(unsigned int j=j_ini; j<j_fi; j++) {
			result[n*6+m] = UMatrix[i][j];
			m++;
		}
		n++;
		m = 0;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Adds the matrix Tab to the U matrix in the positions
/// indicated by (a, b)
///
/// \param UMatrix [In/Out]
///
/// \param Tab [In]
/// \param a [In]
/// \param b [In]
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::addMatrixToU(std::vector< std::vector<double> > &UMatrix, double *Tab, int a, int b) {
	unsigned int i_ini = a * 6;
	unsigned int i_fi = i_ini + 6;
	unsigned int j_ini = b * 6;
	unsigned int j_fi = j_ini + 6;

	unsigned int n = 0;
	unsigned int m = 0;
	for(unsigned int i=i_ini; i<i_fi; i++) {
		for(unsigned int j=j_ini; j<j_fi; j++) {
			UMatrix[i][j] += Tab[n*6+m];
			m++;
		}
		n++;
		m = 0;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Substracts the matrix Tab from the U matrix in the
/// positions indicated by (a, b)
///
/// \param UMatrix [In/Out]
///
/// \param Tab [In]
/// \param a [In]
/// \param b [In]
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
void ANMICMath::subsMatrixFromU(std::vector< std::vector<double> > &UMatrix, double *Tab, int a, int b) {
	unsigned int i_ini = a * 6;
	unsigned int i_fi = i_ini + 6;
	unsigned int j_ini = b * 6;
	unsigned int j_fi = j_ini + 6;

	unsigned int n = 0;
	unsigned int m = 0;
	for(unsigned int i=i_ini; i<i_fi; i++) {
		for(unsigned int j=j_ini; j<j_fi; j++) {
			UMatrix[i][j] -= Tab[n*6+m];
			m++;
		}
		n++;
		m = 0;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Calculates the cross products of two vectors
///
/// \param a [In] point a
/// \param b [In] point b
///
/// \return Point with the result of the cross product
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
Point ANMICMath::crossProduct(Point a, Point b) {
	return Point((a.getY()*b.getZ() - a.getZ()*b.getY()), (a.getZ()*b.getX() - a.getX()*b.getZ()), (a.getX()*b.getY() - a.getY()*b.getX()));
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Calculates the cross products of two vectors
///
/// \param a [In] point a
/// \param b [In] point b
///
/// \return A double with the result of the dot product
///
/// \author arincon
/// \date 22/11/2012
///////////////////////////////////////////////////////////////
double ANMICMath::dotProduct(Point a, Point b) {
	return a.getX()*b.getX() + a.getY()*b.getY() + a.getZ()*b.getZ();
}

void ANMICMath::printMatrix(double matrix[6][6]) {
	cout << "----------" << endl;
	for(unsigned int i=0; i<6; i++) {
		for(unsigned int j=0; j<6; j++) {
			cout <<setprecision(8)<<setw(15)<< matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << "----------"  << endl;
}

void ANMICMath::printMatrix(double matrix[3][3]) {
	cout << "----------" << endl;
	for(unsigned int i=0; i<3; i++) {
		for(unsigned int j=0; j<3; j++) {
			cout <<setprecision(18)<<setw(25)<< matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << "----------"  << endl;
}

void ANMICMath::printMatrix(std::vector<std::vector<double> >& matrix){
	cout << "----------" << endl;
	for(unsigned int i=0; i<matrix.size(); i++) {
		for(unsigned int j=0; j<matrix[i].size(); j++) {
			cout <<setprecision(18)<<setw(25)<< matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << "----------"  << endl;
}

void ANMICMath::printMatrix(double *matrix, int length) {
	cout << "----------"  << endl;
	for(int i=0; i<length; i++) {
		for(int j=0; j<length; j++) {
			cout <<setprecision(18)<<setw(2)<< matrix[i*length+j] << " ";
		}
		cout << endl;
	}
	cout << "----------"  << endl;
}

void ANMICMath::printVector(double *vector, int length) {
	for(int i=0; i<length; i++) {
		cout <<setprecision(18)<<setw(2)<< vector[i] << " ";
	}
	cout << endl;
}

void ANMICMath::convertMatrixToArray(double **matrix, double *array, int rows, int columns) {
	for(int i=0; i<rows; i++) {
		for(int j=0; j<columns; j++) {
			array[i*columns+j] = matrix[i][j];
		}
	}
}

Unit * ANMICMath::getDummyUnit(double x, double y, double z) {
	Atom *atom = new Atom();
	double *coordinates = new double[3];
	coordinates[0] = x;
	coordinates[1] = y;
	coordinates[2] = z;
	atom->setCoordinates(coordinates);
	atom->setCoordinatesIndexes(0, 1, 2);
	atom->setMass(1);

	vector<Atom*> dummy_atoms;
	dummy_atoms.push_back(atom);
	Unit* u = new Unit(dummy_atoms, dummy_atoms, NULL, NULL, NULL, NULL,"");

	return u;
}


void ANMICMath::transpose(std::vector<std::vector<double> >& in,
			std::vector<std::vector<double> >& out){
	out.clear();

	for(unsigned int i = 0; i < in[0].size(); ++i){
		vector<double> new_row;
		for(unsigned int j = 0; j < in.size(); ++j){
			new_row.push_back(in[j][i]);
		}
		out.push_back(new_row);
	}
}
