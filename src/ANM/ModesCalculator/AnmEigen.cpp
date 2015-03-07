/////////////////////////////////////////////////////////////////////////////
/// AnmEigen.cpp
///
/// Implementation of AnmEigen class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author mrivero
/// \author atarraco
/// \author xoro
/// \date 03/09/2012
/////////////////////////////////////////////////////////////////////////////

#include "AnmEigen.h"

#include "../../Tools/Utils.h"
#include "../Tools/AnmNormalizer.h"

#include "../../Tools/Math/MathTools.h"

#include <sstream>
#include <iomanip>
#include "../../PELE/PeleTasks/Output/LogUtils.h"

#include "../../Tools/vectorTools.h"

#include <sstream>

using namespace std;

AnmEigen::AnmEigen()
{}

AnmEigen::~AnmEigen()
{}

///////////////////////////////////////////////////////////////
/// \remarks
/// It averages the eigenvectors using a given list of weights
///
/// \param weights [In] List of weights
///
/// \param averageEigenVector [Out] Resulting average eigenVector
///
/// \author atarraco
/// \author mrivero
/// \date 05/09/2012
///////////////////////////////////////////////////////////////
void AnmEigen::computeAverageEigenVector(const std::vector<double> & weights, std::vector<double> & averageEigenVector)
{
	unsigned int number_of_elements = vectors[0].size();
	unsigned numberOfModes = values.size();

	averageEigenVector.clear();
	averageEigenVector.resize(number_of_elements);

	for (unsigned int i = 0; i < numberOfModes; ++i) {
		for(unsigned int j=0;j < number_of_elements; ++j) {
			averageEigenVector[j] +=  vectors[i][j] * weights[i];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It initializes the eigenvectors and eigenvalues.
///
/// \param eigenValues [In] Array of eigenvalues
/// \param eigenVectors [In] Array of eigenvectors
/// \param numberOfModes [In] Number of modes
/// \param numberOfNodes [In] Number of nodes
///
/// \author atarraco
/// \author mrivero
/// \author vgil
/// \date 05/09/2012
///////////////////////////////////////////////////////////////
void AnmEigen::initialize(const double * const eigenValues, const double * const eigenVectors,
							unsigned int numberOfModes, unsigned int numberOfNodes, bool usingCartesian)
{
	vectors.clear();
	values.clear();

	Utils::convertArrayOfDoubleToVector(eigenValues, values, numberOfModes);
	unsigned int eigenvector_length;
	if(usingCartesian){
		eigenvector_length = numberOfNodes*3;
	}
	else{
		// In this case it is not the number of nodes (which would be numberOfNodes+1) but
		// the real length of the vector.
		eigenvector_length = numberOfNodes;
	}

	Utils::convertArrayOfDoubleToVectorOfVectors(eigenVectors, vectors, numberOfModes, eigenvector_length);

	this->numberOfModes = numberOfModes;
	this->numberOfNodes = numberOfNodes;
}


///////////////////////////////////////////////////////////////
/// \remarks
/// Initialization with vectors. NO CHECKS PERFORMED.
///
/// \param eigenValues [In] vector containing the eigenvalues
/// \param eigenVectors [In] vector containing the eigenvectors
///
/// \author vgil
/// \date 06/11/2014
///////////////////////////////////////////////////////////////

void AnmEigen::initialize(std::vector<double>& eigenValues,
						std::vector<std::vector<double> >& eigenVectors,
						bool usingCartesian){

	vectors.clear();
	values.clear();

	VectorTools::copy(values, eigenValues);
	for(unsigned i = 0; i < eigenVectors.size(); ++i){
		vectors.push_back(eigenVectors[i]);
	}

	this->numberOfModes = values.size();

	if(usingCartesian){
		this->numberOfNodes = vectors[0].size()/3;
	}
	else{
		this->numberOfNodes = vectors[0].size();
	}

	cout<<"DBG: AnmEigen initialized from vector<double>. "<< this->numberOfModes <<" modes read."<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////
/// \remarks
///	The eigenVector is normalized using the inverse of the highest direction module
///	For example:
///	If the eigenvector is {1,3,4,5,6,7,2,8,2}, the three directions will be:
///  direction1: (1,3,4) -> norma = 5,099
///  direction2: (5,6,7) -> norma = 10,488
///  direction3: (2,8,2) -> norma = 8,485
/// Then we choose the normalization factor 1/10,488 = 0,0954198,
/// because we want that the direction module is bigger than one.
///	Finally each value of the eigenvector is multiplied by this factor
/// Eigenvector after normalization:
/// {0.0953463,0.286039,0.381385,0.476731,0.572078,0.667424,0.190693,0.76277,0.190693}
///
/// \author atarraco
/// \date 08/07/2012
///////////////////////////////////////////////////////////////////////////////////
void AnmEigen::normalizeByInverseLargestNorm()
{
	unsigned int numberOfModes = getNumberOfModes();

	for(unsigned int i=0; i<numberOfModes; ++i) {
		AnmNormalizer::normalizeByInverseLargestNorm(vectors[i]);
	}
	cout<<"DBG: Normalized by inverse largest norm"<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////
/// \remarks
/// It normalizes the eigenvectors using the largest inverse norm divided by the
/// root of the hessian constant multiplied by the corresponding eigenvalues
/// For instance:
/// If the eigenvector is {1,3,4,5,6,7,2,8,2} with an eigenvalue equal to 4
/// it computes the of the three directions:
///  norma = sqrt(208) = 14,4222
/// and uses as normalization factor the inverse of the norm divided by the hessian
/// constant multiplied by the eigenvalue.
///
/// Assuming a hessian constant equal to 1.667, the inverse norm would be:
///  inverseNorm = 1/14,4222 = 0.0693375245
/// And then:
///  normalizationFactor = inverseNorm / sqrt(1.667*4) = 0.0268516
///
/// Finally, every element of the eigenvector is multiplied by normalization factor to
/// get the normalized eigenvector:
///  {0.0268516,0.0805549,0.107406,0.134258,0.16111,0.187961,0.0537032,0.214813,0.0537032}
///
/// \param constantForHessian [In] Constant for Hessian
///
/// \author atarraco
/// \date 08/07/2012
///////////////////////////////////////////////////////////////////////////////////
void AnmEigen::normalizeByInverseNormWithThermal(double constantForHessian){

	unsigned int numberOfModes = getNumberOfModes();

	for(unsigned int i=0; i<numberOfModes; ++i) {
		AnmNormalizer::normalizeByInverseLargestNormWithThermalScaling(vectors[i], values[i], constantForHessian);
	}
}

std::vector<double> & AnmEigen::getEigenVectorOfMode(unsigned int mode)
{
	return vectors.at(mode);
}

double AnmEigen::getEigenValueOfMode(unsigned int mode)
{
	return values.at(mode);
}

unsigned int AnmEigen::getNumberOfModes() const
{
	return values.size();
}

unsigned int AnmEigen::getEigenVectorsDimension() const
{
	return vectors.at(0).size();
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It returns a string showing the list of eigenvalues and
/// eigenVectors
///
/// \return A string showing the list of eigenvalues and eigenVectors
///
/// \author mrivero
/// \date 18/09/2012
///////////////////////////////////////////////////////////////
std::string AnmEigen::toString() const
{
	 ostringstream oss (ostringstream::out);

	 unsigned int numberOfModes = getNumberOfModes();
	 oss<<"Number of modes "<<numberOfModes<<endl;
	 for(unsigned int i=0; i<numberOfModes; ++i)
	 {
		 oss << setprecision(10) << values[i] <<endl;
	 }
	 oss << endl;

	 unsigned int dimension = getEigenVectorsDimension();
	 oss<<"Eigen vectors dimension "<<dimension<<endl;
	 for(unsigned int i=0; i<numberOfModes; ++i)
	 {
		 oss << setprecision(10) << LogUtils::showVector(vectors[i]);
		 oss << endl;
	 }
	 return oss.str();
}

AnmEigen& AnmEigen::operator=(const AnmEigen& other){

	VectorTools::copy<double>(values, other.values);

	unsigned int numberOfModes = other.vectors.size();
	for (unsigned int i = 0; i < numberOfModes; i++){
		vector<double> newMode;
		vectors.push_back(newMode);
		VectorTools::copy<double>(vectors[i], other.vectors[i]);
	}

	return *this;
}

void AnmEigen::normalizeByLargestValue()
{

	cout<<"DBG: Normalizing modes (Modes: "<<getNumberOfModes()<<", Nodes:"<<numberOfNodes<<")"<<endl;
	for(unsigned int i=0; i < getNumberOfModes(); ++i) {
		AnmNormalizer::normalizeByLargestValue(vectors[i]);
	}
}
