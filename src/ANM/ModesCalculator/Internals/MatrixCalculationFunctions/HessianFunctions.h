/////////////////////////////////////////////////////////////////////////////
/// \file ANMFunctions.h
///
/// \brief Functions used in the computation of the hessian matrix
///
/// \author arincon
/// \author vgil
/// \date 17/11/2012
/////////////////////////////////////////////////////////////////////////////

#include <vector>
#include "TriangularMatrices.h"
class Point;
class Unit;
class ElasticConstantCalculator;

namespace ANMICHessianCalculator{

	void calculateDij(double k, Point& r_i, Point& r_j, double *result);

	void calculateTab(Unit* unit_a, Unit* unit_b, double sq_cutoff,  ElasticConstantCalculator* force_constant_calculator, double *result, bool skip_OXT);

	std::vector< std::vector<double> > calculateU(double sq_cutoff, ElasticConstantCalculator* force_constant_calculator, std::vector<Unit*>& units, bool skip_OXT);

	void calculateRab(int a, int b, std::vector< std::vector<double> > &UMatrix, double *Rab);

	double calculateHab(Point* e_a,Point* e_b, Point* r_a,Point* r_b, int a_dihedral, int b_dihedral, std::vector< std::vector<double> >& U);

	TriangularMatrix* calculateH(std::vector<Unit*>& units, std::vector< std::vector<double> > &U);

	void modifyHessianWithExtraTorsion(TriangularMatrix* H);

};
