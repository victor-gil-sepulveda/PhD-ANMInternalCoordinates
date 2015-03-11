/////////////////////////////////////////////////////////////////////////////
/// ANMICMath.cpp
///
/// Implementation of KinecticMatrixFunctions
///
/// \author vgil
/// \author arincon
/// \date 11/08/2013
/////////////////////////////////////////////////////////////////////////////

#include "KineticMatrixFunctions.h"

#include "../CoarseGrainModel/Unit.h"
#include "../../../../Molecules/Atom.h"
#include "../../../../Tools/vectorTools.h"
#include "ANMICMath.h"
#include "../../../../Tools/Utils.h"
#include "TriangularMatrices.h"
#include "../../../../PELE/PeleTasks/Sensors/Metrics/Tools/CenterOfMass.h"
#include <vector>
#include <cstring>
#include <utility>
#include <string>
#include "../../../../Molecules/AtomNames.h"
#include "../../../../Molecules/AtomSetsTree/Links/LinkNames.h"
#include "../../../../Tools/Math/Point.h"
#include "../CoarseGrainModel/UnitTools.h"

using namespace std;



///////////////////////////////////////////////////////////////
/// \remarks
///	Calculates the inertia tensor using different algorithms.
/// \param p [In]
/// \param result [In]
///
/// \author vgil
///
/// \date 7/10/2014
///////////////////////////////////////////////////////////////
void ANMICKineticMatrixCalculator::calculateI(	double I[3][3],
														std::vector<Unit*>& units,
														const std::pair<unsigned int, unsigned int>& range,
														ICalcType calc_type){
	// I = 0
	memset(I, 0, 3*3*sizeof(double));

	switch(calc_type){

		case BRAUN:
		{
			// Ec. 33 Braun et al. 1984
			double C[3][3];
			for (unsigned int i = range.first; i <= range.second; ++i){
				vector<Atom*>& atoms = units[i]->atoms;
				for (unsigned int j = 0; j < atoms.size(); ++j){
					double mass = atoms[j]->getMass();
					calculateC(atoms[j]->toPoint(),C);
					ANMICMath::multiplyMatrixByScalar(mass,C);
					ANMICMath::addMatrices(I,C,I);
				}
			}
		}
		break;

		case INMA:
		{
			// Ec 22 Noguti & Go 1983
			// Based in INMA Code by Jose Ramon Lopez Sanchez
			for (unsigned int i = range.first; i <= range.second; ++i){
				vector<Atom*>& atoms = units[i]->atoms;
				for (unsigned int j = 0; j < atoms.size(); ++j){
					double mass = atoms[j]->getMass();
					double coords[3];
					coords[0] = atoms[j]->getX();
					coords[1] = atoms[j]->getY();
					coords[2] = atoms[j]->getZ();

					double tmp1 = coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
					for(unsigned int m = 0; m < 3; ++m){
						for(unsigned int n = 0; n < 3; ++n){
							double tmp2 = -coords[m] * coords[n];
							if ( m == n) tmp2 += tmp1;
							I[m][n] += mass*tmp2;
						}
					}
				}
			}

		}
		break;

		case INMA_NO_OXT:
		{
			// Ec 22 Noguti & Go 1983
			// Based in INMA Code by Jose Ramon Lopez Sanchez
			for (unsigned int i = range.first; i <= range.second; ++i){
				vector<Atom*>& atoms = units[i]->atoms;
				for (unsigned int j = 0; j < atoms.size(); ++j){
					if (atoms[j]->name != " OXT"){
						double mass = atoms[j]->getMass();
						double coords[3];
						coords[0] = atoms[j]->getX();
						coords[1] = atoms[j]->getY();
						coords[2] = atoms[j]->getZ();

						double tmp1 = coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
						for(unsigned int m = 0; m < 3; ++m){
							for(unsigned int n = 0; n < 3; ++n){
								double tmp2 = -coords[m] * coords[n];
								if ( m == n) tmp2 += tmp1;
								I[m][n] += mass*tmp2;
							}
						}
					}
				}
			}

		}
		break;

		case LU:
		{
			// Lu et al. 2006 Eq. 11
			double P[3][3];
			double Pt[3][3];
			double PtP[3][3];

			// Sum over atoms
			for (unsigned int i = range.first; i <= range.second; ++i){
				vector<Atom*>& atoms = units[i]->atoms;
				for (unsigned int j = 0; j < atoms.size(); ++j){
					if (atoms[j]->name != " OXT"){
						double mass = atoms[j]->getMass();
						calculateP(atoms[j]->toPoint(), P);
						ANMICMath::transposeMatrix(P, Pt);
						ANMICMath::multiplyMatrixByMatrix(Pt,P, PtP);
						ANMICMath::multiplyMatrixByScalar(mass, PtP);
						ANMICMath::addMatricesInPlace(I, PtP);
					}
				}
			}

		}
		break;
	}

}

///////////////////////////////////////////////////////////////
/// \remarks
/// ...
///
/// \param p [In] Point p
///
/// \param result [Out]
///
/// \author vgil
/// \author arincon
///
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
// Braun et al. 1984 Eq. 33, analytic form
void ANMICKineticMatrixCalculator::calculateC(const Point& p, double result[3][3]){
	vector<double> y = p.getCoordinates();

	result[0][0]  = y[1]*y[1] + y[2]*y[2];
	result[0][1]  = -y[0]*y[1];
	result[0][2]  = -y[0]*y[2];

	result[1][0]  = -y[1]*y[0];
	result[1][1]  = y[0]*y[0] + y[2]*y[2];
	result[1][2]  = -y[1]*y[2];

	result[2][0]  = -y[2]*y[0];
	result[2][1]  = -y[2]*y[1];
	result[2][2]  = y[0]*y[0] + y[1]*y[1];
}

///////////////////////////////////////////////////////////////
/// \remarks
///
/// \param p [In]
/// \param result [In]
///
/// \author vgil
/// \author arincon
///
/// \date 11/08/2013
///////////////////////////////////////////////////////////////
// Lu et al. 2006 Eq. 11
void ANMICKineticMatrixCalculator::calculateP(const Point& p, double result[3][3]) {
	result[0][0] = 0;
	result[0][1] = -p.getZ();
	result[0][2] = p.getY();

	result[1][0] = p.getZ();
	result[1][1] = 0;
	result[1][2] = -p.getX();

	result[2][0] = -p.getY();
	result[2][1] = p.getX();
	result[2][2] = 0;
}

///////////////////////////////////////////////////////////////
/// \remarks
///
/// \param units [In]
/// \param range [In]
///
/// \return double
///
/// \author vgil
/// \author arincon
///
/// \date 15/08/2013
///////////////////////////////////////////////////////////////
double ANMICKineticMatrixCalculator::calculateM(vector<Unit*>& units, const pair<int,int>& range, bool skip_OXT){
	double M = 0;

	if (!skip_OXT){
		for (int i = range.first; i <= range.second; ++i) {
			M += units[i]->getMass();
		}
	}
	else{
		for (int i = range.first; i <= range.second; ++i) {
			vector<Atom*> unit_atoms = units[i]->atoms;
			for (unsigned j = 0; j < unit_atoms.size(); ++j){
				Atom* atom = unit_atoms[j];
				if(atom->name != AtomNames::OXT){
					M += unit_atoms[j]->getMass();
				}
			}
		}
	}

	return M;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Calculates the center of mass of a rigid macrounit
/// (rigid unit formed of other units)
///
/// \param units [In]
/// \param range [In]
///
/// \param CenterOfMass [Out]
///
/// \author vgil
///
/// \date 15/08/2013
///////////////////////////////////////////////////////////////
CenterOfMass ANMICKineticMatrixCalculator::calculateMobileBodyCOM(vector<Unit*>& units, const pair<int,int>& range, bool skip_OXT){
	vector<Atom*> macrounit_atoms;
	if (!skip_OXT){
		for (int i = range.first; i <= range.second; ++i) {

			VectorTools::add(macrounit_atoms, units[i]->atoms);
		}
	}
	else{
		for (int i = range.first; i <= range.second; ++i) {
			for(unsigned int j = 0;j < units[i]->atoms.size(); ++j){
					Atom* atom = units[i]->atoms[j];
					if (atom->name != AtomNames::OXT){
						macrounit_atoms.push_back(atom);
					}
			}
		}
	}
	return CenterOfMass(CenterOfMass::compute(macrounit_atoms));
}

//----------------------------------------------
// Precalculations (to do). For each dihedral, the right body can be precalculated using left:
// M3
// -------
// M3 = M - M1
// -------
// COM3
// -------
// Point com3 = Point(-M1*com1.getX()/M3, -M1*com1.getY()/M3, -M1*com1.getZ()/M3);
// equal to : CenterOfMass com3 = calculateMacrounitCOM(units, pair<int,int>(alpha_dihedral+1,units.size()-1));
// -------
// I3
// -------
// I3 = I -I1
// ANMICMath::subtractMatrices(I,I1,I3);
//----------------------------------------------
///////////////////////////////////////////////////////////////
/// \remarks
/// Calculates the center of mass of a rigid macrounit
/// (rigid unit formed of other units)
///
/// \param alpha_dihedral [In]
/// \param beta_dihedral [In]
/// \param rigid_zone_starting_unit [In]
/// \param rigid_zone_ending_unit [In]
/// \param r_a [In]
/// \param r_b [In]
///
/// \return double
///
/// \author vgil
/// \date 15/08/2013
///////////////////////////////////////////////////////////////
double ANMICKineticMatrixCalculator::calculateKab(
		int alpha_dihedral, int beta_dihedral,
		Point* r_a, Point* r_b,
		Point* e_a, Point* e_b,
		vector<Unit*>& units,
		double M,
		double const I[3][3],
		double const I_inv[3][3],
		ICalcType tensor_calc_type) {

	bool skip_OXT = tensor_calc_type == INMA_NO_OXT? true: false;

	CenterOfMass com1 = calculateMobileBodyCOM(units, pair<int,int>(0,alpha_dihedral), skip_OXT);
	CenterOfMass com3 = calculateMobileBodyCOM(units, pair<int,int>(beta_dihedral+1,units.size()-1), skip_OXT);

	double M1 = com1.getMass();
	//TODO: M3 = M-M1
	double M3 = com3.getMass();

	Point r_1_0 = Point(com1);
	Point r_3_0 = Point(com3);

	// Calculation of first term
	Point sub_r_a_r_1_0 = Point::subtract(*r_a, r_1_0); // (ra - r10)
	Point sub_r_b_r_3_0 = Point::subtract(*r_b, r_3_0); // (rb - r30)

	Point cross_e_a_sub_r_a_r_1_0 = ANMICMath::crossProduct(*e_a, sub_r_a_r_1_0); // ea x (ra - r10)
	Point cross_e_b_sub_r_b_r_3_0 = ANMICMath::crossProduct(*e_b, sub_r_b_r_3_0); // eb x (rb - r30)

	double dot_product = ANMICMath::dotProduct(cross_e_a_sub_r_a_r_1_0, cross_e_b_sub_r_b_r_3_0);// (ea x (ra - r10)) . (eb x (rb - r30))
	double term1 = dot_product;// * M1 * M3 / M; // (M1 M3/ M) * (ea x (ra - r10)) . (eb x (rb - r30))

	// Calculation of second term
	double I1[3][3], I3[3][3];
	calculateI(I1,units, pair<int,int>(0,alpha_dihedral), tensor_calc_type);

	//TODO: I3 = I - I1
	calculateI(I3,units, pair<int,int>(beta_dihedral+1,units.size()-1), tensor_calc_type);

	Point M1_r_1_0 = Point::multiplyByScalar(r_1_0,M1); // r10 = M1*r10
	Point cross_e_a_r_a = ANMICMath::crossProduct(*e_a, *r_a);
	Point crossEM1 = ANMICMath::crossProduct(M1_r_1_0, cross_e_a_r_a);

	Point M3_r_3_0 = Point::multiplyByScalar(r_3_0,M3); // r30 = M3*r30
	Point cross_e_b_r_b = ANMICMath::crossProduct(*e_b, *r_b);
	Point crossEM3 = ANMICMath::crossProduct(M3_r_3_0, cross_e_b_r_b);

	vector<double> e_a_coords = e_a->getCoordinates();
	vector<double> e_b_coords = e_b->getCoordinates();
	Point I1e_a = ANMICMath::multiplyIMatrixByEVector(I1, Utils::vectorToPointer<double>(e_a_coords));
	Point I3e_b = ANMICMath::multiplyIMatrixByEVector(I3, Utils::vectorToPointer<double>(e_b_coords));

	Point resta1P = Point::subtract(crossEM1, I1e_a);
	Point resta3P = Point::subtract(crossEM3, I3e_b);

	vector<double> resta1P_coords = resta1P.getCoordinates();
	vector<double> resta3P_coords = resta3P.getCoordinates();

	// Multiply resta1P_coords * I_inv
	vector<double> Resta1PImult_coords;
	Resta1PImult_coords.push_back(Point::dotProduct(resta1P, Point(I_inv[0][0], I_inv[0][1], I_inv[0][2])));
	Resta1PImult_coords.push_back(Point::dotProduct(resta1P, Point(I_inv[1][0], I_inv[1][1], I_inv[1][2])));
	Resta1PImult_coords.push_back(Point::dotProduct(resta1P, Point(I_inv[2][0], I_inv[2][1], I_inv[2][2])));

	double term2 =ANMICMath::multiplyRowByColumn(Utils::vectorToPointer<double>(Resta1PImult_coords),
			Utils::vectorToPointer<double>(resta3P_coords), 3);

//	cout<<"** i= "<< alpha_dihedral<<  "  j= "<< beta_dihedral
//				<< "  M1= " << M1<<"  M3= "<<M3
//				<<"  mtot= "<< M <<" temp1= "<< term1<<" term2= " << term2
//				<< " "<<(M1*M3/M) * term1 + term2 << endl;
//	cout<<"T " << Resta1PImult_coords[0] << " "<<Resta1PImult_coords[1]<<" "<<Resta1PImult_coords[2] <<endl;
//	cout<<"emMY " << resta3P_coords[0] << " "<<resta3P_coords[1]<<" "<<resta3P_coords[2] <<endl;

	return (M1*M3/M) * term1 + term2;

}


///////////////////////////////////////////////////////////////
/// \remarks
///
///
/// \param units [In]
///
/// \return double 2D matrix
///
/// \author vgil
/// \date 15/08/2013
///////////////////////////////////////////////////////////////
TriangularMatrix* ANMICKineticMatrixCalculator::calculateK( vector<Unit*>& units,
																		ICalcType tensor_calc_type){

	double I[3][3], I_inv[3][3];

	bool skip_OXT = tensor_calc_type == INMA_NO_OXT? true: false;

	double M = calculateM(units, pair<int,int>(0, units.size()-1), skip_OXT); // sum m_i for i in [0,number_of_units-1]

	calculateI(I, units, pair<int,int>(0, units.size()-1), tensor_calc_type); // calculate I for all units i in [0,number_of_units-1]
	ANMICMath::invertIMatrix(I,I_inv);

	int number_of_dihedrals = units.size() - 1;

	TriangularMatrix* K = new TriangularMatrix(number_of_dihedrals, number_of_dihedrals);

	for(int alpha_dihedral = 0; alpha_dihedral < number_of_dihedrals; ++alpha_dihedral){
		for(int beta_dihedral = alpha_dihedral; beta_dihedral < number_of_dihedrals; ++beta_dihedral){
			(*K)(alpha_dihedral, beta_dihedral) = calculateKab(
					alpha_dihedral, beta_dihedral,
					units[alpha_dihedral]->r_right, units[beta_dihedral]->r_right,
					units[alpha_dihedral]->e_right, units[beta_dihedral]->e_right,
					units, M, I, I_inv,
					tensor_calc_type);

		}
	}

	return K;
}

// based in ... paper
void ANMICKineticMatrixCalculator::Jacobi(vector<Unit*>& units, vector< vector<double> >& J){

	unsigned int number_of_torsions = units.size() - 1;

	// Total mass
	bool doNotSkipOXT = false;
	double M = calculateM(units, pair<int,int>(0, units.size()-1), doNotSkipOXT);

	// Intertia tensor
	double I[3][3], I_inv[3][3];
	calculateI(I, units, pair<int,int>(0, number_of_torsions), INMA);
	ANMICMath::invertIMatrix(I,I_inv);

	// Atoms
	vector<Atom*> atoms;
	bool onlyHeavyAtoms = true;
	UnitTools::getAllAtomsFromUnits(units, atoms, onlyHeavyAtoms);
	unsigned int number_of_atoms = atoms.size();

	// Give J the correct size
	J.clear();
	J.resize(number_of_atoms*3);
	for (unsigned int i =0; i <J.size(); ++i){
		J[i].resize(number_of_torsions);
	}

	for (unsigned int alpha = 0; alpha < number_of_torsions; ++alpha){
		// I_a
		double Ia[3][3];
		calculateI(Ia, units, pair<int,int>(0, alpha), INMA);

		// r1_0
		CenterOfMass r_a_0_c = calculateMobileBodyCOM(units, pair<int,int>(0,alpha), doNotSkipOXT);
		Point ra0(r_a_0_c.getX(), r_a_0_c.getY(), r_a_0_c.getZ());
		double Ma = r_a_0_c.getMass();

		Point* ea = units[alpha]->e_right;
		Point* ra = units[alpha]->r_right;

		Point ra0_ra = Point::subtract(ra0, *ra);
		Point eaxra0_ra = ANMICMath::crossProduct(*ea, ra0_ra);

		// ta calculation
		Point t = Point::multiplyByScalar(eaxra0_ra, -Ma/M);

		// Aa calculation
		Point Iaea = ANMICMath::multiplyIMatrixByEVector(Ia, *ea);
		Point ra0eaxra0_ra = ANMICMath::crossProduct(ra0,eaxra0_ra);
		Point _Mara0eaxra0_ra = Point::multiplyByScalar(ra0eaxra0_ra, -Ma);
		Point p = Point::subtract(_Mara0eaxra0_ra, Iaea);
		Point A = ANMICMath::multiplyIMatrixByEVector(I_inv, p);

		// When getting the atoms we assume same ordering every time
		bool onlyHeavyAtoms = true;
		unsigned int num_left_atoms = UnitTools::getNumberOfAtomsOfUnitRange(units, 0, alpha, onlyHeavyAtoms);

		// Calculate Jia for each atom (i)
		for (unsigned int i = 0; i < number_of_atoms; ++i){
			Point Jia;
			Point ri = atoms[i]->toPoint();
			if (i < (unsigned int) num_left_atoms){
				Point ri_ra = Point::subtract(ri,*ra);
				Point eaxri_ra = ANMICMath::crossProduct(*ea, ri_ra);
				Jia = Point::add(eaxri_ra, Point::add(ANMICMath::crossProduct(A,ri),t));
			}
			else{
				Jia = Point::add(ANMICMath::crossProduct(A,ri),t);
			}

			unsigned int offset = i*3;
			J[offset][alpha]   = Jia.getX();
			J[offset+1][alpha] = Jia.getY();
			J[offset+2][alpha] = Jia.getZ();
		}
	}
}


// Uses dr dq descriptions of JRLB thesis
void ANMICKineticMatrixCalculator::Jacobi2(std::vector<Unit*>& units, std::vector< std::vector<double> >& J){
	J.clear();

	unsigned int number_of_torsions = units.size()-1;
	double I[3][3], I_inv[3][3];
	calculateI(I, units, pair<int,int>(0, number_of_torsions), INMA);
	ANMICMath::invertIMatrix(I,I_inv);

	vector<Atom*> atoms;
	bool onlyHeavyAtoms = true;
	UnitTools::getAllAtomsFromUnits(units, atoms, onlyHeavyAtoms);

	unsigned int number_of_atoms = atoms.size();
	J.resize(number_of_atoms*3);
	for (unsigned int i =0; i <J.size(); ++i){
		J[i].resize(number_of_torsions);
	}

//  .----- d --->
//	|
//	|
//	n
//	|
//  \/
//

	std::vector<Point> terms1l,terms2l, terms1r, terms2r;
	for (unsigned int alpha = 0; alpha < number_of_torsions; ++alpha){
			dr1dq(units, alpha, I_inv, terms1l, terms2l);
			dr2dq(units, alpha, I_inv, terms1r, terms2r);
	}

	for (unsigned int alpha = 0; alpha< number_of_torsions; ++alpha){
		bool onlyHeavyAtoms = true;
		unsigned int num_left_atoms = UnitTools::getNumberOfAtomsOfUnitRange(units, 0, alpha, onlyHeavyAtoms);

		for (unsigned int i = 0; i < number_of_atoms; ++i){
			Point drdq;
			Point ri = atoms[i]->toPoint();

			if (i < (unsigned int) num_left_atoms){
				drdq = Point::subtract(terms1l[alpha], ANMICMath::crossProduct(terms2l[alpha],ri));
			}
			else{
				drdq = Point::subtract(ANMICMath::crossProduct(terms2r[alpha],ri), terms1r[alpha]);
			}

			int offset = i*3;
			J[offset][alpha]   = drdq.getX();
			J[offset+1][alpha] = drdq.getY();
			J[offset+2][alpha] = drdq.getZ();
		}
	}

}
/////////////////
/// Functions used for the CI -> CC  modes conversion (rewriting of K functions).
///
///////////////
// Eq 4.1.12 de JR, calculates "left" contributions
void  ANMICKineticMatrixCalculator::dr1dq(vector<Unit*>& units,
												unsigned int alpha,
												double const I_inv[3][3],
												std::vector<Point>& term1_v,
												std::vector<Point>& term2_v){

	Unit* unit = units[alpha];
	double M = calculateM(units, pair<int,int>(0, units.size()-1), false);
	CenterOfMass ri_0 = calculateMobileBodyCOM(units, pair<int,int>(0,alpha), false);
	double M1 = ri_0.getMass();
	double M2 =  M - M1;

	Point* r_a = unit->r_right;
	Point* e_a = unit->e_right;
	vector<double> r_a_coords = unit->r_right->getCoordinates();
	vector<double> e_a_coords = unit->e_right->getCoordinates();
	double I2[3][3];
	calculateI(I2, units, pair<int,int>(alpha+1, units.size()-1), INMA);

	// e_a x{(M2/M)*ra+(M1/M)*ri0]}
	Point M2Mr_a = Point::multiplyByScalar(*r_a, M2/M);
	Point M1Mri_0 = Point::multiplyByScalar(ri_0, M1/M);
	Point term1 = Point::add(M2Mr_a, M1Mri_0);
	Point e_a_x_term1 = ANMICMath::crossProduct(*e_a, term1);

	//{[M1ri_0 x(eaxra)+I2ea]I^-1}xr1
	Point Mri_0 = Point::multiplyByScalar(ri_0, M1);
	Point ea_ra = ANMICMath::crossProduct(*e_a, *r_a);
	Point Mri_0_aa_x_r_a = ANMICMath::crossProduct(Mri_0, ea_ra);

	Point I2ea = ANMICMath::multiplyIMatrixByEVector(I2, Utils::vectorToPointer<double>(e_a_coords));
	Point term2 =  Point::add(Mri_0_aa_x_r_a, I2ea);

	Point term2Iinv = ANMICMath::multitplyPointByMatrix( term2, I_inv);

	term1_v.push_back(e_a_x_term1);
	term2_v.push_back(term2Iinv);
}

// Eq 4.1.13 de JR, calculates "right" contributions
void  ANMICKineticMatrixCalculator::dr2dq(vector<Unit*>& units,
												unsigned int beta,
												double const I_inv[3][3],
												std::vector<Point>& term1_v,
												std::vector<Point>& term2_v){

	Unit* unit = units[beta];
	double M = calculateM(units, pair<int,int>(0, units.size()-1), false);
	CenterOfMass r2_0 = calculateMobileBodyCOM(units, pair<int,int>(beta+1, units.size()-1), false);
	double M2 = r2_0.getMass();
	double M1 =  M - M2;

	Point* r_a = unit->r_right;
	Point* e_a = unit->e_right;
	vector<double> r_a_coords = unit->r_right->getCoordinates();
	vector<double> e_a_coords = unit->e_right->getCoordinates();
	double I1[3][3];
	calculateI(I1, units, pair<int,int>(0, beta), INMA);

	// e_a x{(M1/M)*ra+(M2/M)*r20]}
	Point M1Mr_a = Point::multiplyByScalar(*r_a, M1/M);
	Point M2Mri_0 = Point::multiplyByScalar(r2_0, M2/M);
	Point term1 = Point::add(M1Mr_a, M2Mri_0);
	Point e_a_x_term1 = ANMICMath::crossProduct(*e_a, term1);

	//{[M2r2_0 x(eaxra)+I1ea]I^-1}
	Point Mri_0 = Point::multiplyByScalar(r2_0, M2);
	Point ea_ra = ANMICMath::crossProduct(*e_a, *r_a);
	Point Mr2_0_aa_x_r_a = ANMICMath::crossProduct(Mri_0, ea_ra);

	Point I1ea = ANMICMath::multiplyIMatrixByEVector(I1, Utils::vectorToPointer<double>(e_a_coords));
	Point term2 =  Point::add(Mr2_0_aa_x_r_a, I1ea);

	Point term2Iinv = ANMICMath::multitplyPointByMatrix(term2, I_inv);

	term1_v.push_back(e_a_x_term1);
	term2_v.push_back(term2Iinv);
}

