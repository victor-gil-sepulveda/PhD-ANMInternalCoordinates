/////////////////////////////////////////////////////////////////////////////
///
/// \author vgil
/// \date 04/02/2015
/////////////////////////////////////////////////////////////////////////////

#include "ANMICMovement.h"
#include "../../../../Tools/Math/Point.h"
#include "../../../../Molecules/Atom.h"
#include "../../../ModesCalculator/Internals/CoarseGrainModel/Unit.h"
#include "../DualQuaternion/dual_quat_cu.hpp"
#include "../DualQuaternion/point3.hpp"
#include "../DualQuaternion/quat_cu.hpp"
#include "../../../ModesCalculator/Internals/CoarseGrainModel/UnitTools.h"
#include "../../../ModesCalculator/Internals/MatrixCalculationFunctions/ANMICMath.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// Calculates the quaternion representing a rotation in the axis formed from unit i last atom to unit i+1
/// first atom.
/// TODO: Moving it to UnitTools
///
/// \param torsional_disp [In] The actual angular torsion (in radians)
/// \param units [In] A vector holding all units.
/// \param unit_number [In] The unit we are working with inside the units array.
///
/// \return A quaternion.
///
/// \author vgil
/// \date 31/10/2014
///////////////////////////////////////////////////////////////
Quat_cu ANMICMovement::calculate_quaternion_for_unit(double torsional_disp,
		vector<Unit*>& units,
		int unit_number){

	Point P(units[unit_number]->left_torsion_bond_atom->toPoint());
	Point Q(units[unit_number]->right_torsion_bond_atom->toPoint());

	Point axis = Point::subtract(Q, P);
	axis.normalize();

	double sin_w = sin(torsional_disp/2.);
	Quat_cu q(cos(torsional_disp/2.),
				axis.getX()*sin_w,
				axis.getY()*sin_w,
				axis.getZ()*sin_w);
	q.normalize();
	return q;
}

void calculate_axes(vector<Unit*>& units, vector<Point>& axes){
	axes.clear();
	for(unsigned int i = 0; i < units.size()-1; ++i){
		Point P(units[i]->left_torsion_bond_atom->toPoint());
		Point Q(units[i]->right_torsion_bond_atom->toPoint());

		Point axis = Point::subtract(Q, P);
		axis.normalize();

		axes.push_back(axis);
	}
}

void calculate_Qs(vector<Unit*>& units, vector<Point>& Qs){
	Qs.clear();
	for(unsigned int i = 0; i < units.size()-1; ++i){
		Qs.push_back(units[i]->right_torsion_bond_atom->toPoint());
	}
}

Quat_cu calculate_quaternion(double torsional_disp, Point axis){

	double sin_w = sin(torsional_disp/2.);
	Quat_cu q(cos(torsional_disp/2.),
				axis.getX()*sin_w,
				axis.getY()*sin_w,
				axis.getZ()*sin_w);
	q.normalize();
	return q;
}

void ANMICMovement::apply_rotations_to_molecule_units(vector<Unit*>& units,
		vector<double>& angular_increments,
		double increment_divider){

	Quat_cu q_i;

	vector<vector<double> > M_prev(4, vector<double>(4, 0.0));
	for(unsigned int i = 0; i < 4; ++i)
		M_prev[i][i] = 1.0;

	vector<Point> axes, Q;
	calculate_Qs(units, Q);
	calculate_axes(units, axes);

	// Do transformations
	for(unsigned int i = 0; i < units.size()-1; ++i){
		Point axis = Point::subtract(units[i]->right_torsion_bond_atom->toPoint(),
				units[i]->left_torsion_bond_atom->toPoint());
		axis.normalize();
		Point Qi = units[i]->right_torsion_bond_atom->toPoint();

		double w = angular_increments[i] * increment_divider;

		double x = axis.getX();
		double y = axis.getY();
		double z = axis.getZ();

		double s = sin(w);
		double c = cos(w);
		double t = 1 - c;

		vector<vector<double> > Ri(3, vector<double>(3,0));
		// MAL EN CHOI!? -> http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/
		Ri[0][0] = t*x*x + c; 	  Ri[0][1] = t*x*y - z*s; Ri[0][2] = t*x*z + y*s;
		Ri[1][0] = t*x*y + z*s;   Ri[1][1] = t*y*y + c;   Ri[1][2] = t*y*z - x*s;
		Ri[2][0] = t*x*z - y*s;   Ri[2][1] = t*y*z + x*s; Ri[2][2] = t*z*z + c;

		// Get Qi
		vector<vector<double> > mQi(3,vector<double>(1,0));
		mQi[0][0] = Qi.getX(); mQi[1][0] = Qi.getY(); mQi[2][0] = Qi.getZ();

		// Perform the calculations
		// Ti -> Ri(Tprev) + Qi - Ri(Qi)
		vector<vector<double> > RiQi,Ti(3, vector<double>(3,0));
		// Calc. Ri(Qi)
		ANMICMath::multiplyMatrixByMatrix(Ri, mQi, RiQi);
		// Calc. Ri(Tprev) + Qi - Ri(Qi)
		Ti[0][0] =  mQi[0][0] - RiQi[0][0];
		Ti[1][0] =  mQi[1][0] - RiQi[1][0];
		Ti[2][0] =  mQi[2][0] - RiQi[2][0];

		// Mount transf. matrix for this bond rotation
		vector<vector<double> > Mi(4, vector<double>(4,0));
		Mi[0][0] = Ri[0][0]; Mi[0][1] = Ri[0][1]; Mi[0][2] = Ri[0][2]; Mi[0][3] = Ti[0][0];
		Mi[1][0] = Ri[1][0]; Mi[1][1] = Ri[1][1]; Mi[1][2] = Ri[1][2]; Mi[1][3] = Ti[1][0];
		Mi[2][0] = Ri[2][0]; Mi[2][1] = Ri[2][1]; Mi[2][2] = Ri[2][2]; Mi[2][3] = Ti[2][0];
		Mi[3][0] = 0.0; 	 Mi[3][1] = 0.0;   	  Mi[3][2] = 0.0;      Mi[3][3] = 1.0;

//		vector<vector<double> > tmp;
//		ANMICMath::multiplyMatrixByMatrix(Mi, M_prev, tmp);
//
//		M_prev[0][0] = tmp[0][0]; M_prev[0][1] = tmp[0][1]; M_prev[0][2] = tmp[0][2]; M_prev[0][3] = tmp[0][3];
//		M_prev[1][0] = tmp[1][0]; M_prev[1][1] = tmp[1][1]; M_prev[1][2] = tmp[1][2]; M_prev[1][3] = tmp[1][3];
//		M_prev[2][0] = tmp[2][0]; M_prev[2][1] = tmp[2][1]; M_prev[2][2] = tmp[2][2]; M_prev[2][3] = tmp[2][3];
//		M_prev[3][0] = tmp[3][0]; M_prev[3][1] = tmp[3][1]; M_prev[3][2] = tmp[3][2]; M_prev[3][3] = tmp[3][3];

		// Apply
		vector<Atom*> atoms;
		UnitTools::getAllAtomsFromUnitRange(units, atoms, i+1,units.size()-1, false);
//		vector<Atom*> atoms = units[i+1]->getAllAtoms();
		for (unsigned int j = 0; j < atoms.size(); ++j){
			Atom* atom = atoms[j];
			vector<vector<double> > p(4,vector<double>(1,0));
			p[0][0] = atom->getX(); p[1][0] = atom->getY(); p[2][0] = atom->getZ(); p[3][0] = 1.0;
			vector<vector<double> > result;
			ANMICMath::multiplyMatrixByMatrix(Mi, p, result);

			atom->setX(result[0][0]); atom->setY(result[1][0]); atom->setZ(result[2][0]);
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Moves the torsion angles of a molecule (delimited by its units) using the algorithm
/// described in "On Updating Torsion Angles of Molecular Conformations, Choi, 2006"
/// Variables follow the nomenclature used in the paper.
/// Note that this function will apply the increments PARTIALLY.
///
/// \param units [In] A vector holding all molecule units.
/// \param angular_increments [In] The angular increment for each of the dihedrals of the molecule.
/// \param increment_divider [In] A multiplicative factor for all angular increments.
///
/// \author vgil
/// \date 31/10/2014
///////////////////////////////////////////////////////////////
//void ANMICMovement::apply_rotations_to_molecule_units(vector<Unit*>& units,
//		vector<double>& angular_increments,
//		double increment_divider){
//
//	Quat_cu q_i;
//	Quat_cu s_im1;
//	Point3 T;
//	Point3 T_im1;
//	vector<Point> axes, Q;
//	calculate_axes(units, axes);
//	calculate_Qs(units, Q);
//
//	for(unsigned int i = 0; i < units.size()-1; ++i){
//		// Accumulate rotations and translations
//		double torsion_increment = angular_increments[i] * increment_divider;
//		cout<<"DBG: applying increment "<<torsion_increment<<endl;
//		q_i = calculate_quaternion(torsion_increment, axes[i]);
//		Point3 Q_i(Q[i].getX(), Q[i].getY(), Q[i].getZ());
//
//		if (i == 0){
//			Point3 q_iQ_iq_i_conj = q_i.rotate(Q_i);
//			T = Q_i +( - q_iQ_iq_i_conj);
//			cout<<T<<endl;
//
//			// Result accumulation for next iter
//			T_im1 = T;
//			s_im1 = q_i;
//			s_im1.normalize();
//		}
//		else{
//			Point3 T_im1mQ_i = T_im1 + (- Q_i);
//			Point3 q_iT_im1mQ_i_q_i_conj = q_i.rotate(T_im1mQ_i);
//			T = q_iT_im1mQ_i_q_i_conj + Q_i;
//
//			// Result accumulation for next iter
//			T_im1 = T;
//			s_im1 = q_i*s_im1;
//			s_im1.normalize();
//		}
//
//		// Rotate unit i+1
//		vector<Atom*> atoms = units[i+1]->getAllAtoms();
////		vector<Atom*> atoms;
////		UnitTools::getAllAtomsFromUnitRange(units, atoms, i+1,units.size()-1, false);
//
//		// Apply with M multiplication
//		Mat3 _M = s_im1.to_matrix3();
//		vector<vector<double> > M(4, vector<double>(4,0));
//		M[0][0] = _M.a; M[0][1] = _M.b; M[0][2] = _M.c; M[0][3] = T.x;
//		M[1][0] = _M.d; M[1][1] = _M.e; M[1][2] = _M.f; M[1][3] = T.y;
//		M[2][0] = _M.g; M[2][1] = _M.h; M[2][2] = _M.i; M[2][3] = T.z;
//		M[3][0] = 0.0;  M[3][1] = 0.0;  M[3][2] = 0.0;  M[3][3] = 1.0;
//		for (unsigned int j = 0; j < atoms.size(); ++j){
//			Atom* atom = atoms[j];
//			vector<vector<double> > p(4,vector<double>(1,0));
//			p[0][0] = atom->getX(); p[1][0] = atom->getY(); p[2][0] = atom->getZ(); p[3][0] = 1.0;
//
//			vector<vector<double> > result;
//
//			ANMICMath::multiplyMatrixByMatrix(M, p, result);
////			cout<<result[0][0]<<" "<<result[1][0]<<" "<<result[2][0]<<endl;
////			cout<<result.size()<<" ... "<<result[0].size()<<endl;
//
//			atom->setX(result[0][0]);
//			atom->setY(result[1][0]);
//			atom->setZ(result[2][0]);
//		}
//
//	}
//}

