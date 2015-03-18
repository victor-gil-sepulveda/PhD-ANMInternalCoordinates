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

	Point pivot(*(units[unit_number]->r_right));
	Point Q(units[unit_number]->right_torsion_bond_atom->toPoint());

	Point axis = Point::subtract(Q, pivot);
	return Quat_cu(cos(torsional_disp/2.),
					axis.getX()*sin(torsional_disp/2.),
					axis.getY()*sin(torsional_disp/2.),
					axis.getZ()*sin(torsional_disp/2.)
	);
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
void ANMICMovement::apply_rotations_to_molecule_units(vector<Unit*>& units,
		vector<double>& angular_increments,
		double increment_divider){
	Quat_cu q_i;
	Quat_cu s_im1;
	Point3 T;
	Point3 T_im1;

	for(unsigned int i = 0; i < units.size()-1; ++i){
		// Accumulate rotations and translations
		double torsion_increment = angular_increments[i] * increment_divider;
		q_i = ANMICMovement::calculate_quaternion_for_unit(torsion_increment, units, i);

		Point Q_p(units[i]->right_torsion_bond_atom->toPoint());
		Point3 Q_i(Q_p.getX(), Q_p.getY(), Q_p.getZ());

		if (i == 0){
			Point3 q_iQ_iq_i_conj = q_i.rotate(Q_i);
			T = Q_i +( - q_iQ_iq_i_conj);

			// Result accumulation for next iter
			T_im1 = T;
			s_im1 = q_i;
		}
		else{
			Point3 T_im1mQ_i = T_im1 + (- Q_i);
			Point3 q_iT_im1mQ_i_q_i_conj = q_i.rotate(T_im1mQ_i);
			T = q_iT_im1mQ_i_q_i_conj + Q_i;

			// Result accumulation for next iter
			T_im1 = T;
			s_im1 = q_i*s_im1;
		}

		// Rotate unit i+1
		Dual_quat_cu M(s_im1, T);
		vector<Atom*> atoms = units[i+1]->getAllAtoms();
		for (unsigned int j = 0; j < atoms.size(); ++j){
			Atom* atom = atoms[j];
			Point3 atom_v3(	atom->getX(), atom->getY(), atom->getZ());
			Point3 new_atom_p3 = M.transform(atom_v3);
			atom->setX(new_atom_p3.x);
			atom->setY(new_atom_p3.y);
			atom->setZ(new_atom_p3.z);
		}
	}
}
