/////////////////////////////////////////////////////////////////////////////
///
/// \author vgil
/// \date 04/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef ANMICMOVER_H_
#define ANMICMOVER_H_

#include <vector>
#include "../DualQuaternion/quat_cu.hpp"

class Unit;

namespace ANMICMovement {

	Quat_cu calculate_quaternion_for_unit(double , std::vector<Unit*>& , int );
	void apply_rotations_to_molecule_units(std::vector<Unit*>& , std::vector<double>& , double );

}

#endif /* MOVER_H_ */
