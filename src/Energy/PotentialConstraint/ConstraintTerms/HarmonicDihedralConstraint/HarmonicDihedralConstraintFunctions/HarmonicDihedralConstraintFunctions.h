/////////////////////////////////////////////////////////////////////////////
/// \file HarmonicDihedralConstraintFunctions.h
///
/// \brief It contains harmonic dihedral constraint functions
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author arincon
/// \author vgil
/// \date 19/09/2012
/////////////////////////////////////////////////////////////////////////////

#include <vector>

/////////////////////////////////////////////////////////////////////////////
/// \brief It contains the harmonic dihedral constraint functions derived (and generated) with wxMaxima
/// scripts
///
/// \author arincon
/// \author vgil
/// \date 19/09/2012
/////////////////////////////////////////////////////////////////////////////
namespace HarmonicDihedralConstraintFunctions {
	double put_in_pi_minus_pi_range(double angle);

	double to_0_2PI_range(double angle);

	double rad_subtraction(double ang1, double ang2);

	double calculateDihedralAngleWithArcTanFunction(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd);

	double calculateEnergyWithArcTanFunction(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);

	std::vector<double> calculateGradient(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);

	std::vector<double> calculateHessian(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);

	// Partial derivatives for gradient calculation (autogenerated with maxima)
	double derivativeXa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeYa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeZa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeXb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double derivativeZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);

	// Partial second derivatives for hessian calculation (autogenerated with maxima)
	double hessianXaXa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaYa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaZa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaXb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXaZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaYa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaZa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaXb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYaZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaZa(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaXb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZaZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbXb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXbZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbYb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYbZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbZb(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZbZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcXc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXcZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYcYc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYcZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYcXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYcYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYcZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZcZc(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZcXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZcYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZcZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXdXd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXdYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianXdZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYdYd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianYdZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
	double hessianZdZd(double k, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, double xd, double yd, double zd, double eq);
}
