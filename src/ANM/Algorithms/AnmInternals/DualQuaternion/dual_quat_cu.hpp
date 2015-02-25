#ifndef TOOL_BOX_DUAL_QUAT_CU_HPP__
#define TOOL_BOX_DUAL_QUAT_CU_HPP__

#include "quat_cu.hpp"


/** @class Dual_quat_cu
    @brief Representation of a dual quaternion to express rotation and translation

    A dual quaternion (DQ) is based on the algebra of dual numbers. Dual numbers
    are somewhat close to complex numbers (a + ib) as they are writen :
    nd + eps * d where nd is the non-dual part and d the dual part and
    (eps)^2=0.

    Dual quaternion are represented as follow : q0 + eps * qe where q0
    is the non-dual part (a quaternion) and qe the dual part (another quaternion)

    With dual quaternion we can express a rotation and a translation. This
    enables us to substitute rigid transformations matrices to dual quaternions
    and transform a point with the method 'transform()'

    As a dual quaternions is the sum of two quaternions we have to store eight
    coefficients corresponding to the two quaternions.

    To move a point with a rigid transformation (i.e. solely composed
    of a translation and a rotation) you just need to construct the DQ with a
    quaternion wich express the rotation and a translation vector. You can
    now translate and rotate the point at the same time with 'transform()'.

    Linear blending of dual quaternions (DLB) is possible (dq0*w0 + dq1*w1 ...)
    where w0, w1 ... wn are scalar weights whose sum equal one. The weights
    defines the influence of each transformations expressed by the dual
    quaternions dq0, dq1 ... dqn.
    N.B : this is often used to compute mesh deformation for animation systems.

    You can compute DLB with the overloaded operators (+) and (*) and use
    the method transform() of the resulting dual quaternion to deform a point
    according to the DLB.

    @note Article: "Geometric skinning with approximate dual quaternion blending"
 */
class Dual_quat_cu {
    public:

    // -------------------------------------------------------------------------
    /// @name Constructors
    // -------------------------------------------------------------------------

    /// Default constructor generates a dual quaternion with no translation
    /// and no rotation either
    Dual_quat_cu()
    {
        Dual_quat_cu res = dual_quat_from(Quat_cu(), Vec3(0.f, 0.f, 0.f));
        *this = res;
    }


    /// Fill directly the dual quaternion with two quaternion for the non-dual
    /// and dual part
    Dual_quat_cu(const Quat_cu& q0, const Quat_cu& qe)
    {
        _quat_0 = q0;
        _quat_e = qe;
    }

    /// Construct a dual quaternion with a quaternion 'q' which express the
    /// rotation and a translation vector
    Dual_quat_cu(const Quat_cu& q, const Vec3& t)
    {
        Dual_quat_cu res = dual_quat_from(q, t);
        *this = res;
    }

    /// Construct from rigid transformation 't'
    Dual_quat_cu(const Transfo& t)
    {
        Quat_cu q(t);
        Vec3 translation(t.m[3], t.m[7], t.m[11]);
        Dual_quat_cu res = dual_quat_from(q, translation);
        *this = res;
    }


    // -------------------------------------------------------------------------
    /// @name Methods
    // -------------------------------------------------------------------------

    void normalize()
    {
        float norm = _quat_0.norm();
        _quat_0 = _quat_0 / norm;
        _quat_e = _quat_e / norm;
    }

    /// Transformation of point p with the dual quaternion
    Point3 transform(const Point3& p ) const
    {
        // As the dual quaternions may be the results from a
        // linear blending we have to normalize it :
        float norm = _quat_0.norm();
        Quat_cu qblend_0 = _quat_0 / norm;
        Quat_cu qblend_e = _quat_e / norm;

        // Translation from the normalized dual quaternion equals :
        // 2.f * qblend_e * conjugate(qblend_0)
        Vec3 v0 = qblend_0.get_vec_part();
        Vec3 ve = qblend_e.get_vec_part();
        Vec3 trans = (ve*qblend_0.w() - v0*qblend_e.w() + v0.cross(ve)) * 2.f;

        // Rotate
        return qblend_0.rotate(p) + trans;
    }

    /// Rotate a vector with the dual quaternion
    Vec3 rotate(const Vec3& v) const
    {
        Quat_cu tmp = _quat_0;
        tmp.normalize();
        return tmp.rotate(v);
    }

    Dual_quat_cu dual_quat_from(const Quat_cu& q, const Vec3& t) const
    {
        float w = -0.5f*( t.x * q.i() + t.y * q.j() + t.z * q.k());
        float i =  0.5f*( t.x * q.w() + t.y * q.k() - t.z * q.j());
        float j =  0.5f*(-t.x * q.k() + t.y * q.w() + t.z * q.i());
        float k =  0.5f*( t.x * q.j() - t.y * q.i() + t.z * q.w());

        return Dual_quat_cu(q, Quat_cu(w, i, j, k));
    }

    /// Convert the dual quaternion to a homogenous matrix
    /// N.B: Dual quaternion is normalized before conversion
    Transfo to_transformation()
    {
        Vec3 t;
        float norm = _quat_0.norm();

        // Rotation matrix from non-dual quaternion part
        Mat3 m = (_quat_0 / norm).to_matrix3();

        // translation vector from dual quaternion part:
        t.x = 2.f*(-_quat_e.w()*_quat_0.i() + _quat_e.i()*_quat_0.w() - _quat_e.j()*_quat_0.k() + _quat_e.k()*_quat_0.j()) / norm;
        t.y = 2.f*(-_quat_e.w()*_quat_0.j() + _quat_e.i()*_quat_0.k() + _quat_e.j()*_quat_0.w() - _quat_e.k()*_quat_0.i()) / norm;
        t.z = 2.f*(-_quat_e.w()*_quat_0.k() - _quat_e.i()*_quat_0.j() + _quat_e.j()*_quat_0.i() + _quat_e.k()*_quat_0.w()) / norm;

        return Transfo(m, t);
    }

    // -------------------------------------------------------------------------
    /// @name Operators
    // -------------------------------------------------------------------------

    Dual_quat_cu operator+(const Dual_quat_cu& dq) const
    {
        return Dual_quat_cu(_quat_0 + dq._quat_0, _quat_e + dq._quat_e);
    }

    /* Multiplication table:
	*
	* Q1*Q2 | Q2.1 Q2.i Q2.j Q2.k Q2.ei Q2.ej Q2.ek Q2.e
	* ------+---------------------------------------------------
	* Q1.1 | 1 i j k ei ej ek e
	* Q1.i | i -1 k -j -e ek -ej ei
	* Q1.j | j -k -1 i -ek -e ei ej
	* Q1.k | k j -i -1 ej -ei -e ek
	* Q1.ei | ei -e ek -ej 0 0 0 0
	* Q1.ej | ej -ek -e ei 0 0 0 0
	* Q1.ek | ek ej -ei -e 0 0 0 0
	* Q1.e | e ei ej ek 0 0 0 0
	*
	* We can also decomopose the problem into quaternion multiplication:
	*
	* Q = q + \epsilon q0
	* P = p + \epsilon p0
	* Q*P = q*p + \epsilon(q*p0 + q0*p)
	*
	* We can treat quaternion multiplication as:
	*
	* q1 = (r1, v1)
	* q2 = (r2, v2)
	* q1*q2 = (r1*r2-v1.v2, r1*v2 + r2*v1 + v1 x v2)
	*/
    Dual_quat_cu operator*(const Dual_quat_cu& dq) const
	{
    	// From bobbens libdq @github

    	double T[8];
    	double P [] = { 	this->_quat_0.w(),
    						this->_quat_0.i(),
    						this->_quat_0.j(),
    						this->_quat_0.k(),
    						this->_quat_e.i(),
    						this->_quat_e.j(),
    						this->_quat_e.k(),
    						this->_quat_e.w()
    						};
    	double Q [] = {		dq._quat_0.w(),
							dq._quat_0.i(),
							dq._quat_0.j(),
							dq._quat_0.k(),
							dq._quat_e.i(),
							dq._quat_e.j(),
							dq._quat_e.k(),
							dq._quat_e.w()};

    	T[0] = P[0]*Q[0] - P[1]*Q[1] - P[2]*Q[2] - P[3]*Q[3];
    	T[1] = P[0]*Q[1] + P[1]*Q[0] + P[2]*Q[3] - P[3]*Q[2];
    	T[2] = P[0]*Q[2] + P[2]*Q[0] - P[1]*Q[3] + P[3]*Q[1];
    	T[3] = P[0]*Q[3] + P[3]*Q[0] + P[1]*Q[2] - P[2]*Q[1];

		/* Dual unit Quaternion. */
    	T[4] = P[4]*Q[0] + P[0]*Q[4] + P[7]*Q[1] + P[1]*Q[7] - P[6]*Q[2] + P[2]*Q[6] + P[5]*Q[3] - P[3]*Q[5];
    	T[5] = P[5]*Q[0] + P[0]*Q[5] + P[6]*Q[1] - P[1]*Q[6] + P[7]*Q[2] + P[2]*Q[7] - P[4]*Q[3] + P[3]*Q[4];
    	T[6] = P[6]*Q[0] + P[0]*Q[6] - P[5]*Q[1] + P[1]*Q[5] + P[4]*Q[2] - P[2]*Q[4] + P[7]*Q[3] + P[3]*Q[7];
    	T[7] = P[7]*Q[0] + P[0]*Q[7] - P[1]*Q[4] - P[4]*Q[1] - P[2]*Q[5] - P[5]*Q[2] - P[3]*Q[6] - P[6]*Q[3];

    	Quat_cu new_quat_0(T[0],T[1],T[2],T[3]);
    	Quat_cu new_quat_e(T[7],T[4],T[5],T[6]);

		return Dual_quat_cu(new_quat_0, new_quat_e);
	}

    Dual_quat_cu operator*(float scalar) const
    {
        return Dual_quat_cu(_quat_0 * scalar, _quat_e * scalar);
    }

    /// Return a dual quaternion with no translation and no rotation
    static Dual_quat_cu identity()
    {
        return Dual_quat_cu(Quat_cu(1.f, 0.f, 0.f, 0.f),
                            Vec3(0.f, 0.f, 0.f) );
    }

    // -------------------------------------------------------------------------
    /// @name Getters
    // -------------------------------------------------------------------------

    Quat_cu get_dual_part() const { return _quat_e; }

    Quat_cu get_non_dual_part() const { return _quat_0; }

    Quat_cu translation() const { return _quat_e; }

    Quat_cu rotation() const { return _quat_0; }

    void set_rotation( const Quat_cu& q ){ _quat_0 = q; }

    /// Non-dual part of the dual quaternion. It also represent the rotation.
    /// @warning If you want to compute the rotation with this don't forget
    /// to normalize the quaternion as it might be the result of a
    /// dual quaternion linear blending
    /// (when overloaded operators like '+' or '*' are used)
    Quat_cu _quat_0;

    /// Dual part of the dual quaternion which represent the translation.
    /// translation can be extracted by computing
    /// 2.f * _quat_e * conjugate(_quat_0)
    /// @warning don't forget to normalize quat_0 and quat_e :
    /// quat_0 = quat_0 / || quat_0 || and quat_e = quat_e / || quat_0 ||
    Quat_cu _quat_e;
};


#endif // TOOL_BOX_DUAL_QUAT_CU_HPP__
