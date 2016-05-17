/*
 * Quaternion.h
 *
 *  Created on: 17.09.2015
 *      Author: miessen
 */

#include "mymath.h"

#ifndef QUATERNION_H_
#define QUATERNION_H_

class mathMethods;

//Quaternionen class
class Quaternion{
private:
	double q0;
	double q1;
	double q2;
	double q3;
public:
	Quaternion(void);
	Quaternion(double q0, double q1, double q2, double q3);
	Quaternion(double alpha,double x,double y,double z, bool radiant);
	~Quaternion(void);
	Quaternion& operator= ( const Quaternion &rhs);
	//Quaternion& Quaternion::operator+ (Quarternion const& lhs, Quarternion const& rhs);
	//Quaternion& Quaternion::operator- (Quarternion const& lhs, Quarternion const& rhs);
	Quaternion operator* (const Quaternion &rhs);

	bool operator == (Quaternion const& rhs);
	bool operator != (Quaternion const& rhs);

	Quaternion Inverse();
	double getNorm(void);
	void Normalize();
	void Sort(void);
	Quaternion misorientationQuaternionCubic( Quaternion* p );
	double rotationAngleBetweenQuaternions( Quaternion *p );
	void euler2quaternion(double *euler);
	double* quaternion2EulerConst(void) const;
	double* quaternion2Euler(void);
	void randomOriShoemakeQuat(mathMethods* math);
	double misorientationCubicQxQ(Quaternion* p);

	inline double get_q0() {return q0;};
	inline double get_q1() {return q1;};
	inline double get_q2() {return q2;};
	inline double get_q3() {return q3;};
};


#endif /* QUATERNION_H_ */
