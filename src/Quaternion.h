/*
 GraGLeS 2D A grain growth simulation utilizing level set approaches
 Copyright (C) 2015  Christian Miessen, Nikola Velinov
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
