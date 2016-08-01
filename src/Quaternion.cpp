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
#include "Quaternion.h"

Quaternion::Quaternion(void) :
	q0(1), q1(0), q2(0), q3(0) {
}

Quaternion::Quaternion(double q0, double q1, double q2, double q3) :
	q0(q0), q1(q1), q2(q2), q3(q3) {
}

Quaternion& Quaternion::operator=(const Quaternion & rhs) {
	if (this != &rhs) //oder if (*this != rhs)
	{
		q0 = rhs.q0;
		q1 = rhs.q1;
		q2 = rhs.q2;
		q3 = rhs.q3; //Copy-Konstruktor
	}
	return *this; //Referenz auf das Objekt selbst zurueckgeben
}

Quaternion::Quaternion(double alpha, double x, double y, double z, bool radiant) {
	double _normVector = 1.0 / sqrt(SQR(x) + SQR(y) + SQR(z));
	if (!radiant) alpha = alpha * PI/180;

	// alpha in radiant
	// Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors,	James Diebel
	q0 = cos(alpha / 2);
	q1 = x * sin(alpha / 2) * _normVector;
	q2 = y * sin(alpha / 2) * _normVector;
	q3 = z * sin(alpha / 2) * _normVector;
}
Quaternion::~Quaternion(void) {
}

double Quaternion::getNorm(void) {
	double qNorm = sqrt(SQR(q1) + SQR(q2) + SQR(q3) + SQR(q0));
	return qNorm;
}
void Quaternion::Normalize() {
	double qNorm = getNorm();
	q0 /= qNorm;
	q1 /= qNorm;
	q2 /= qNorm;
	q3 /= qNorm;
}

Quaternion Quaternion::Inverse() {
	Quaternion r(q0,-q1,-q2,-q3);
	return r;
}

Quaternion Quaternion::operator*(const Quaternion & rhs) {
	//MK::ok, mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	//mind vector cross product notation qp = q0p0 - qbar*pbar + q0 pbar + p0 qbar + qbar cross pbar with bar vector quantities parts of the quaternion
	Quaternion r;
	r.q0 = +q0 * rhs.q0 - q1 * rhs.q1 - q2 * rhs.q2 - q3 * rhs.q3;
	r.q1 = +q1 * rhs.q0 + q0 * rhs.q1 - q3 * rhs.q2 + q2 * rhs.q3;
	r.q2 = +q2 * rhs.q0 + q3 * rhs.q1 + q0 * rhs.q2 - q1 * rhs.q3;
	r.q3 = +q3 * rhs.q0 - q2 * rhs.q1 + q1 * rhs.q2 + q0 * rhs.q3;
	return r;
}


void Quaternion::Sort() {
	double arr[4] = { q0, q1, q2, q3 };
	int last = 4 - 2;
	int isChanged = 1;

	while (last >= 0 && isChanged) {
		isChanged = 0;
		for (int k = 0; k <= last; k++)
			if (arr[k] > arr[k + 1]) {
				SWAP(arr[k], arr[k + 1]);
				isChanged = 1;
			}
		last--;
	}
	q3 = arr[0];
	q2 = arr[1];
	q1 = arr[2];
	q0 = arr[3];
}

Quaternion Quaternion::misorientationQuaternionCubic(Quaternion* p) {
	//MK::ok
	Quaternion qm1; //Inverse of quaternion q

	//Inverse of quaternion q is the same like the conjugate for unit quaternions
	qm1 = (*this).Inverse();

	Quaternion result; //Resulting misorientation quaternion, m = pq-1

	result = qm1*(*p); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3


	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688
	double r0[6][4]; //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = result.q0;
	double b = result.q1;
	double c = result.q2;
	double d = result.q3;

	double fac = 1.0 / sqrt(2.0); //0.70710678;

	//six fundamental quaternions
	r0[0][0] = (a + b) * fac;
	r0[0][1] = (a - b) * fac;
	r0[0][2] = (c + d) * fac;
	r0[0][3] = (c - d) * fac;
	r0[1][0] = (a + c) * fac;
	r0[1][1] = (a - c) * fac;
	r0[1][2] = (b + d) * fac;
	r0[1][3] = (b - d) * fac;
	r0[2][0] = (a + d) * fac;
	r0[2][1] = (a - d) * fac;
	r0[2][2] = (b + c) * fac;
	r0[2][3] = (b - c) * fac;
	r0[3][0] = (a + b + c + d) * 0.5;
	r0[3][1] = (a + b - c - d) * 0.5;
	r0[3][2] = (a - b + c - d) * 0.5;
	r0[3][3] = (a - b - c + d) * 0.5;
	r0[4][0] = (a + b + c - d) * 0.5;
	r0[4][1] = (a + b - c + d) * 0.5;
	r0[4][2] = (a - b + c + d) * 0.5;
	r0[4][3] = (a - b - c - d) * 0.5;
	r0[5][0] = a;
	r0[5][1] = b;
	r0[5][2] = c;
	r0[5][3] = d;


	int mi=0;
	double max = 0.0;

	for (int i = 0; i < 6; i++) //Determing the quaternion with the maximal component and the component itself
		for (int j = 0; j < 4; j++) {
			if (fabs(r0[i][j]) > max) {
				max = fabs(r0[i][j]);
				mi = i;
			}
		}

	result.q0 = fabs(r0[mi][0]); //Disorientation requires all components positive
	result.q1 = fabs(r0[mi][1]);
	result.q2 = fabs(r0[mi][2]);
	result.q3 = fabs(r0[mi][3]);

	result.Sort(); //Sorting into ascending order, because a desorientation in the SST
	//requires a quaternion with q0>=q1>=q2>=q3 which represents a minimal
	return result;
	//additionally it is required that rq3 >= sum of all others and rq3 * (sqrt(2)-1) >= rq2
}

void Quaternion::euler2quaternion(double *euler) {
	/*20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006
	 Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	 //cos(a+b) = c(a+b) = cacb-sasb
	 //cos(a-b) = c(a-b) = cacb+sasb
	 //sin(a+b) = s(a+b) = sacb+casb
	 //sin(a-b) = s(a-b) = sacb-casb
	 */

	double p1 = euler[0];
	double t = euler[1];
	double p2 = euler[2];

	Real co1 = cos(t / 2);
	Real s1 = sin(t / 2);

	//double test[4]={0};
	//quaternion2Euler( p,test);

	q0 = co1 * cos((p1 + p2) / 2);
	q1 = s1 * cos((p1 - p2) / 2);
	q2 = s1 * sin((p1 - p2) / 2);
	q3 = co1 * sin((p1 + p2) / 2);
}

double* Quaternion::quaternion2Euler(void) {
	//convention: Bunge, ZXZ, equal to case (3,1,3) as
	//  analyzed in Diebel, James, 2006:
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401  413
	//dont forget to define QUAT2EUL_ETA 1e-20

	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);

	Real y0 = 2 * q1 * q3 - 2 * q0 * q2;
	Real x0 = 2 * q2 * q3 + 2 * q0 * q1;
	Real y1 = 2 * q1 * q3 + 2 * q0 * q2;
	Real x1 = -2 * q2 * q3 + 2 * q0 * q1;

	if (cosPHI > 1.)
		cosPHI = 1.;

	if (SQR(1. - cosPHI) <= QUAT2EUL_ETA)
		PHI = 0.;
	else
		PHI = acos(cosPHI);

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a Quaternion does not define a Bunge Euler angle uniquely

	if (sP != 0) {
		phi2 = atan2(y0 / sP, x0 / sP);
		phi1 = atan2(y1 / sP, x1 / sP);
	} else {
		phi1 = atan2((2 * q1 * q2 + 2 * q0 * q3),
				SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3));
		phi2 = 0.;
	}

	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * PI;
	if (phi2 < 0.0)
		phi2 += 2 * PI;

	double *euler = new double[3];
	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
	return euler;
}

double* Quaternion::quaternion2EulerConst(void) const {
	//convention: Bunge, ZXZ, equal to case (3,1,3) as
	//  analyzed in Diebel, James, 2006:
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401  413
	//dont forget to define QUAT2EUL_ETA 1e-20

	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3) - SQR(q2) - SQR(q1) + SQR(q0);

	Real y0 = 2 * q1 * q3 - 2 * q0 * q2;
	Real x0 = 2 * q2 * q3 + 2 * q0 * q1;
	Real y1 = 2 * q1 * q3 + 2 * q0 * q2;
	Real x1 = -2 * q2 * q3 + 2 * q0 * q1;

	if (cosPHI > 1.)
		cosPHI = 1.;

	if (SQR(1. - cosPHI) <= QUAT2EUL_ETA)
		PHI = 0.;
	else
		PHI = acos(cosPHI);

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a quarternion does not define a Bunge Euler angle uniquely

	if (sP != 0) {
		phi2 = atan2(y0 / sP, x0 / sP);
		phi1 = atan2(y1 / sP, x1 / sP);
	} else {
		phi1 = atan2((2 * q1 * q2 + 2 * q0 * q3),
				SQR(q0) + SQR(q1) - SQR(q2) - SQR(q3));
		phi2 = 0.;
	}

	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * PI;
	if (phi2 < 0.0)
		phi2 += 2 * PI;

	double *euler = new double[3];
	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;
	return euler;
}

void Quaternion::randomOriShoemakeQuat(mathMethods* math) {
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//##MK::mind order
	double s = math->r.parkMiller();
	double sigma1 = sqrt(1 - s);
	double sigma2 = sqrt(s);
	double theta1 = 2 * _PI_ * math->r.parkMiller();
	double theta2 = 2 * _PI_ * math->r.parkMiller();

	q0 = fabs(sigma1 * sin(theta1));
	q1 = fabs(sigma1 * cos(theta1));
	q2 = fabs(sigma2 * sin(theta2));
	q3 = fabs(sigma2 * cos(theta2));

	//normalize
	Normalize();
}

double Quaternion::misorientationCubicQxQ(Quaternion* p) {
	int i;
	Quaternion qm1 = (*this).Inverse(); //Inverse of quaternion this
	Quaternion r = qm1*(*p); //Resulting quaternion, rotation of the two previous quaternions pq-1

	//Now, we have to determine the smallest angle.

	Real r0[6][4]; //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r.q0;
	double b = r.q1;
	double c = r.q2;
	double d = r.q3;

	Real fac = 1.0 / sqrt(2.0); //0.70710678;

	//six fundamental quaternions
	r0[0][0] = (a + b) * fac;
	r0[0][1] = (a - b) * fac;
	r0[0][2] = (c + d) * fac;
	r0[0][3] = (c - d) * fac;
	r0[1][0] = (a + c) * fac;
	r0[1][1] = (a - c) * fac;
	r0[1][2] = (b + d) * fac;
	r0[1][3] = (b - d) * fac;
	r0[2][0] = (a + d) * fac;
	r0[2][1] = (a - d) * fac;
	r0[2][2] = (b + c) * fac;
	r0[2][3] = (b - c) * fac;
	r0[3][0] = (a + b + c + d) * 0.5;
	r0[3][1] = (a + b - c - d) * 0.5;
	r0[3][2] = (a - b + c - d) * 0.5;
	r0[3][3] = (a - b - c + d) * 0.5;
	r0[4][0] = (a + b + c - d) * 0.5;
	r0[4][1] = (a + b - c + d) * 0.5;
	r0[4][2] = (a - b + c + d) * 0.5;
	r0[4][3] = (a - b - c - d) * 0.5;
	r0[5][0] = a;
	r0[5][1] = b;
	r0[5][2] = c;
	r0[5][3] = d;

	Real omega = 0.0;

	for (i = 0; i < 6; i++)
		for (int j = 0; j < 4; j++)
			if (fabs(r0[i][j]) > omega)
				omega = fabs(r0[i][j]);

	if (omega > 1.0) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega = 2 * acos(omega);
	return omega;
}

