/*
	Utility math library
    Copyright (C) 2015  Markus Kuehbach

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
#pragma once
#ifndef _mymath_h_
#define _mymath_h_

#include <stdlib.h>
#include <math.h>
#include "random.h"
#include <stdint.h>
//#include <vector>

#define QUAT2EUL_ETA 	1e-20
#define RGBRANGE 		255

class randomClass;
typedef double Real;

struct Point{
	double x;
	double y;
	double z;
};


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );
	
	

	//additional functions from vertex model L. Barrales-Mora and helper functions utilized in CORe
	Real areaPolygon( int *poly, int counter );
	int pnpoly(int nvert, Real *vertx, Real *verty, Real testx, Real testy);
	void bubbleSort ( Real arr [ ], int size );
	void sortInt( long arr [], long size );
	void swap ( Real& x, Real& y );
	void swapInt ( long& x, long& y );
	Real convertConc(Real conc,int j);
	
	void factorizeIn3( long n, long * factors );
	char isPrime( long n );
	void counting_sort(int *A, int size, int range);
	void sort(int n, double *ra);
	
	//statistical tests
	void K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob );
	double probks( double alam );
	
	
	//quaternion algebra
	void multiplyQuaternions( double *q, double* p, double* r );
	void multiplyQuaternions2( double *q, double* p, double* r ); //this function is what has been implemented in COReV2
	double distanceBetweenQuaternions( double * q, double * p );
	
	
	//convert orientations among parametrizations, such as Bunge (3,1,3) rotation convention and quaternion
	void euler2quaternion( double * euler, double * q );
	void quaternion2Euler( const double * quat, double * euler );
	
	
	//calculate disorientation among two orientations in various parametrizations
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	void misorientationQuaternionCubic( double* p, double* q, double* quat  );
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );
	
	
	//generate new orientations with some scatter
	void randomOrientation( double * result );
	void randomOriShoemakeEuler( double * result );
	void randomOriShoemakeQuat( double * result );
	//void randomOri ( double * result );
	void randomMisorientation( double theta, double* qr  );
	void randomMisorientationAxisConsidered(  double * qref, double * qr, double maxDev  );
	
	void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
	void newOrientationFromReference( double *oriOri, double deviation, double *newOri );
	
	//identify orientation by RGB scheme
	void devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb); //blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	
	
	//an own Park-Miller PRNG
	randomClass r;
	
	
	//currently not implemented functions
	void devtoaxisEuler2RGB( double *bunge, double *uvw, double maxDev, unsigned char *rgb); //blue channel stretch from 0.0 to maxDev in radian, disorientation to axis of misorientation
	void patalaQuat2RGB ( double *q, unsigned char *rgb); //Patala, Schuh unique coloring of misorientations via HSV2RGB

	
	//###MK:: currently not save to use functions...
	//currently not save to be used, function should bring a defined angular scatter around an axis angle representation of some arbitrary orientation quat
	void newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri);
	
	//currently not save to be used, function converts quaternion representation into Rodrigues-Frank parametrization
	void quaternion2rodrigues ( double * q, double * rodrigues );

	//compatibility with COReV2
	double misorientationCubicOrigInv( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	double misorientationCubicCOReV2( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	double misorientationCubicCorrectCOReV2InvertAndMult( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
		
};

#endif
