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
#include "mymath.h"
#include "applic.h"
#include "random.h"


//###LB: Unless indicated, all functions in radians.
//###MK: Quaternion algebra: q = -q define an equivalent rotation.

mathMethods::mathMethods( void )
{
	r.init( -4256 );
}

mathMethods::~mathMethods( void )
{
}

char mathMethods::isPrime( long n )
{
	long n_root = (long) ( sqrt((double)n) + 0.5 );

	for( long i = 2; i <= n_root; i++ )
		if( n % i == 0 ) return 0;

	return 1;
}

void mathMethods::factorizeIn3( long n, long * factors )
{
	long f[3]={0};

	if( isPrime( n ) )
	{
		factors[0]=n;
		factors[1]=1;
		factors[2]=1;
		return;
	}

	long nhalf = (long) (0.5*n);
	long *Fac = (long *) calloc( nhalf, sizeof( long ) );

	int count=0;

	for( long i = nhalf; i > 0; i-- )
	{
		if( n % i==0 )
		{
			Fac[count]=i;
			count++;
		}
	}

	long ft = count;
	count=0;
	long sum = 3*n;

	for( long i=0;i<ft;i++ )
		for( long j=0;j<ft;j++ )
		{
			long p = Fac[i]*Fac[i]*Fac[j];
			long s = Fac[i]+Fac[i]+Fac[j];

			if( p == n && s < sum )
			{
				f[0]=Fac[i];
				f[1]=Fac[i];
				f[2]=Fac[j];
				sum = s;
			}
		}

		for( long i=0;i<ft;i++ )
			for( long j=i+1;j<ft;j++ )
				for( long k=j+1;k<ft;k++ )
				{
					long p = Fac[i]*Fac[j]*Fac[k];
					long s = Fac[i]+Fac[j]+Fac[k];
					if( p == n && s < sum )
					{
						f[0]=Fac[i];
						f[1]=Fac[j];
						f[2]=Fac[k];
						sum = s;
					}
				}

		sortInt( f,3 );
		factors[0]=f[2];
		factors[1]=f[1];
		factors[2]=f[0];
		free(Fac);
}

double mathMethods::areaPolygon( int *poly, int counter )
{
        double sum=0;

        for( int i=0;i<(2*counter-2);i+=2 )
        {
                sum += (poly[i]+poly[i+2])*(poly[i+3]-poly[i+1]);
        }
        if( sum < 0 ) sum *= -1;
        return 0.5 * sum;
}

int mathMethods::pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
        int i, j, c = 0;
        for (i = 0, j = nvert-1; i < nvert; j = i++)
        {
                if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
                        c = !c;
        }
        return c;
}

void mathMethods::bubbleSort ( Real arr [ ], int size ) // Sort components of a quaternion in ASECENDING order
 { 
    int last = size - 2; 
    int isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    { 
            isChanged = 0; 
            for ( int k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swap ( arr[k], arr[k+1] ); 
                    isChanged = 1; 
                } 
            last--; 
    } 
 }

void mathMethods::sortInt( long arr [ ], long size ) // Sort integers
 { 
    long last = size - 2; 
    long isChanged = 1; 

    while ( last >= 0 && isChanged ) 
    { 
            isChanged = 0; 
            for ( long k = 0; k <= last; k++ ) 
                if ( arr[k] > arr[k+1] ) 
                { 
                    swapInt ( arr[k], arr[k+1] ); 
                    isChanged = 1; 
                } 
            last--; 
    } 
 }

void mathMethods::swapInt( long& x, long& y ) //Required for the sorting
{
    long temp;
    temp = x;
    x = y;
    y = temp;
}
 
void mathMethods::swap( Real& x, Real& y ) //Required for the sorting
{
    Real temp;
    temp = x;
    x = y;
    y = temp;
}

void mathMethods::counting_sort(int *A, int size, int range) 
{
	int *B = new int[size];
	int *C = new int[range + 1];

	for(int i = 0; i <= range; ++i)
		C[i] = 0;

	for(int i = 0; i < size; ++i)
		++C[A[i]];

	for(int i = 1; i <= range; ++i)
		C[i] += C[i - 1];

	for (int i = size - 1; i >= 0; --i)
	{
		B[C[A[i]] - 1] = A[i];
		--C[A[i]];
	}

	for (int i = 0; i < size; ++i)
		A[i] = B[i];

  delete[] B;
  delete[] C;
}

////////////////////////////- Convert Weightpercent to Atompercent/100 -////////////////////////////
Real mathMethods::convertConc(Real conc,int elemTyp)
{
/*     int k,kMax;        //Al elemTyp= Cr-1 Cu-2 Fe-3 Mg-4 Mn-5 Si-6 Ti-7 Zn-8 (elemTyp=8)
     kMax = 9;         //number of chemical elements considered in ClaNG
     Real ElemC[ORICLASSMAX];
     Real SumC=0.0;
     for (k = 0 ;k<kMax ;k++ )
     {
            if(k == (elemTyp-1)) ElemC[k] = conc/100.0;
            else ElemC[k]=parm.ElemC[k]/100.0;
            if(k == 0) ElemC[k]= 1 - 1/100.0 * (parm.ElemC[k+1]+parm.ElemC[k+2]+parm.ElemC[k+3]+parm.ElemC[k+4]+parm.ElemC[k+5]+parm.ElemC[k+6]+parm.ElemC[k+7]+parm.ElemC[k+8]);
//            printf("%f\n", ElemC[k]);
            SumC+= ElemC[k]/parm.ElemMolM[k];
     }
     QUICKASSERT(SumC > 0.0);
     conc = ElemC[elemTyp-1]/(parm.ElemMolM[elemTyp-1]*SumC);

     return conc;*/
	return 0;
}

void mathMethods::K_S_Test( double * data1, unsigned long n1, double * data2, unsigned long n2, double * d, double * prob )
{
	unsigned long j1 = 1;
	unsigned long j2 = 1;
	double d1, d2, dt, en1, en2, en;
	double fn1 = 0.0;
	double fn2 = 0.0;

	sort( n1, data1 );
	sort( n2, data2 );

	en1 = n1;
	en2 = n2;
	*d = 0.0;

	while( j1 <= n1 && j2 <= n2 ) {
		if( (d1=data1[j1]) <= (d2=data2[j2]) ) fn1 = j1++ / en1;
		if( d2 <= d1 ) fn2 = j2++/en2;
		if( (dt=fabs(fn2-fn1)) > *d ) *d = dt;
	}
	en = sqrt(en1*en2/(en1+en2));
	*prob = probks( (en+0.12+0.11/en)*(*d) );
}

double mathMethods::probks( double alam )
{
	int j;
	double a2, term;
	double fac = 2.0;
	double sum = 0.0;
	double termbf = 0.0;

	a2 = -2.0 * SQR(alam);
	for( j=1; j<=100; j++ )
	{
		term = fac * exp( a2 * SQR( j ) );
		sum += term;
		if( fabs(term) <= EPS1 * termbf || fabs(term) <= EPS2 * sum ) return sum;
		fac = -fac;
		termbf = fabs( fac );
	}
	return 1.0;
}

void mathMethods::sort(int n, double *ra)
{
	int l, j, ir, i;
	float rra;

	l = (n >> 1) + 1;
	ir = n;

	for ( ; ; )
	{
		if (l > 1)					/* still in hiring phase */
			rra = ra[--l];
		else						/* in retirement-and-promotion phase */
		{
			rra = ra[ir];           /* clear a space at end of array */
			ra[ir]=ra[1];			/* retire the top of the heap into it */
			if (--ir == 1) 			/* done with last promotion */
			{
				ra[1] = rra;
				return;
			}						/* end if */
		}							/* end else */
		i = l;						/* whether we are in the hiring phase */
		j = l << 1;					/* or promotion phase, we here set up */
		while ( j <= ir )
		{
			if ( (j < ir) && (ra[j] < ra[j + 1]) )
				++j;				/* compare to the better underling */
				if ( rra < ra[j] )	/* demote rra */
				{
					ra[i] = ra[j];
					j += (i = j);
				}
				else
					j = ir + 1;		/* this is rra's level; set j to */
		}                           /* terminate the sift-down */
		ra[i] = rra;				/* put rra into its slot */
	}
}

double mathMethods::misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 )
{
	int i;

	Real p[4] = {q01,q11,q21,q31};
	Real q[4] = {q02,q12,q22,q32};

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++) {
		qm1[i]=q[i];
		if( i>0 ) qm1[i] *= -1;
	}
	//qm1[0] *= -1;

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions(qm1, p, r );

	//Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

    //Note: The notation r0 is due to the definition of the quaternion which lie
    //in the fundamental zone, this vector possesses the smallest angle, in such a way
    //that r0 is actually the scalar part of this quaternion

	double a = r[0]; 		
	double b = r[1]; 		
	double c = r[2]; 		
	double d = r[3];

	Real fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	Real omega = 0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 ) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}

void mathMethods::randomMisorientationAxisConsidered(  double * qref, double * qr, double maxDev  )
{
        Real theta = cos( 0.5 * _PI_ );

        double q[4] = {0};

        double dev = cos(0.5 * maxDev);

        double refNorm = sqrt( SQR(qref[0]) + SQR(qref[1]) + SQR(qref[2]) + SQR(qref[3]) );

        qref[0] /= refNorm;
        qref[1] /= refNorm;
        qref[2] /= refNorm;
        qref[3] /= refNorm;

        while( theta < dev  )
        {
                double s = r.parkMiller();
                double sigma1 = sqrt(1-s);
                double sigma2 = sqrt(s);
                double theta1 = 2 * _PI_ * r.parkMiller();
                double theta2 = 2 * _PI_ * r.parkMiller();

                q[0]=fabs(sigma2*cos(theta2));
                q[1]=fabs(sigma1*sin(theta1));
                q[2]=fabs(sigma1*cos(theta1));
                q[3]=fabs(sigma2*sin(theta2));

                double norm = sqrt(SQR(q[0])+SQR(q[1])+SQR(q[2])+SQR(q[3]));

                q[0] /= norm;
                q[1] /= norm;
                q[2] /= norm;
                q[3] /= norm;

                bubbleSort( q, 4 );

                theta = q[3]*qref[0] + q[2]*qref[1] + q[1]*qref[2] + q[0]*qref[3];
        }
        qr[0]=q[3];
        qr[1]=q[2];
        qr[2]=q[1];
        qr[3]=q[0];

}

double mathMethods::distanceBetweenQuaternions( double * q, double * p )
{
        double _qNorm = 1 / sqrt(SQR(q[1])+SQR(q[2])+SQR(q[3])+SQR(q[0]));
        double _pNorm = 1 / sqrt(SQR(p[1])+SQR(p[2])+SQR(p[3])+SQR(p[0]));
        return (2 * acos( _qNorm * _pNorm * ( q[0]*p[0] + q[1]*p[1] + q[2]*p[2] + q[3]*p[3] ) ) );
}

void mathMethods::randomOrientation( double * result )
{
	double q[4]={0,0,0,0};

	double s = r.parkMiller();
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2 * _PI_ * r.parkMiller();
	double theta2 = 2 * _PI_ * r.parkMiller();

	q[0]=fabs(sigma2*cos(theta2));
	q[1]=fabs(sigma1*sin(theta1));
	q[2]=fabs(sigma1*cos(theta1));
	q[3]=fabs(sigma2*sin(theta2));

	bubbleSort( q,4 );

	Real qr[4];

	qr[0]=q[3];
	qr[1]=q[2];
	qr[2]=q[1];
	qr[3]=q[0];

	Real angles[3] = {0};

	quaternion2Euler( qr, angles );

	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}

/* this intuitive way is not successful!!
void mathMethods::randomOri( double * result )
{
	double q[4];
	double qnorm;

	q[0] = r.parkMiller();
	q[1] = r.parkMiller();
	q[2] = r.parkMiller();
	q[3] = r.parkMiller();

	bubbleSort( q, 4); //sorts ascending

	//normalize
	qnorm = sqrt(SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]));
	QUICKASSERT (qnorm > 0.0);
		
	q[0] = q[3] / qnorm;
	q[1] = q[2] / qnorm;
	q[2] = q[1] / qnorm;
	q[3] = q[0] / qnorm;

	double angles[3] = {0.0, 0.0, 0.0};
	quaternion2Euler ( q, angles);
	
	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}
*/

void mathMethods::randomOriShoemakeQuat( double * quat )
{
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//##MK::mind order 
	double q[4]={0,0,0,0};
	double qnorm;

	double s = r.parkMiller();
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2 * _PI_ * r.parkMiller();
	double theta2 = 2 * _PI_ * r.parkMiller();

	q[0]=fabs(sigma1*sin(theta1));
	q[1]=fabs(sigma1*cos(theta1));
	q[2]=fabs(sigma2*sin(theta2));
	q[3]=fabs(sigma2*cos(theta2));
	
	qnorm = sqrt(SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]));
	QUICKASSERT (qnorm > 0.0);

	/*
	bubbleSort( q,4 );
	Real qr[4];
	qr[0]=q[3];
	qr[1]=q[2];
	qr[2]=q[1];
	qr[3]=q[0];
	*/

	//normalize
    q[0] = q[0] / qnorm;
	q[1] = q[1] / qnorm;
	q[2] = q[2] / qnorm;
	q[3] = q[3] / qnorm;
	
	quat[0] = q[0];
	quat[1] = q[1];
	quat[2] = q[2];
	quat[3] = q[3];

}

void mathMethods::randomOriShoemakeEuler( double * result )
{
	//K. Shoemake, Graphic Gems III (editor D. Kirk) CalTech pp124-134
	//##MK::mind order 
	double q[4]={0,0,0,0};
	double qnorm;

	double s = r.parkMiller();
	double sigma1 = sqrt(1-s);
	double sigma2 = sqrt(s);
	double theta1 = 2 * _PI_ * r.parkMiller();
	double theta2 = 2 * _PI_ * r.parkMiller();

	q[0]=fabs(sigma1*sin(theta1));
	q[1]=fabs(sigma1*cos(theta1));
	q[2]=fabs(sigma2*sin(theta2));
	q[3]=fabs(sigma2*cos(theta2));
	
	qnorm = sqrt(SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]));
	QUICKASSERT (qnorm > 0.0);

	/*
	bubbleSort( q,4 );
	Real qr[4];
	qr[0]=q[3];
	qr[1]=q[2];
	qr[2]=q[1];
	qr[3]=q[0];
	*/

	//normalize
    q[0] = q[0] / qnorm;
	q[1] = q[1] / qnorm;
	q[2] = q[2] / qnorm;
	q[3] = q[3] / qnorm;

	Real angles[3] = {0};

	quaternion2Euler( q, angles );

	result[0] = angles[0];
	result[1] = angles[1];
	result[2] = angles[2];
}

void mathMethods::randomMisorientation( double theta, double* qr  ) // theta in radians
{
        double q[4]={0,0,0,0};
        double qcrit = cos(0.5*theta);

        while( q[3]<qcrit )
        {
                double s = r.parkMiller();
                double sigma1 = sqrt(1-s);
                double sigma2 = sqrt(s);
                double theta1 = 2 * _PI_ * r.parkMiller();
                double theta2 = 2 * _PI_ * r.parkMiller();
                q[0]=fabs(sigma2*cos(theta2));
                q[1]=fabs(sigma1*sin(theta1));
                q[2]=fabs(sigma1*cos(theta1));
                q[3]=fabs(sigma2*sin(theta2));
				//warum nicht wie in randomOriShoemake?
                bubbleSort( q,4 );
        }
        qr[0]=q[3];
        qr[1]=q[2];
        qr[2]=q[1];
        qr[3]=q[0];
}




void mathMethods::multiplyQuaternions( double *q, double* p, double* r )
{
        //mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
		//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
		//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
		r[0] = + q[0] *	p[0]	- q[1] * p[1]	- q[2] *	p[2]	- q[3] *	p[3];
        r[1] = + q[1] *	p[0]	+ q[0] * p[1]	- q[3] *	p[2]	+ q[2] *	p[3];
        r[2] = + q[2] *	p[0]	+ q[3] * p[1]	+ q[0] *	p[2]	- q[1] *	p[3];
        r[3] = + q[3] *	p[0]	- q[2] * p[1]	+ q[1] *	p[2]	+ q[0] *	p[3];
}

void mathMethods::multiplyQuaternions2( double *p, double *q, double* r)
{
		//what has been implemented in misorientationCubic in COReV2 is not a the standard definition of a quaternion multiplication!!!
		//but with qm1 as -q0, q1, q2, q3 applied to p it is the same as MTex
		r[0] = + p[0] * q[0] 	- p[1] * q[1] 		- p[2] * q[2] 			- p[3] * q[3];
		r[1] = + p[1] * q[0]	+ p[0] * q[1]		+ p[3] * q[2]			- p[2] * q[3];
		r[2] = + p[2] * q[0] 	- p[3] * q[1]		+ p[0] * q[2]			+ p[1] * q[3];
		r[3] = + p[3] * q[0]	+ p[2] * q[1]		- p[1] * q[2]			+ p[0] * q[3];
}

void mathMethods::misorientationQuaternionCubic( double* p, double* q, double* quat  )
{

	Real qm1[4];    //Inverse of quaternion q

	for( int i=0; i<4;i++) {
		qm1[i]=q[i];
        if( i > 0 ) qm1[i]*=-1;
    }
	//qm1[0] *= -1;

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1
	
	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3

	//Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688

	Real r0[6][4];    //There are 12 possible angles

    //Note: The notation r0 is due to the definition of the quaternion which lie
    //in the fundamental zone, this vector possesses the smallest angle, in such a way
    //that r0 is actually the scalar part of this quaternion

	double a = r[0]; 		
	double b = r[1]; 		
	double c = r[2]; 		
	double d = r[3];

	Real fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	Real rq[4];

	int mi;
	Real max=0.0;

	for( int i=0;i<6;i++ )						//Determing the quaternion with the maximal component and the component itself
		for( int j=0;j<4;j++ )
		{
			if( fabs(r0[i][j]) > max )
			{
				max=fabs(r0[i][j]);
				mi=i;
				//mj=j;
			}
		}

	rq[0] = fabs( r0[mi][0] );					//Desorientation requires all components positive
	rq[1] = fabs( r0[mi][1] );
	rq[2] = fabs( r0[mi][2] );
	rq[3] = fabs( r0[mi][3] );

	bubbleSort( rq,4 );						//Sorting into ascending order, because a desorientation in the SST
									//requires a quaternion with q0>=q1>=q2>=q3 which represents a minimal 
	quat[0] = rq[3];						//rotation angle and an axis fulfilling h>k>l
	quat[1] = rq[2];
	quat[2] = rq[1];
	quat[3] = rq[0];
}

double mathMethods::misorientationCubicOrigInv( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	int i;
	Real oria[3] = { pa1, Pa, pa2 };
	Real orib[3] = { pb1, Pb, pb2 };

	Real p[4];
	Real q[4];

	euler2quaternion( oria, p );
	euler2quaternion( orib, q );

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++)               //Inverting unit quaternion
        {
	        qm1[i]=q[i];
            if( i>0 ) qm1[i]*=-1;
        }
	//qm1[0] *= -1;

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

        multiplyQuaternions( p, qm1, r );

        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion

	Real fac = 1 / sqrt( 2.0 ); //0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];

	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 ) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}

double mathMethods::misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	Real oria[3] = { pa1, Pa, pa2 };
	Real orib[3] = { pb1, Pb, pb2 };

	Real p[4];
	Real q[4];

	euler2quaternion( oria, p );
	euler2quaternion( orib, q );

	Real qm1[4];    //Inverse of quaternion q

	for( int i=0; i<4;i++) {
		qm1[i]=q[i];
        if( i > 0 ) qm1[i]*=-1;
    }
	//qm1[0] *= -1;

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1
	
	multiplyQuaternions( qm1, p, r ); //MK:: was p post applied to qm1 but that was not consistent with MTex, in particular not with -(q^-1) = -q0,q1,q2,q3

    //Now, we have to determine the smallest angle, following Grimmer H, Acta Cryst 1974, A30 pp685-688

	Real r0[6][4];    //There are 12 possible angles

    //Note: The notation r0 is due to the definition of the quaternion which lie
    //in the fundamental zone, this vector possesses the smallest angle, in such a way
    //that r0 is actually the scalar part of this quaternion

	double a = r[0]; 		
	double b = r[1]; 		
	double c = r[2]; 		
	double d = r[3];

	Real fac = 1.0 / sqrt( 2.0 ); //0.70710678;

	//six fundamental quaternions
	r0[0][0] =(a+b)*fac; 		r0[0][1]=(a-b)*fac; 		r0[0][2] = (c+d)*fac; 		r0[0][3] = (c-d)*fac;
	r0[1][0] =(a+c)*fac; 		r0[1][1]=(a-c)*fac;	 		r0[1][2] = (b+d)*fac; 		r0[1][3] = (b-d)*fac;
	r0[2][0] =(a+d)*fac; 		r0[2][1]=(a-d)*fac; 		r0[2][2] = (b+c)*fac; 		r0[2][3] = (b-c)*fac;
	r0[3][0] =(a+b+c+d)*0.5; 	r0[3][1]=(a+b-c-d)*0.5; 	r0[3][2] = (a-b+c-d)*0.5; 	r0[3][3] = (a-b-c+d)*0.5;
	r0[4][0] =(a+b+c-d)*0.5; 	r0[4][1]=(a+b-c+d)*0.5; 	r0[4][2] = (a-b+c+d)*0.5; 	r0[4][3] = (a-b-c-d)*0.5;
	r0[5][0] = a;				r0[5][1] = b;				r0[5][2] = c;				r0[5][3] = d;

	Real omega = 0.0;

	for(int i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 ) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}


double mathMethods::misorientationCubicCOReV2( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	int i;

/*	pa1*=_PI_/180; Pa*=_PI_/180; pa2*=_PI_/180;
	pb1*=_PI_/180; Pb*=_PI_/180; pb2*=_PI_/180;*/

	Real p1 = pa1;  Real p12 = pb1;
	Real t = Pa;    Real t2 = Pb;
	Real p2 = pa2;  Real p22 = pb2;

	Real e1[3] = {p1, t, p2};
	Real e2[3] = {p12, t2, p22};
	
	// Quaternions
	Real p[4], q[4];

	euler2quaternion( e1, p );
	euler2quaternion( e2, q );

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++) qm1[i]=q[i]; //Copy quaternion

	qm1[0] *= -1;   //Inverting unit quaternion

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

	r[0]=p[0]*qm1[0]-p[1]*qm1[1]-p[2]*qm1[2]-p[3]*qm1[3];
	r[1]=p[1]*qm1[0]+p[0]*qm1[1]-p[2]*qm1[3]+p[3]*qm1[2];
	r[2]=p[2]*qm1[0]+p[0]*qm1[2]-p[3]*qm1[1]+p[1]*qm1[3];
	r[3]=p[3]*qm1[0]+p[0]*qm1[3]-p[1]*qm1[2]+p[2]*qm1[1];

        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion


	Real fac=0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];
	//these are the right fundamental quaternions of equivalent misorientations as shown by Grimmer, H 1974 Acta Cryst A30 685-688


	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 )
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}

double mathMethods::misorientationCubicCorrectCOReV2InvertAndMult( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 )
{
	int i;

/*	pa1*=_PI_/180; Pa*=_PI_/180; pa2*=_PI_/180;
	pb1*=_PI_/180; Pb*=_PI_/180; pb2*=_PI_/180;*/

	Real p1 = pa1;  Real p12 = pb1;
	Real t = Pa;    Real t2 = Pb;
	Real p2 = pa2;  Real p22 = pb2;

	Real e1[3] = {p1, t, p2};
	Real e2[3] = {p12, t2, p22};
	
	// Quaternions
	Real p[4], q[4];

	euler2quaternion( e1, p );
	euler2quaternion( e2, q );

	Real qm1[4];    //Inverse of quaternion q

	for(i=0;i<4;i++) {
		qm1[i]=q[i];
		if ( i > 0) qm1[i] *= -1;
	}

	//qm1[0] *= -1;

	Real r[4];     //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions(qm1, p , r);
	/*
	r[0]=p[0]*qm1[0]-p[1]*qm1[1]-p[2]*qm1[2]-p[3]*qm1[3];
	r[1]=p[1]*qm1[0]+p[0]*qm1[1]-p[2]*qm1[3]+p[3]*qm1[2];
	r[2]=p[2]*qm1[0]+p[0]*qm1[2]-p[3]*qm1[1]+p[1]*qm1[3];
	r[3]=p[3]*qm1[0]+p[0]*qm1[3]-p[1]*qm1[2]+p[2]*qm1[1];
	*/
        //Now, we have to determine the smallest angle.

	Real r0[6][4];    //There are 12 possible angles

        //Note: The notation r0 is due to the definition of the quaternion which lie
        //in the fundamental zone, this vector possesses the smallest angle, in such a way
        //that r0 is actually the scalar part of this quaternion

	Real fac=0.70710678;

	r0[0][0]=(r[0]+r[1])*fac; r0[0][1]=(r[0]-r[1])*fac; r0[0][2]=(r[2]+r[3])*fac; r0[0][3]=(r[2]-r[3])*fac;
	r0[1][0]=(r[0]+r[2])*fac; r0[1][1]=(r[0]-r[2])*fac; r0[1][2]=(r[1]+r[3])*fac; r0[1][3]=(r[1]-r[3])*fac;
	r0[2][0]=(r[0]+r[3])*fac; r0[2][1]=(r[0]-r[3])*fac; r0[2][2]=(r[1]+r[2])*fac; r0[2][3]=(r[1]-r[2])*fac;
	r0[3][0]=(r[0]+r[1]+r[2]+r[3])*0.5; r0[3][1]=(r[0]+r[1]-r[2]-r[3])*0.5; r0[3][2]=(r[0]-r[1]+r[2]-r[3])*0.5; r0[3][3]=(r[0]-r[1]-r[2]+r[3])*0.5;
	r0[4][0]=(r[0]+r[1]+r[2]-r[3])*0.5; r0[4][1]=(r[0]+r[1]-r[2]+r[3])*0.5; r0[4][2]=(r[0]-r[1]+r[2]+r[3])*0.5; r0[4][3]=(r[0]-r[1]-r[2]-r[3])*0.5;
	r0[5][0]=r[0];r0[5][1]=r[1];r0[5][2]=r[2];r0[5][3]=r[3];
	//these are the right fundamental quaternions of equivalent misorientations as shown by Grimmer, H 1974 Acta Cryst A30 685-688


	Real omega=0.0;

	for(i=0;i<6;i++)
		for( int j=0;j<4;j++ )
			if( fabs(r0[i][j]) > omega )
				omega=fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if( omega > 1.0 )
		omega = (Real) (int) omega;

	omega=2*acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}




void mathMethods::newOrientationFromReferenceFixedAngularCone(double * oriOri, double maxDev,double angle, double u, double v, double w,double * newOri)
{
        Real qr[4];
        Real ori[4], qideal[4];

        Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

        euler2quaternion( oriOri, qideal );
        Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

        randomMisorientationAxisConsidered(  qref, qr, maxDev  );
        multiplyQuaternions( qr,qideal,ori );

        Real euler[3];
        quaternion2Euler( ori, euler );

        newOri[0] = euler[0];
        newOri[1] = euler[1];
        newOri[2] = euler[2];
}

void mathMethods::rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri ) //Angle in radians
{
        Real ori[4], qideal[4];
        Real _norm = 1.0 / sqrt(SQR(u)+SQR(v)+SQR(w));

        euler2quaternion( oriOri, qideal );

        Real qref[4] = { cos( 0.5* angle ), u * _norm * sin( 0.5 * angle ), v * _norm * sin( 0.5 * angle ), w * _norm * sin( 0.5 * angle ) };

        multiplyQuaternions( qref, qideal, ori );

        double euler[3] = {0};

        quaternion2Euler( ori, euler );

        newOri[0] = euler[0];
        newOri[1] = euler[1];
        newOri[2] = euler[2];
}

void mathMethods::euler2quaternion( double * euler, double * q )
{
	/*20130326MK convention: Bunge ZXZ, which represents the (3,1,3) case analyzed in: Diebel 2006 
	Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, James Diebel, Stanford University, Stanford, California 94301-9010, Email: diebel@stanford.edu
	//cos(a+b) = c(a+b) = cacb-sasb
	//cos(a-b) = c(a-b) = cacb+sasb
	//sin(a+b) = s(a+b) = sacb+casb
	//sin(a-b) = s(a-b) = sacb-casb
	 */

	double p1 = euler[0];
	double t  = euler[1];
	double p2 = euler[2];

	Real co1 = cos(t/2);
	Real s1 = sin(t/2);

	Real p[4] = {co1*cos((p1+p2)/2),s1*cos((p1-p2)/2),s1*sin((p1-p2)/2),co1*sin((p1+p2)/2)}; //applying sin, cos addition theorems

	//double test[4]={0};
	//quaternion2Euler( p,test);

	q[0] = p[0];
	q[1] = p[1];
	q[2] = p[2];
	q[3] = p[3];
}

void mathMethods::quaternion2Euler( const double * quat, double * euler )
{
	//convention: Bunge, ZXZ, equal to case (3,1,3) as 
	//  analyzed in Diebel, James, 2006: 
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	//Gimbal lock situation analyzed following the line of Melcher et. al., Conversion of EBSD data by a quaternion based algorithm....
	//TECHNISCHE MECHANIK, 30, 4, (2010), 401 � 413
	//dont forget to define QUAT2EUL_ETA 1e-20

	Real q0 = quat[0];
	Real q1 = quat[1];
	Real q2 = quat[2];
	Real q3 = quat[3];
	Real PHI, sP, phi1, phi2;

	Real cosPHI = SQR(q3)-SQR(q2)-SQR(q1)+SQR(q0);

	Real y0 =	2*q1*q3	-	2*q0*q2;
	Real x0 =	2*q2*q3	+	2*q0*q1;
	Real y1 =	2*q1*q3	+	2*q0*q2;
	Real x1 = -	2*q2*q3	+	2*q0*q1;

	if( cosPHI > 1. ) cosPHI = 1.;

	if( SQR( 1. - cosPHI ) <= QUAT2EUL_ETA ) 
		PHI = 0.;
	else
		PHI = acos( cosPHI );

	sP = sin(PHI); //handle the gimbal lock situation that arouses as a quarternion does not define a Bunge Euler angle uniquely

	if( sP != 0 )
	{
		phi2 = atan2( y0 / sP, x0 / sP );
		phi1 = atan2( y1 / sP, x1 / sP );
	}else
	{
		phi1 = atan2( (2*q1*q2+2*q0*q3),SQR(q0)+SQR(q1)-SQR(q2)-SQR(q3) );
		phi2 = 0.;
	}
	
	//without additional sample and crystal symmetry the Euler space is symmetric to 0 <= phi1 <= 2*_PI_ , 0 <= PHI <= _PI_, 0 <= phi2 <= 2*_PI_
	if (phi1 < 0.0)
		phi1 += 2 * _PI_;
	if (phi2 < 0.0)
		phi2 += 2 * _PI_;

	euler[2] = phi2; //following the notation order used in Diebel, James 2006
	euler[1] = PHI;
	euler[0] = phi1;  
	
}


void mathMethods::quaternion2rodrigues( double * q, double * rodrigues )
{
	//convention: Bunge, ZXZ, equal to case (3,1,3) as 
	//  analyzed in Diebel, James, 2006: 
	//  Representing Attitude: Euler Angles, Unit Quaternions and Rotation Vectors
	
	//###implementation still due
	
	/*
	double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );
	double qq[4] = {q[0] / qnorm, q[1] / qnorm, q[2] / qnorm, q[3] / qnorm};

	double angleomega = 2 * acos( qq[0] );
	double s2 = sin(angle/2);
	*/
	
	rodrigues[0] = q[0]; //tan(angle/2);
	rodrigues[1] = q[1]; // / qnorm;
	rodrigues[2] = q[2]; // / qnorm;
	rodrigues[3] = q[3]; // / qnorm;
	
	//please contact me, if at some point this representation becomes necessary! kuehbach@imm.rwth-aachen.de
}

void mathMethods::newOrientationFromReference( double *bunge, double deviation, double *newOri )
{
	Real qrndmisori[4];
	Real rotated[4], qbunge[4];

	euler2quaternion( bunge, qbunge );

	randomMisorientation( deviation, qrndmisori  ); //this is really a problem when deviation approaches zero
	multiplyQuaternions( qrndmisori , qbunge, rotated );

	Real newEuler[3];
	quaternion2Euler( rotated, newEuler );

	newOri[0] = newEuler[0];
	newOri[1] = newEuler[1];
	newOri[2] = newEuler[2];
}

void mathMethods::devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb)
{
	//blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}

void mathMethods::devtoaxisEuler2RGB( double *bunge, double *uvw, double maxDev, unsigned char *rgb)
{
	//blue channel stretch from 0.0 to maxDev in radian, disorientation to axis of misorientation
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}

#define X	0
#define Y	1
#define Z	2

void mathMethods::patalaQuat2RGB ( double *q, unsigned char *rgb)
{
/*
	//Patala, Schuh unique coloring of misorientations via HSV2RGB
	//Diebel, James Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, 2006 Stanford University
	//Patala, Schuh Acta Materialia 59 (2011) 554�562
	//HSV2RGB by Gonzalez 
	
	//first get axis-angle parameterization of unit quaternion
	double p[4];
	double n[3]; //nx, ny, nz
	double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );
	//renormalization
	p[0] = q[0] / qnorm;
	p[1] = q[1] / qnorm;
	p[2] = q[2] / qnorm;
	p[3] = q[3] / qnorm;
	double _pnorm123 = (1.0 / sqrt ( 1.0 - SQR(p[0]) ) );
	
	double omega = 2*acos(p[0]); //no singularities in arccos(angle)
	n[X] = p[1] * _pnorm123;
	n[Y] = p[2] * _pnorm123;
	n[Z] = p[3] * _pnorm123;
	double tanomegahalf = tan(omega / 2.0);
	//coloring scheme is for misorientations only, thus q is a MISORIENTATION quaternion
	//then for m-3m fabs(omega) is bounded by 62.8/180*_PI_
	
	//step 1
	double xyz[3];
	xyz[X] = tanomegahalf * n[X]; xyz[Y] = tanomegahalf * n[Y]; xyz[Z] = tanomegahalf * n[Z];
	
	//step 2
	double xyz1[3];
	QUICKASSERT ( xyz[X] > 0.0 );
	QUICKASSERT ( xyz[Y] > 0.0 );
	//check boundness of yz
	if ( ( xyz[X] >= (1.0 / 3.0) ) && ( tan(xyz[Z] / xyz[Y]) >= ( ( 1 - 2*xyz[X] ) / xyz[X] ) ) {
		xyz1[X] = xyz[X]; 	
		xyz1[Y] = ( ( xyz[X] * (xyz[Y] + xyz[Z]) ) / (1.0 - xyz[X]) );
		xyz1[Z] = ( ( xyz[X] * xyz[Z] * (xyz[Y] + xyz[Z]) ) / ( xyz[Y] * (1.0 - xyz[X]) ) );
	}
	else {
		xyz1[X] = xyz[X];
		xyz1[Y] = xyz[Y];
		xyz1[Z] = xyz[Z];
	}
		
	//###MK::construction site, please let me know if this becomes necessary	
*/ 
	
	rgb[0] = RGBRANGE;
	rgb[1] = RGBRANGE;
	rgb[2] = RGBRANGE;
}



