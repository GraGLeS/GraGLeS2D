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
#ifndef _random_h_
#define _random_h_

#include <math.h>
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-30
#define RNMX (1.0-EPS)
#define ISNAN(a) ((a) != (a))
#define _PI_ 3.1415926535897932384626433832795
#define EPSILON 0.01
#define EPS1 0.001
#define EPS2 1.0e-8
#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MYHASH(a,b) ( (((uint_fast64_t) max(a,b)) << 32) | ((uint_fast64_t) min(a,b)) )
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)
#define XOR(a,b) ((!!a) != (!!b))

class randomClass
{
public:
	randomClass(long seed = 0);

	void init(long sd) { seed = sd; }
	//void init(long sd) { seed = (long) time(NULL);} hat nicht funktioniert

	long integer(void);
	double real(void);
/*	double gauss(void);
	double wls(void);*/
	double parkMiller( void );
	double randomInterval( double rmin, double rmax );

  private:
		long seed;
};

typedef randomClass *randomClassP;

#endif
