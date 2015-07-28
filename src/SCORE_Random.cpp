/*
    GraGLeS 2D A grain growth simulation utilizing level set approach

 	Copyright (C) 2015  Markus K{\"u"}bach, Luis Baralles-Mora

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

   	PRNG as implemented by L. A. Barrales-Mora from references (see functions), some explicit 32-bit datatypes where necessary
	SCORE automaton developed by M K{\"u}hbach, 2014/2015, for questions and details contact markus.kuehbach@rwth-aachen.de
 */

#include "SCORE_Random.h"

randomClass::randomClass(int32_t sd) 
{
	//Constructor calls always init but init can be called outside to re-initialize the seeds if desired
	init( sd );
}

void randomClass::init( int32_t sd )
{
	int32_t seed;
	
	if( sd == DEFAULT_SEED )	seed = (int32_t) time( NULL );
	
	mti = NMT+1;

	seedOPM = seed;
	seedMT = (uint32_t) seed;
	
	seedOLE1 = seedMT;
	seedOLE2 = (uint32_t) time(NULL);

	if( seedOPM < 0 ) seedOPM += MPM;

	sgenrand();
}

long randomClass::randomIntervalLongInclusive_parkMiller( long lmin, long lmax ) 
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = parkMiller();

	return (long) ( lmin + r * ( lmax - lmin ) );
}


long randomClass::randomIntervalLongInclusive_leEcuyer( long lmin, long lmax )
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = leEcuyer();

	return lmin + ( long ) (r * ( lmax - lmin  ) );
}


double randomClass::randomInterval( double rmin, double rmax ) 
{
	//MK::ok
	if( rmax < rmin ) return 0.0;
	double r = parkMiller();
	return rmin + r * ( rmax - rmin );
}


void randomClass::generateShuffledLongArray( long lmin, long lmax, long * result ) 
{
	//MK::ok
	if( lmax <= lmin ) { result = NULL; return; }

	for( long i = 0; i < ( lmax - lmin ); i++ ) result[i] = i + lmin;

	for( long i = 0; i < ( lmax - lmin ); i++ )
	{
		//swop element at a random pos with the one at i successively
		long pos = i + randomIntervalLongInclusive_parkMiller( 0,  lmax - lmin - i - 1 );
		long temp = result[i];
		result[i] = result[pos];
		result[pos] = temp;
	}
}

double randomClass::parkMiller( void )
{


	int32_t lo, hi, test;

	hi =  seedOPM / QPM;
	lo =  seedOPM - QPM * hi;

	test = APM * lo - RPM * hi;

	if( test > 0 )
		seedOPM = test;
	else
		seedOPM = test + MPM;
	
	double randr = ((double) seedOPM) / MPM;
	
	if( randr > RANDMAX ) return RANDMAX;

	return randr;
}

void randomClass::sgenrand( void )	//Initialization of Mersenne twister
{
	mt[0] = seedMT & 0xffffffff;
	for( mti=1; mti<NMT; mti++ )
		mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

double randomClass::MersenneTwister( void )
{
/* Saito M, Matsumoto M. SIMD-oriented fast Mersenne twister: a 128-bit
   pseudorandom number generator, Monte-Carlo and Quasi-Monte Carlo Methods
   2006;2:607 */

//Original algorithm interval [0,1]; here modified for [0,1)
	
	uint32_t y;
	static uint32_t mag01[2]={0x0, MATRIX_A};

	if( mti >= NMT )
	{
		int kk;

		if( mti == NMT+1 )
		{
			seedMT = (uint32_t) time(NULL);
			sgenrand();
		}

		for( kk=0; kk<NMT-MMT; kk++ )
		{
			y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
			mt[kk] = mt[kk+MMT] ^ ( y >> 1 ) ^ mag01[y & 0x1];
		}

		for(;kk<NMT-1;kk++)
		{
			y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
			mt[kk] = mt[kk+(MMT-NMT)] ^ (y >> 1) ^ mag01[y & 0x1];
		}

		y = ( mt[NMT-1] & UPPER_MASK ) | ( mt[0] & LOWER_MASK );
		mt[NMT-1] = mt[MMT-1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);
	
	double randr = ( (double) y / (uint32_t)0xffffffff );
	
	if( randr > RANDMAX )	return RANDMAX;

	return randr;
}

double randomClass::leEcuyer( void )
{

/* L'Ecuyer P. Efficient and portable combined random number generators, Communications of the ACM 1998;31:742 */
// Interval (0,1)

	int32_t Z, k;

	while( seedOLE1 > M1-1 )
	{
		seedOLE1 = seedOLE1 / 2;
	}

	while( seedOLE2 > M2-1 )
	{
		seedOLE2 = seedOLE2 / 2;
	}

	k = seedOLE1 / D1RE;
	seedOLE1 = A1 * (seedOLE1 - k * D1RE) - k * M1RE;

	if( seedOLE1 < 0 )	seedOLE1 += M1;

	k = seedOLE2 / D2RE;
	seedOLE2 = A2 * (seedOLE2 - k * D2RE) - k * M2RE;

	if( seedOLE2 < 0 )	seedOLE2 += M2;

	Z = seedOLE1 - seedOLE2;

	if( Z < 1 ) Z += (M1-1);

	return Z * (1.0/M1);
	
}
