#include "random.h"
#include "mymath.h"

randomClass::randomClass(long sd)
{
	if( sd )	seed = sd;
		else	seed = (long) time(NULL);
}

double randomClass::real(void)
{
	seed *= 22695477L ;
	seed += 1L;

	return ( (65536.0 * 32768.0) + seed) * (1.0 / ( 65536.0 * 65536.0 ));

/*	long ls = seed & 0xffff;
	long hs = seed >> 16;

	return (double)ls*(1.0/65536.0) + (double)hs*(1.0/(65536.0*65536.0));   */

}

double randomClass::randomInterval( Real rmin, Real rmax )
{
	if( rmax < rmin ) return 0;
	Real r = parkMiller();
	return rmin + r * ( rmax - rmin );
}

double randomClass::parkMiller( void )
{
// Random Number Generator Park-Miller-Bays-Durham "Numerical Recipes in C"
        int j;
        long k;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if( seed <= 0 || !iy )
        {
                if( -seed < 1 ) seed=1;
                else seed = -seed;
                for( j=NTAB+7;j>=0;j-- )
                {
                        k=seed/IQ;
                        seed=IA*(seed-k*IQ)-IR*k;
                        if( seed < 0 ) seed += IM;
                        if( j < NTAB ) iv[j] = seed;
                }
                iy=iv[0];
        }
        k=seed/IQ;
        seed=IA*(seed-k*IQ)-IR*k;
        if( seed < 0 ) seed += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = seed;
        temp = AM*iy;
/*  if( (temp=AM*iy) > RNMX ) return RNMX;
  else*/ return temp;
}

long randomClass::integer(void)
{
	seed *= 22695477L ;
	seed += 1L;
	return seed & 0xffff;
}


/*double randomClass::gauss(void)		// 1/sqrt(_2PI_) * exp(-r*r) ;
{
	double u,v,w,x,sum;
 
	u = real();

	if( u <= 0.8638 )
	{
		v = - real() * 2.0;
		u *= 2.3153508;
		v += real()*2.0 -1.0;
		return u + v;
	}

	if( u <= 0.9745 )
	{
		v = -1.5 * real();
		u -= 0.8638;
		u *= (9.0334237 * 1.5);
		return v + u;
	}

	if( u <= 0.9973002 )
	{
		do
		{
			x = 6.0 * real() - 3.0;
			v = fabs(x);
			u = real();
			w = 3.0 - v;
			w *= w;
			w *= 6.6313339;
			sum = 0;
			if( v < 1.5 ) sum = 6.0432809 * (1.5 - v);
			if( v < 1.0  ) sum += 13.2626678 * (3.0 - v*v) - w;
			v *= -0.5 * v;
			u += w + sum;
		}
		 while( u > 49.0024445 * exp(v) );
		return x;
	}
   
	do 
	{
		do
		{
			w = real() ;
		} 
		 while( w == 0 );
 
		x = 4.5 - log(w);
		v = real();
		v *= v * x;
	}while( v > 4.5 ); 
	  
	x = sqrt( x * 2.0 );   
	if(u > 0.9986501 ) return x;
	if(u < 0.9986501 ) return -x;
	return 0;
}*/

/*double randomClass::wls(void)
{
	double r,z,w;
	const double max = 2.16;
	
	do
	{
		r = real() * 1.5;
		w = WLS( r );
		z = real() * max;
			
	} while( z > w );
	 
	return r;
}*/

