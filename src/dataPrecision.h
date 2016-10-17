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
#ifndef __DATAPRECISION__
#define __DATAPRECISION__

#define PRECISION 0	//  0 = double 1 = single
#if PRECISION >0
  #define dataprecision float
  #define fftwp_malloc fftwf_malloc
  #define PooledDimensionalBuffer PooledDimensionalBufferReal
  #define fftwp_complex fftwf_complex
  #define fftwp_plan fftwf_plan
  #define fftwp_free fftwf_free
  #define fftw_destroy_planp fftwf_destroy_plan
#elif PRECISION <1
  #define dataprecision double
  #define fftwp_complex fftw_complex
  #define fftwp_plan fftw_plan
  #define fftwp_free fftw_free
  #define fftwp_malloc fftw_malloc
  #define fftw_destroy_planp fftw_destroy_plan
  #define PooledDimensionalBuffer PooledDimensionalBufferDouble
#endif

#endif
