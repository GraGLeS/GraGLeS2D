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

#ifndef GGLS_h
#define GGLS_h



// #define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define INTERIMVAL -sqrt(1. / PARTICLES)

#define IMAGEOUT true
#define mymin(x,y) ( (x)>(y) ? (y) : (x) )

// // switch between double or single precision!


#pragma once
#include <map>
#include <vector>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <list>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <fftw3.h>
#include <time.h>  
#include <omp.h>

// g++ -I/home/user/bin/R-devel/include -L/home/user/bin/R-devel/lib -lR embed.cpp

#include "voro++/include/voro++/voro++.hh"
#include "dataPrecision.h"
#include "box.h"
#include "grainhdl.h"
#include "mymath.h"
#include "SCORE_Random.h"
#include "dimensionalBuffer.h"
#include "marchingSquares.h"




#endif

