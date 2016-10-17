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
#ifndef		__GRAIN_JUNCTION__
#define		__GRAIN_JUNCTION__

#include "spoint.h"
#include "Settings.h"

/*!
 * \enum E_JUNCTION_TYPE
 * \brief Enumeration that is used to identify the type of junction described by the
 * structure.
 */
enum E_JUNCTION_TYPE
{
	E_TRIPPLE_JUNCTION = 3,
	E_QUADRUPLE_JUNCTION,
	E_INVALID_JUNCTION
};

#define MAXIMUM_JUNCTION_ORDER	4
/*!
 * \struct GrainJunction
 * \brief Structure that encapsulates a grain junction. It carries information about the
 * order of the junction as well as all the grains active at that junction.
 */

class grainhdl;
class LSbox;
struct GrainJunction
{
	E_JUNCTION_TYPE	junction_type;
	SPoint coordinates;
	unsigned int grains[MAXIMUM_JUNCTION_ORDER];
	int	contourSegment;
	double lambda;
	static grainhdl* handler;
	GrainJunction(unsigned int* boxes, int count, SPoint coords) : coordinates(coords), contourSegment(-1), lambda(0)
	{
		if(count == 3)
			junction_type = E_TRIPPLE_JUNCTION;
		else if(count == 4)
			junction_type = E_QUADRUPLE_JUNCTION;
		else
			junction_type = E_INVALID_JUNCTION;
		for(int i=0; i < count && i < MAXIMUM_JUNCTION_ORDER; i++)
		{
			grains[i] = boxes[i];
		}
	}
	double getWeight(LSbox* me);
};

#endif		//__GRAIN_JUNCTION__
