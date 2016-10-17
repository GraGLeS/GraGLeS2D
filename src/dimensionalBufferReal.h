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
#ifndef		__DIMENSIONAL_BUFFER_REAL__
#define		__DIMENSIONAL_BUFFER_REAL__

#include "dimensionalBuffer.h"
#include "dataPrecision.h"

/*!
 * \class DimensionalBufferReal
 * \brief Class containing the distance value to the contour at each point managed by the LSbox.<br>
 * The buffer may contain either <b>float</b> or <b>double</b> values depending on the definition of <b>dataprecision</b>.
 */
class DimensionalBufferReal : public DimensionalBuffer<dataprecision>
{
public:

	DimensionalBufferReal() : DimensionalBuffer<dataprecision>()

		{}
	DimensionalBufferReal(unsigned int upperLeftX, unsigned int upperLeftY,
					  unsigned int lowerRightX, unsigned int lowerRightY) :
						  DimensionalBuffer(upperLeftX, upperLeftY, lowerRightX, lowerRightY)
	{}
	/*!
	 * \brief This method clamps all contained values in the given range.
	 */
	void clampValues(double minimumValue, double maximumValue)
	{
		for(unsigned int i=0; i<m_values.size(); i++)
		{
			if ( m_values[i] < minimumValue )
			{
				m_values[i] = minimumValue; continue;
			}
			if (m_values[i] > maximumValue )
			{
				m_values[i] = maximumValue; continue;
			}
		}
	}
};
#endif		//__DIMENSIONAL_BUFFER_REAL__
