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
#ifndef __POOLED_DIMENSIONAL_BUFFER_DOUBLE__
#define __POOLED_DIMENSIONAL_BUFFER_DOUBLE__


/*!
 * \class PooledDimensionalBufferDouble
 * \brief Class that can create a dimensional buffer from an already allocated memory. The
 * class is not responsible for freeing the memory.
 */
class PooledDimensionalBufferDouble
{
public:
	PooledDimensionalBufferDouble(char* pool, unsigned int size,
			unsigned int upperLeftX, unsigned int upperLeftY,
			unsigned int lowerRightX, unsigned int lowerRightY) :
				m_xMin(upperLeftX), m_xMax(lowerRightX), m_yMin(upperLeftY), m_yMax(lowerRightY),
				m_pool(pool), m_poolSize(size)
	{
	}
	double getValueAt(unsigned int row, unsigned int column)
	{
		double* pointer = (double*) m_pool;
		return pointer[(row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)];
	}

	void setValueAt(unsigned int row, unsigned int column, double value)
	{
		double* pointer = (double*) m_pool;
		pointer[(row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin)] = value;
	}
	void clearValues(double value)
	{
		for(unsigned int i = 0; i < m_poolSize / sizeof(double); i++)
		{
			m_pool[i] = value;
		}
	}

private:
	int 	m_xMin;
	int 	m_xMax;
	int 	m_yMin;
	int 	m_yMax;
	char*	m_pool;
	int 	m_poolSize;

};
#endif 	//__POOLED_DIMENSIONAL_BUFFER_DOUBLE__
