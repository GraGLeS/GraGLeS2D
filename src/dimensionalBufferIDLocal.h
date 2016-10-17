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
#ifndef __DIMENSIONAL_BUFFER_IDLOCAL__
#define __DIMENSIONAL_BUFFER_IDLOCAL__

#include <vector>
#include <map>
#include <iostream>

#include "dimensionalBuffer.h"
using namespace std;

#define ID_MAX_CHUNK_COUNT 1

class LSbox;

enum E_POSITIONS
{
	E_FIRST_POSITION = 1,
	E_SECOND_POSITION,
	E_THIRD_POSITION,
	E_FOURTH_POSITION,
	E_LAST_POSITION = -1
};
/*!
 * \struct IDChunk
 * \brief Deprecated structure that holds the pointers to the grains that are closest to a given point.
 */
struct IDChunk
{
	LSbox* local_chunks[ID_MAX_CHUNK_COUNT];
	unsigned int total_elements;
	IDChunk() : total_elements(0)
	{}
	IDChunk& operator=(const IDChunk& other)
	{
	    if (this != &other)
	    {
	    	total_elements = other.total_elements;
	    }
	    return *this;
	}
	bool insertAtPosition(E_POSITIONS position, LSbox* element)
	{
		if(position == E_LAST_POSITION)
		{
			if(total_elements >= ID_MAX_CHUNK_COUNT)
			{
				return false;
			}
			else
			{
				local_chunks[total_elements] = element;
				total_elements++;
			}
		}
		else
		{
			int actual_position = position - E_FIRST_POSITION;
			if(actual_position < 0 || actual_position >= ID_MAX_CHUNK_COUNT)
			{
				return false;
			}
			else
			{
				LSbox* element_to_write = element;
				for(unsigned int i=actual_position; i<total_elements+1 && i < ID_MAX_CHUNK_COUNT; i++)
				{
					LSbox* swap = local_chunks[i];
					local_chunks[i] = element_to_write;
					element_to_write = swap;
				}
			}
			total_elements++;
			if(total_elements > ID_MAX_CHUNK_COUNT)
				total_elements = ID_MAX_CHUNK_COUNT;
		}
		return true;
	}
	LSbox* getElementAt(unsigned int pos)
	{
		if(pos < total_elements)
			return local_chunks[pos];
		else
			return NULL;
	}
	inline void clear()
	{
		total_elements = 0;
	}
};

/*!
 * \struct IDChunkMinimal
 * \brief Structure used to store the ID of the closest grain to a specific point. Usually holds the ID of the
 * grain in which the point lies.
 */
struct IDChunkMinimal
{
	unsigned int grainID;
	void clear() {}
};
/*!
 * \class DimensionalBufferIDLocal
 * \brief Class containing the IDLocal information for each point of the grid of an
 * LSbox object. <br> For each point this class contains the ID of the closest grain i.e.
 * the grain in which this point belongs.
 */
class DimensionalBufferIDLocal	: public DimensionalBuffer<IDChunkMinimal>
{
public:
	/*!
	 * \brief This method clears all the data stored in the buffer.
	 */
	void clear()
	{
		for(auto& iterator : this->m_values)
			iterator.clear();
	}
};



#endif //__DIMENSIONAL_BUFFER_IDLOCAL__
