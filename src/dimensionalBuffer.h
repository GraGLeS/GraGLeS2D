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
#ifndef		__DIMENSIONAL_BUFFER__
#define		__DIMENSIONAL_BUFFER__

#include <vector>
#include <iostream>
#include <cmath>
#include "ExpandingVector.h"
#include <iostream>
/*!
 * \class DimensionalBuffer
 * \brief Class that encapsulates a dimensional buffer
 *
 * Dimensional buffer for ease in management of dimensional data. The data exists in space, where
 * the minimum x and y coordinates are 0 and the maximal coordinate value is bounded by the size
 * of signed 32 bit integer.
 */
template<class T>
class DimensionalBuffer {
public:
	/*!
	 * \brief Basic constructor. Just initializes the object with some default values.
	 */
	DimensionalBuffer() :
		m_xMin(0), m_xMax(1), m_yMin(0), m_yMax(1) {
	}
	/*!
	 * \brief Constructor, which receives the coordinates of the managed region.
	 *
	 * This constructor takes the boundaries of the managed region as input, and properly
	 * initializes the internal memory.
	 */
	DimensionalBuffer(unsigned int upperLeftX, unsigned int upperLeftY,
			unsigned int lowerRightX, unsigned int lowerRightY) :
		m_xMin(upperLeftX), m_xMax(lowerRightX), m_yMin(upperLeftY),
				m_yMax(lowerRightY) {
		resize(m_xMin, m_yMin, m_xMax, m_yMax);
	}
	/*!
	 * \brief Default destructor.
	 */
	~DimensionalBuffer() {
	}
	/*!
	 * \brief Method that returns the value at the given coordinates.
	 *
	 * This method retrieves the value at the specified coordinates. If such a value does not exist
	 * this method throws an out of bound exception.
	 * \param row the y coordinate of the element.
	 * \param column the x coordinate of the element.
	 */
	T& getValueAt(unsigned int row, unsigned int column) {
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		return m_values.at(
				(row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin));
	}
	/*!
	 * \brief Method that sets the value at the given coordinates.
	 *
	 * This method sets the value at the specified coordinates. If such a value does not exist
	 * this method throws an out of bound exception.
	 * \param row the y coordinate of the element.
	 * \param column the x coordinate of the element.
	 * \param value the value to be set.
	 */
	void setValueAt(unsigned int row, unsigned int column, T value) {
		//Will throw exception if accessed out of bound.
		//TODO: Analyze performance and replace with [] if needed.
		m_values.at((row - m_yMin) * (m_xMax - m_xMin) + (column - m_xMin))
				= value;
	}
	/*!
	 *\brief Method that checks whether a point lies within the dimensional buffer.
	 *
	 * \param row the y coordinate of the point to test.
	 * \param column the x coordinate of the point to test.
	 */
	bool isPointInside(unsigned int row, unsigned int column) const {
		return row >= m_yMin && row < m_yMax && column >= m_xMin && column
				< m_xMax;
	}
	/*!
	 * \brief This method resizes the dimensions.
	 *
	 * This method resizes the dimensions and properly manages the internal data.
	 * \param upperLeftX the desired new minimal x coordinate.
	 * \param upperLeftY the desired new minimal y coordinate.
	 * \param lowerRightX the desired new maximal x coordinate.
	 * \param lowerRightY the desired new maximal y coordinate.
	 */
	void resize(unsigned int upperLeftX, unsigned int upperLeftY,
			unsigned int lowerRightX, unsigned int lowerRightY) {
		m_xMin = upperLeftX;
		m_xMax = lowerRightX;
		m_yMin = upperLeftY;
		m_yMax = lowerRightY;
		m_values.resize((m_xMax - m_xMin) * (m_yMax - m_yMin));
	}

	/*!
	 * \brief This method resizes the dimensions.
	 *
	 * Resizes the current area to a square area and manages the internal memory.
	 * \param maximumLength The maximal value for the x or y coordinates.
	 */
	void resizeToSquare(unsigned int maximumLength) {
		int grid_size = maximumLength - 1;
		unsigned int height = m_yMax - m_yMin;
		unsigned int width = m_xMax - m_xMin;
		//First resize the rectangle to a square
		if (width == height)
			return;
		else if (width > height) {
			int diff = width - height;
			m_yMin -= diff / 2;
			m_yMax += diff / 2 + diff % 2;
		} else {
			int diff = height - width;
			m_xMin -= diff / 2;
			m_xMax += diff / 2 + diff % 2;
		}
		//Now move the square in the bounds if it has left them
		if (m_xMin < 0) {
			int delta = -(int) m_xMin;
			m_xMax += delta;
			m_xMin += delta;
		} else if (m_xMax > grid_size) {
			int delta = grid_size - m_xMax;
			m_xMax += delta;
			m_xMin += delta;
		}
		if (m_yMin < 0) {
			int delta = -(int) m_yMin;
			m_yMax += delta;
			m_yMin += delta;
		} else if (m_yMax > grid_size) {
			int delta = grid_size - m_yMax;
			m_yMax += delta;
			m_yMin += delta;
		}

		resize(m_xMin, m_yMin, m_xMax, m_yMax);
	}
	/*!
	 * \brief This method fills the area with the provided value.
	 */
	void clearValues(T value) {
		std::fill(m_values.begin(), m_values.end(), value);
	}
	/*!
	 * \brief This method returns the total memory allocated by this buffer.
	 */
	int getTotalMemoryUsed() const {
		return m_values.capacity() * sizeof(T);
	}
	/*!
	 * \brief This method returns the left boundary of the managed region.
	 */
	inline int getMinX() const {
		return m_xMin;
	}
	/*!
	 * \brief This method returns the right boundary of the managed region.
	 */
	inline int getMaxX() const {
		return m_xMax;
	}
	/*!
	 * \brief This method returns the bottom boundary of the managed region.
	 */
	inline int getMinY() const {
		return m_yMin;
	}
	/*!
	 * \brief This method returns the top boundary of the managed region.
	 */
	inline int getMaxY() const {
		return m_yMax;
	}
	/*!
	 * \brief This method returns a pointer to the actual data stored in the buffer.
	 */
	inline T* getRawData() {
		return &m_values[0];
	}

private:
	int m_xMin;
	int m_xMax;
	int m_yMin;
	int m_yMax;
protected:
	ExpandingVector<T> m_values;
};
#endif		//__DISTANCE_BUFFER__
