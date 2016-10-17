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
#ifndef		__MARCHING_SQUARES_ALGORITHM__
#define		__MARCHING_SQUARES_ALGORITHM__

 
#include "dimensionalBufferReal.h"
#include "dimensionalBufferIDLocal.h"
#include "spoint.h"
#include "junction.h"
#include "dimensionalBuffer.h"
#include "ggLS.h"



class LSbox;

/*!
 * \enum E_CORNER_VALUES
 * \brief Enumeration that is used to mark the edges of a square that are inside the area
 * of the object whose contour is being traced by the marching squares algorithm.
 */
enum E_CORNER_VALUES
{
	TOP_LEFT = 1,
	TOP_RIGHT = 1<<1,
	BOTTOM_LEFT = 1<<2,
	BOTTOM_RIGHT = 1<<3
};
/*!
 * \enum E_MOVEMENT_DIRECTIONS
 * \brief Enumeration that describes the directions in which the marching square shall advance.
 */
enum E_MOVEMENT_DIRECTIONS
{
	UP,
	DOWN,
	LEFT,
	RIGHT,
	NO_DIRECTION
};
/*!
 * \class MarchingSquaresAlgorithm
 * \brief Class that implements the marching squares algorithm.
 *
 * This class handles 2 dimensional patterns. It uses the standard marching squares algorithm in
 * counter-clock-wise direction. Currently it needs to be improved for ambiguity points. The class has
 * virtual methods that can be overridden to customize the way the algorithm works. The basic
 * implementation is suited to the Levelset project.
 */
class MarchingSquaresAlgorithm
{
public:
	/*!
		 * \brief Basic constructor. Requires a dimensional that defines the shape of the object.
	*/
	MarchingSquaresAlgorithm(DimensionalBufferReal& distance_buffer, DimensionalBufferIDLocal& id_local, LSbox* current_grain);
	/*!
		 * \brief Basic destructor. Virtual so that the class can be easily inherited.
	*/
	~MarchingSquaresAlgorithm(){}
	/*!
		 * \brief This method generates the contour as a vector of SPoint objects. In this vector
		 * the first and the last point will be the same. Although this method can be overridden
		 * The core of the algorithm resides here.
	*/
	bool	generateContour(std::vector<SPoint>& contour_output, std::vector<GrainJunction>& junction_output);
	/*!
		 * \brief This method is used to determine whether a point at the specified row and column
		 * is actually inside the object.
	*/
	bool	isInside(int row, int column);
	/*!
		 * \brief This method generates a point by the current movement direction.
	*/
	SPoint	generatePoint(E_MOVEMENT_DIRECTIONS dir);
	/*!
		 * \brief This method generates a junction from the current square being inspected.
		 * Junctions always contain all grains including the current grain whose contour is
		 * being constructed.
	*/
	void	generateJunction(std::vector<GrainJunction>& junctions, int state_mask) const;
private:

	/*!
		 * \brief This method inserts a point in the output by first checking if the point
		 * is not already there. Prevents point duplication.
	*/
	void	insertPoint(std::vector<SPoint>& output, SPoint p);
	/*!
		 * \brief This method inserts an LSbox pointer in the specified array, by maintaining only distinct
		 * pointers. This is used to detect the order of the junction.
	*/
	inline void	insertDistinctInteger(unsigned int intger, unsigned int* array, int& elem_count) const;

	DimensionalBufferReal& 		m_DistanceBuffer;
	DimensionalBufferIDLocal&	m_IDLocal;
	LSbox*	m_CurrentGrain;

	int m_top;
	int m_bottom;
	int m_left;
	int m_right;

};
#endif		//__MARCHING_SQUARES_ALGORITHM__
