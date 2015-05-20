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
#ifndef		__CONTOUR_SECTOR__
#define		__CONTOUR_SECTOR__
#include "junction.h"
#include "dimensionalBufferIDLocal.h"
#include <fstream>
#include <vector>
using namespace std;

/*!
 * \enum E_CONTOUR_TYPE
 * \brief Enumeration that is used to distinguish between the types of contour sectors.
 * Shall be replaced by a polymorphic implementation.
 */
enum E_CONTOUR_TYPE
{
	E_INTERPOLATION_SECTOR,
	E_CONSTANT_JUNCTION_SECTOR,
	E_INVALID_TYPE_SECTOR
};
class LSbox;

/*!
 * \class ContourSector
 * \brief Class that encapsulates a Contour Sector.
 *
 * This class handles the different  types of sectors defined i the convolution correction
 * method. Each sector is represented by a start point and an end point on the contour line.
 * The winding of the contour line is counter clockwise and thus a left end point and a right
 * end point can be defined. The Constant Sectors are an area around triple and quadruple junctions
 * and can contain on or more of them. Their weight is an average of the weights of all the junctions
 * in the sector. The Interpolating Sectors are areas that are immediately following the constant sectors
 * and provide a smooth transition between the weight of a constant sector and the weight of the
 * neighboring region (another sector or a grain boundary without any junctions).
 */
class ContourSector
{
public:
	/*!
	* \brief Basic constructor. Requires an initial junction. Can be used only for constant sectors.
	*/
	ContourSector(GrainJunction* initialJunction);
	/*!
	* \brief Constructor for arbitrary sectors.
	*/
	ContourSector(E_CONTOUR_TYPE type = E_INTERPOLATION_SECTOR);
	/*!
	* \brief Adds another junction to a constant sector.
	*/
	bool			mergeWith(GrainJunction* junction);
	/*!
	* \brief Checks if the given segment (an ID of a point on a contour line) of a grain boundary is within this sector.
	* \param contour_grain The contour line of the grain to which this contour belongs.
	* \param segment The id of the point on the contour line (segment) to be checked.
	*/
	bool			isSegmentWithinSector(vector<SPoint>& contour_grain, int segment) const;
	/*!
	* \brief Returns the weight for a given point. The point is given as an segment and an offset.
	*/
	double			getWeight(LSbox* owner, vector<SPoint>* contour_grain = NULL, int segment=0, double lambda=0.0) const;
	/*!
	* \brief Sets the left end of the contour.
	*/
	void 			setLeftContourPoint(int ID);
	/*!
	* \brief Sets the right end of the contour.
	*/
	void 			setRightContourPoint(int ID);
	/*!
	* \brief Sets the right end of the contour (only for interpolating sectors)
	*/
	void			setRightSectorBoundary(ContourSector* sector);
	/*!
	* \brief Sets the left end of the contour (only for interpolating sectors)
	*/
	void			setLeftSectorBoundary(ContourSector* sector);
	/*!
	* \brief Gets the right end of the contour (only for interpolating sectors)
	*/
	ContourSector*  getRightSectorBoundary() const {return m_rightSectorBoundary;}
	/*!
	* \brief Gets the left end of the contour (only for interpolating sectors)
	*/
	ContourSector*  getLeftSectorBoundary() const {return m_leftSectorBoundary;}
	/*!
	* \brief Gets the right end of the contour.
	*/
	inline int		getRightContourPoint() const {return m_rightContourPointID;}
	/*!
	* \brief Gets the left end of the contour.
	*/
	inline int		getLeftContourPoint() const {return m_leftContourPointID;}
	/*!
	* \brief Returns an interpolating sector on the left side of the sector (Only for constant sectors).
	*/
	ContourSector 	generateLeftInterpolatingContour(vector<SPoint>& contour_grain, double h);
	/*!
	* \brief Returns an interpolating sector on the right side of the sector (Only for constant sectors).
	*/
	ContourSector 	generateRightInterpolatingContour(vector<SPoint>& contour_grain, double h);
	/*!
	* \brief Computes the total length of the sector.
	*/
	void			recalculateGrainBoundaryLength(vector<SPoint>& contour_grain, double h);
	/*!
	* \brief Returns an interpolating sector that connects two constant sectors.
	*/
	ContourSector	generateInterpolatingSector(ContourSector& other,vector<SPoint>& contour_grain, double h);
	/*!
	* \brief Outputs the sector in a plot-friendly manner
	*/
	void			debugPrintSector(vector<SPoint>& contour_grain, ofstream& ofs);
	/*!
	* \brief Returns an array of the Junctions that are interacting in this sector.
	*/
	const std::vector<GrainJunction*>& getJunctions() {return m_junctions;}

private:
	E_CONTOUR_TYPE				m_sectorType;
	int 						m_leftContourPointID;
	int 						m_rightContourPointID;
	std::vector<GrainJunction*>	m_junctions;
	ContourSector*				m_leftSectorBoundary;
	ContourSector*				m_rightSectorBoundary;
	double						m_grainBoundaryLength;
};

#endif		//__CONTOUR_SECTOR__
