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

#ifndef		__GRAHAM_SCAN_ALGORITHM__
#define		__GRAHAM_SCAN_ALGORITHM__
#include "voro++/include/voro++/voro++.hh"
#include "spoint.h"
#include <vector>
using namespace std;

/*!
 * \class GrahamScan
 * \brief Class that implements the Graham Scan algorithm. Can be used to extract a grain boundary
 * from a voro++ structure. The points in the resulting boundary are sorted counter clockwise.
 */

class GrahamScan
{
	/*!
	 * \struct ComparisonFunctor
	 * \brief Functor that is used to compare the ordering of three points (clockwise or counter-clockwise)
	 */
	struct ComparisonFunctor{
		const SPoint*	compareReference;
		bool operator() (const SPoint& lhs, const SPoint& rhs)
		{
			double val = (lhs.y - compareReference->y) * (rhs.x - lhs.x) -
						  (lhs.x - compareReference->x) * (rhs.y - lhs.y);
			if(val == 0)
			{
				return compareReference->squaredDistanceTo(lhs) < compareReference->squaredDistanceTo(rhs);
			}
			else
			{
				return val < 0;
			}
		}
	};

public:
	/*!
	* \brief This constructor is used when a contour is to be extracted from a voro++ structure.
	*/
	GrahamScan(voro::voronoicell_neighbor& voro_cell, unsigned int cellID, double* partPos);
	/*!
	* \brief This constructor is used when a contour is to be extracted from a point cloud structure.
	*/
	GrahamScan(vector<SPoint>& point_cloud);
	/*!
	* \brief This method performs the Grahm Scan and outputs the result in the provided vector.
	*/
	void generateCovnexHull(vector<SPoint>& output_hull);
private:
	vector<SPoint>				m_sortedPoints;
};
#endif		//__GRAHAM_SCAN_ALGORITHM__
