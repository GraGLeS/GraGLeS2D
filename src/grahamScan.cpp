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

#include "grahamScan.h"
#include <algorithm>
#include <set>

GrahamScan::GrahamScan(voro::voronoicell_neighbor& voro_cell, unsigned int cellID, double* partPos)
{
	vector<double> vv;
	double x1[2], x2[2];
	set<SPoint>	pointset;
	//Iterate over the voro++ structure and produce unique points on the contour.
	voro_cell.vertices(partPos[3*(cellID-1)],partPos[3*(cellID-1)+1],partPos[3*(cellID-1)+2],vv);
	for (int ii = 0; ii < voro_cell.p; ii++)
	{
		for (int jj = 0; jj < voro_cell.nu[ii]; jj++)
		{

			int k = voro_cell.ed[ii][jj];
			x1[0] = vv[3 * ii];
			x1[1] = vv[3 * ii + 1];
			x2[0] = vv[3 * k];
			x2[1] = vv[3 * k + 1];
			pointset.insert(SPoint(x1[1], x1[0], 0));
			pointset.insert(SPoint(x2[1], x2[0], 0));
		}
	}
	std::copy(pointset.begin(), pointset.end(), std::back_inserter(m_sortedPoints));

	ComparisonFunctor comparator;
	comparator.compareReference = &(*pointset.begin());

	vector<SPoint>::iterator start = m_sortedPoints.begin(); start++;
	sort(start, m_sortedPoints.end(), comparator);
}
GrahamScan::GrahamScan(vector<SPoint>& point_cloud)
{
	std::copy(point_cloud.begin(), point_cloud.end(), std::back_inserter(m_sortedPoints));
	sort(m_sortedPoints.begin(), m_sortedPoints.end());

	ComparisonFunctor comparator;
	comparator.compareReference = &m_sortedPoints[0];

	vector<SPoint>::iterator start = m_sortedPoints.begin(); start++;
	sort(start, m_sortedPoints.end(), comparator);
}

void GrahamScan::generateCovnexHull(vector<SPoint>& output_hull)
{
	output_hull.clear();
	output_hull.push_back(m_sortedPoints[0]);
	output_hull.push_back(m_sortedPoints[1]);
	output_hull.push_back(m_sortedPoints[2]);

	ComparisonFunctor compare;
	compare.compareReference = &output_hull[1];

	int curr_top = 2;
	for(unsigned int i=3; i<m_sortedPoints.size(); i++)
	{
		while(compare(output_hull[curr_top], m_sortedPoints[i]) == false)
		{
			curr_top--;
			output_hull.pop_back();
			compare.compareReference = &output_hull[curr_top-1];
		}
		output_hull.push_back(m_sortedPoints[i]);
		curr_top++;
		compare.compareReference = &output_hull[curr_top-1];
	}
	output_hull.push_back(output_hull[0]);
}
