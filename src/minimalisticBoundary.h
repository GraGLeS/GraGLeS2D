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
#ifndef __MINIMALISTIC_BOUNDARY__
#define __MINIMALISTIC_BOUNDARY__

#include "spoint.h"
#include "contourSector.h"
#include <vector>
#include <set>
using namespace std;

class LSbox;

/*!
 * \class MinimalisticBoundary
 * \brief Class that reduces the boundary to the minimal amount of sectors required. Not yet
 * functional.
 */
class MinimalisticBoundary
{
private:
	vector<SPoint>			m_minimalContour;
	vector<ContourSector> 	m_minimalConstantSectors;
	vector<ContourSector> 	m_minimalInterpolatingSectors;
	set<SPoint>				m_uniquePoints;
	LSbox*					m_owningGrain;
	void orderCloudInContour();
	int getIDInMinimalisticContour(SPoint& toFind);
	int getConstantContourID(ContourSector* sector, vector<ContourSector>& constantSectors);
	int projectToMinimalBondary(SPoint point, double& out_lambda) const;
public:
	MinimalisticBoundary();
	~MinimalisticBoundary();
	void initialize(vector<SPoint>& originalContour,
					vector<ContourSector>& constantSectors,
					vector<ContourSector>& interpolatingSectors);
	double getWeight(int i, int j, LSbox* owner);
	vector<SPoint>&		getBoundary() {return m_minimalContour;}
};

#endif 	//__MINIMALISTIC_BOUNDARY__
