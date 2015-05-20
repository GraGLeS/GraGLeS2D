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
#include "minimalisticBoundary.h"
#include "grahamScan.h"
#include "box.h"
#include <algorithm>
MinimalisticBoundary::MinimalisticBoundary()
{

}

MinimalisticBoundary::~MinimalisticBoundary()
{

}

void MinimalisticBoundary::orderCloudInContour()
{
	struct ComparisonFunctor{
		const SPoint*	center;
		bool operator() (const SPoint& lhs, const SPoint& rhs)
		{
			double val = (rhs.x - center->x)*(lhs.y - center->y) -
						 (rhs.y - center->y)*(lhs.x - center->x);
			if(val == 0)
			{
				return center->squaredDistanceTo(lhs) > center->squaredDistanceTo(rhs);
			}
			else
			{
				return val > 0;
			}
		}
	};

	SPoint massCenter(0,0,0);
	for(unsigned int i=0; i<m_minimalContour.size(); i++)
		massCenter = massCenter + m_minimalContour[i];
	massCenter = massCenter * (1.0 / m_minimalContour.size());

	ComparisonFunctor functor;
	functor.center = &massCenter;
	m_uniquePoints.clear();
	m_uniquePoints.insert(m_minimalContour.begin(), m_minimalContour.end());
	m_minimalContour.clear();
	copy(m_uniquePoints.begin(), m_uniquePoints.end(), back_inserter(m_minimalContour));
	sort(m_minimalContour.begin(), m_minimalContour.end(), functor);
	m_minimalContour.push_back(m_minimalContour[0]);

}

int MinimalisticBoundary::getIDInMinimalisticContour(SPoint& toFind)
{
	int id = -1;
	for(unsigned int i=0; i<m_minimalContour.size(); i++)
	{
		if(toFind == m_minimalContour[i])
		{
			id = i;
			break;
		}
	}
	return id;
}

int MinimalisticBoundary::getConstantContourID(ContourSector* sector, vector<ContourSector>& constantSectors)
{
	int id = -1;
	for(unsigned int i=0; i < constantSectors.size(); i++)
	{
		if (sector == &constantSectors[i])
		{
			id = i;
			break;
		}
	}
	return id;
}

void MinimalisticBoundary::initialize(	vector<SPoint>& originalContour,
										vector<ContourSector>& constantSectors,
										vector<ContourSector>& interpolatingSectors)
{
	m_minimalContour.clear();
	m_minimalConstantSectors.clear();
	m_minimalInterpolatingSectors.clear();
	for(unsigned int i=0; i<interpolatingSectors.size(); i++)
	{
		m_minimalContour.push_back(originalContour[interpolatingSectors[i].getLeftContourPoint()]);
		m_minimalContour.push_back(originalContour[interpolatingSectors[i].getRightContourPoint()]);
	}
	for(unsigned int i=0; i<constantSectors.size(); i++)
	{
		m_minimalContour.push_back(originalContour[constantSectors[i].getLeftContourPoint()]);
		m_minimalContour.push_back(originalContour[constantSectors[i].getRightContourPoint()]);

		for(unsigned int j=0; j<constantSectors[i].getJunctions().size(); j++)
		{
			m_minimalContour.push_back(constantSectors[i].getJunctions()[j]->coordinates);
		}
	}

	orderCloudInContour();

	for(unsigned int i=0; i<constantSectors.size(); i++)
	{
		ContourSector new_sector(E_CONSTANT_JUNCTION_SECTOR);
		for(unsigned int j = 0; j < constantSectors[i].getJunctions().size(); j++)
		{
			new_sector.mergeWith(constantSectors[i].getJunctions()[j]);
			new_sector.setLeftContourPoint(getIDInMinimalisticContour(
					originalContour[constantSectors[i].getLeftContourPoint()]));
			new_sector.setRightContourPoint(getIDInMinimalisticContour(
					originalContour[constantSectors[i].getRightContourPoint()]));
		}
		m_minimalConstantSectors.push_back(new_sector);
	}
	for(unsigned int i=0; i<interpolatingSectors.size(); i++)
	{
		ContourSector new_sector(E_INTERPOLATION_SECTOR);
		new_sector.setLeftContourPoint(getIDInMinimalisticContour(
					originalContour[interpolatingSectors[i].getLeftContourPoint()]));
		new_sector.setRightContourPoint(getIDInMinimalisticContour(
					originalContour[interpolatingSectors[i].getRightContourPoint()]));
		new_sector.setRightSectorBoundary( interpolatingSectors[i].getRightSectorBoundary() == NULL ?
		NULL : &m_minimalConstantSectors[
				getConstantContourID(interpolatingSectors[i].getRightSectorBoundary(),
										constantSectors)]);
		new_sector.setLeftSectorBoundary( interpolatingSectors[i].getLeftSectorBoundary() == NULL ?
		NULL : &m_minimalConstantSectors[
				getConstantContourID(interpolatingSectors[i].getLeftSectorBoundary(),
										constantSectors)]);
		new_sector.recalculateGrainBoundaryLength(m_minimalContour, m_owningGrain->get_h());
		m_minimalInterpolatingSectors.push_back(new_sector);
	}
}

int MinimalisticBoundary::projectToMinimalBondary(SPoint point, double& out_lambda) const {
	double minDist = 1000000.0;
	int minimal_segment =0;
	for (unsigned int k = 1, l = 0; k < m_minimalContour.size(); k++, l++) {
		SPoint u = m_minimalContour[k] - m_minimalContour[l];
		double lambda = (point - m_minimalContour[l]).dot(u);
		lambda /= u.dot(u);

		double dist;
		if (lambda < 0)
			dist = (point - m_minimalContour[l]).len();
		else if (lambda > 1)
			dist = (m_minimalContour[k] - point).len();
		else
			dist = (point - (m_minimalContour[l] + u * lambda)).len();

		if (dist < minDist) {
			minimal_segment = l;
			minDist = dist;
			out_lambda = lambda;
			if (out_lambda > 1)
				out_lambda = 1;
			if (out_lambda < 0)
				out_lambda = 0;
		}
	}
	return minimal_segment;
}

double MinimalisticBoundary::getWeight(int i, int j, LSbox* owner)
{
	double lambda = 0;
	int segment = projectToMinimalBondary(SPoint(j, i, 0), lambda);

	for (unsigned int q = 0; q < m_minimalConstantSectors.size(); q++)
		if (m_minimalConstantSectors[q].isSegmentWithinSector(m_minimalContour, segment)) {
			return m_minimalConstantSectors[q].getWeight(owner, &m_minimalContour,
					segment, lambda);

		}

	for (unsigned int q = 0; q < m_minimalInterpolatingSectors.size(); q++)
		if (m_minimalInterpolatingSectors[q].isSegmentWithinSector(m_minimalContour, segment)) {
			return m_minimalInterpolatingSectors[q].getWeight(owner, &m_minimalContour,
					segment, lambda);
		}

	return owner->getGBEnergyTimesGBMobility(i, j);
}
