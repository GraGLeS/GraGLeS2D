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
#include "contourSector.h"
#include "grahamScan.h"
#include "box.h"
#include "spoint.h"
#include "Settings.h"

ContourSector::ContourSector(GrainJunction* initialJunction) :
	m_sectorType(E_CONSTANT_JUNCTION_SECTOR),
	m_leftContourPointID(-1),
	m_rightContourPointID(-1),
	m_leftSectorBoundary(NULL),
	m_rightSectorBoundary(NULL),
	m_grainBoundaryLength(0.0)
{
	m_junctions.push_back(initialJunction);
}

ContourSector::ContourSector(E_CONTOUR_TYPE type) :
		m_sectorType(type),
		m_leftContourPointID(-1),
		m_rightContourPointID(-1),
		m_leftSectorBoundary(NULL),
		m_rightSectorBoundary(NULL),
		m_grainBoundaryLength(0.0)
{
}

bool ContourSector::mergeWith(GrainJunction* junction)
{
	m_junctions.push_back(junction);
	return true;
}


void ContourSector::setLeftContourPoint(int ID)
{
	m_leftContourPointID = ID;
}
void ContourSector::setRightContourPoint(int ID)
{
	m_rightContourPointID = ID;
}

void ContourSector::setRightSectorBoundary(ContourSector* sector)
{
	m_rightSectorBoundary = sector;
}
void ContourSector::setLeftSectorBoundary(ContourSector* sector)
{
	m_leftSectorBoundary = sector;
}

void ContourSector::debugPrintSector(vector<SPoint>& contour_grain, ofstream& ofs)
{
	int realContourSize = contour_grain.size()-1;
	int left = m_leftContourPointID;
	int right = m_rightContourPointID;
	if(left == right && right == -1)
	{
		left = realContourSize;
		right = 0;
	}

	int width = right < left ? left - right + 1: left + realContourSize - right + 1;
	vector<SPoint> temp_pts;
	temp_pts.resize(2*width);
	for(int i=right, cnt=0; cnt< width; i = (i+1+realContourSize)%realContourSize, cnt++)
	{
		int prev_pnt = (i-1+realContourSize)%realContourSize;
		SPoint v1;
		v1.x = -(contour_grain[i] - contour_grain[prev_pnt]).y;
		v1.y = (contour_grain[i] - contour_grain[prev_pnt]).x;
		v1 = v1*(1/v1.len());
		temp_pts[cnt] = contour_grain[i] + v1*2.0;
		temp_pts[temp_pts.size()-cnt-1] = contour_grain[i] + v1*(-2.0);
	}

	for(unsigned int i=0; i<temp_pts.size(); i++)
	{
		ofs<<temp_pts[i].x<<" "<<temp_pts[i].y<<"\n";
	}
	ofs<<temp_pts[0].x<<" "<<temp_pts[0].y<<"\n";
	ofs<<"\n";
}

ContourSector ContourSector::generateLeftInterpolatingContour(vector<SPoint>& contour_grain, double h)
{
	ContourSector result;
	int left;
	double dist=0;
	int j =0;
	int realContourSize = contour_grain.size()-1;

	for(j = m_leftContourPointID; dist <= Settings::InterpolatingSectorRadius;
			j = (j+1+realContourSize)%realContourSize)
	{
		SPoint first = contour_grain[j];
		SPoint second = contour_grain[(j+1+realContourSize)%realContourSize];
		dist += (first - second).len() * h;
	}
	left =(j+1+realContourSize)%realContourSize;
	result.m_leftContourPointID = left;
	result.m_rightContourPointID = m_leftContourPointID;
	result.m_rightSectorBoundary = this;
	result.m_grainBoundaryLength = dist;
	return result;
}

ContourSector ContourSector::generateInterpolatingSector(ContourSector& other, vector<SPoint>& contour_grain, double h)
{
	ContourSector result;
	int realContourSize = contour_grain.size()-1;
	result.m_leftContourPointID = this->m_rightContourPointID;
	result.m_leftSectorBoundary = this;
	result.m_rightContourPointID = other.m_leftContourPointID;
	result.m_rightSectorBoundary = &other;
	double dist=0;
	for(int j = result.m_rightContourPointID; j != result.m_leftContourPointID;
				j = (j+1+realContourSize)%realContourSize)
	{
		SPoint first = contour_grain[j];
		SPoint second = contour_grain[(j-1+realContourSize)%realContourSize];
		dist += (first - second).len() * h;
	}
	result.m_grainBoundaryLength = dist;
	return result;
}

ContourSector ContourSector::generateRightInterpolatingContour(vector<SPoint>& contour_grain, double h)
{
	ContourSector result;
	int right;
	double dist=0;
	int j =0;
	int realContourSize = contour_grain.size()-1;

	for(j = m_rightContourPointID; dist <= Settings::InterpolatingSectorRadius;
			j = (j-1+realContourSize)%realContourSize)
	{
		SPoint first = contour_grain[j];
		SPoint second = contour_grain[(j-1+realContourSize)%realContourSize];
		dist += (first - second).len() * h;
	}
	right =(j-1+realContourSize)%realContourSize;
	result.m_leftContourPointID = m_rightContourPointID;
	result.m_rightContourPointID = right;
	result.m_leftSectorBoundary = this;
	result.m_grainBoundaryLength = dist;

	return result;
}

void ContourSector::recalculateGrainBoundaryLength(vector<SPoint>& contour_grain, double h)
{
	int realContourSize = contour_grain.size()-1;
	double dist=0;
	for(int j = m_rightContourPointID; j != m_leftContourPointID;
				j = (j+1+realContourSize)%realContourSize)
		{
			SPoint first = contour_grain[j];
			SPoint second = contour_grain[(j+1+realContourSize)%realContourSize];
			dist += (first - second).len() * h;
		}
	m_grainBoundaryLength = dist;
}

bool ContourSector::isSegmentWithinSector(vector<SPoint>& contour_grain, int segment) const
{
	int totalP = contour_grain.size()-1;
	if(m_rightContourPointID == m_leftContourPointID && m_leftContourPointID == -1)
		return true;
	if(m_rightContourPointID > m_leftContourPointID)
	{
		return ( segment >=0 && segment <= m_leftContourPointID ) ||
			   ( segment >= m_rightContourPointID && segment < totalP );
	}
	else
	{
		return (segment >= m_rightContourPointID && segment <= m_leftContourPointID);
	}
}

double ContourSector::getWeight(LSbox* owner, vector<SPoint>* contour_grain, int segment, double lambda) const
{
	if( m_sectorType == E_CONSTANT_JUNCTION_SECTOR )
	{
		double weight =0;
		for(unsigned int i=0; i<m_junctions.size(); i++)
		{
			weight += m_junctions[i]->getWeight(owner);
		}
		return weight/m_junctions.size();
	}
	else
	{
		if(contour_grain == NULL)
		{
			return -1.0;
		}

		int real_contour_size = contour_grain->size()-1;
		double left_w, right_w;
		if(m_leftSectorBoundary)
			left_w = m_leftSectorBoundary->getWeight(owner);
		else
		{
			int sample_x = (*contour_grain)[m_leftContourPointID].x;
			int sample_y = (*contour_grain)[m_leftContourPointID].y;
			left_w = owner->getGBEnergyTimesGBMobility(sample_y,sample_x);
		}
		if(m_rightSectorBoundary)
			right_w = m_rightSectorBoundary->getWeight(owner);
		else
		{
			int sample_x = (*contour_grain)[m_rightContourPointID].x;
			int sample_y = (*contour_grain)[m_rightContourPointID].y;
			right_w = owner->getGBEnergyTimesGBMobility(sample_y,sample_x);
		}
		double interpolator=0;
		int i=0;
		for(i = m_rightContourPointID; i != segment; i= (i+1 + real_contour_size)%real_contour_size)
		{
			interpolator += ((*contour_grain)[i+1]-(*contour_grain)[i]).len() * owner->get_h();
		}
		SPoint addition  = ((*contour_grain)[(i+1 + real_contour_size)%real_contour_size] - (*contour_grain)[i])*lambda;
		interpolator += addition.len();
		interpolator /= m_grainBoundaryLength;

		if(interpolator < 0)
			interpolator = 0;
		if(interpolator > 1)
			interpolator = 1;

		if(m_leftSectorBoundary && m_rightSectorBoundary)
		{
			double A,B,C;
			double amplitude;
			SPoint sample = (*contour_grain)[(m_rightContourPointID + 1) % real_contour_size];
			double boundaryWeigth = owner->getGBEnergyTimesGBMobility(sample.y, sample.x);
			double quot = (m_grainBoundaryLength > 2*Settings::InterpolatingSectorRadius ? 2*Settings::InterpolatingSectorRadius : m_grainBoundaryLength);
			amplitude = (boundaryWeigth - (left_w + right_w)/2)/(2*Settings::InterpolatingSectorRadius) * quot + (left_w + right_w)/2;
			//http://isezen.com/2012/01/15/quadratic-interpolation-three-point/
			double a0 = 2*right_w; double a1 = -4*amplitude; double a2 = 2*left_w;
			A = a0 + a1 + a2;
			B = -(a0*1.5 + a1 + a2*0.5);
			C = a0*0.5;
			double result = A*interpolator*interpolator + B*interpolator + C;
			return result;
		}
		else
			return interpolator*left_w + (1-interpolator)*right_w;

	}
	return 0;
}
