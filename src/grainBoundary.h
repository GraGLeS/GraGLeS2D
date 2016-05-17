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

#ifndef		__GRAIN_BOUNDARY__
#define		__GRAIN_BOUNDARY__

#include "dimensionalBufferReal.h"
#include "dimensionalBufferIDLocal.h"
#include "contourSector.h"
#include "junction.h"
#include "spoint.h"
#include "charasteristic.h"

#include <vector>
#include <map>
using namespace std;

class LSbox;

/*!
 * \class ExplicitGrainBoundary
 * \brief Class that holds all information required to describe a grain boundary.
 *
 * This class is a utility class that is used to encapsulate the explicit artifacts of a
 * grain boundary. It is responsible for generating the explicit boundary from a given
 * distance function, identify the neighboring grains to the given grain, construct all sectors
 * on the boundary as well as compute information for the grain like perimeter and volume.
 */

class ExplicitGrainBoundary
{
public:
	ExplicitGrainBoundary(LSbox* owner);
	~ExplicitGrainBoundary();

	bool extractBoundaryAndJunctions(DimensionalBufferReal& distanceBuffer, DimensionalBufferIDLocal& idLocal, int timestep=-1, bool verbose = false);
	void buildBoundarySectors(DimensionalBufferIDLocal& idLocal, int timestep = -1, bool verbose = false);
	void buildDirectNeighbours(DimensionalBufferReal& distanceBuffer, DimensionalBufferIDLocal& idLocal,  int timestep = -1, bool verbose = false);
	double computeVolume(int timestep = -1, bool verbose = false);
	double computeEnergy(int timestep = -1, bool verbose = false);
	double computePerimeter(int timestep = -1, bool verbose = false);
	double getWeight(int i, int j, DimensionalBufferIDLocal& idLocal, int timestep = -1, bool verbose = false);
	double getMeanM() const;
	double getMeanA() const;
	inline vector<SPoint>&			getRawBoundary() { return m_grainBoundary; }
	inline vector<GrainJunction>&	getRawJunctions() { return m_grainJunctions; }
	vector<int>						getDirectNeighbours();
	vector<double>              	getGbSegmentLength();
	inline int	getBoundarySegmentCount() { return m_grainBoundary.size(); }
	inline int  getDirectNeighboursCount() { return m_directNeighbourhood.size(); }
	inline void addDirectNeighbourManual( LSbox* neighbour ) { m_directNeighbourhood.push_back(characteristics(neighbour, 0, 0, 0,0)); }
	characteristics& getDirectNeighbourCaracteristic(LSbox* neighbour);
	map<int, double>& getlocalMODF();
	bool isBoxDirectNeighbour(LSbox* neighbour);
	SPoint calculateCentroid();
	double get_f_StoredElasticEnergy(int neighbor);
	double get_f_magneticEnergy(int neighbor);
	double getMobility(int neighbor);
	SPoint get_GB_Element(int i){return m_grainBoundary[i];};
	double DistanceToGrainBondary(SPoint point) const;
private:

	int projectToGrainBondary(SPoint point, double& out_lambda) const;

	void setPointsOnBoundary(ContourSector& sector);
	void calculateDiscreteEnergyDistribution();

	LSbox*					m_owningGrain;
	vector<ContourSector>	m_constantSectors;
	vector<ContourSector>	m_interpolatingSectors;
	vector<characteristics>	m_directNeighbourhood;
	vector<SPoint>			m_grainBoundary;
	vector<GrainJunction>	m_grainJunctions;
	map<int,double>			m_localMODF;
	int						m_grainBoundaryTimestep;
	int						m_sectorsTimestep;
	int						m_directNeighbourhoodTimestep;
};

#endif		//__GRAIN_BOUNDARY__
