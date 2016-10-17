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

#include "grainBoundary.h"
#include "marchingSquares.h"
#include "Settings.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "Structs.h"

#define PERIODIC(x) (((x) + (m_grainBoundary.size()-1)) % (m_grainBoundary.size() -1 ))

ExplicitGrainBoundary::ExplicitGrainBoundary(LSbox* owner) :
		m_owningGrain(owner) {
}

ExplicitGrainBoundary::~ExplicitGrainBoundary() {
}

bool ExplicitGrainBoundary::extractBoundaryAndJunctions(
		DimensionalBufferReal& distanceBuffer,
		DimensionalBufferIDLocal& idLocal, int timestep, bool verbose) {
	m_grainBoundary.clear();
	m_grainJunctions.clear();

	m_grainBoundaryTimestep = timestep;

	MarchingSquaresAlgorithm marchingSquares(distanceBuffer, idLocal,
			m_owningGrain);
	bool has_contour = marchingSquares.generateContour(m_grainBoundary,
			m_grainJunctions);
	if (false == has_contour)
		return false;

	//Project all found grain junctions to the grain boundary
	double lambda;
	for (unsigned int i = 0; i < m_grainJunctions.size(); i++) {
		int segment = projectToGrainBondary(m_grainJunctions[i].coordinates,
				lambda);
		m_grainJunctions[i].coordinates = m_grainBoundary[segment]
				+ (m_grainBoundary[PERIODIC(segment + 1)]
						- m_grainBoundary[segment]) * lambda;
		m_grainJunctions[i].contourSegment = segment;
		m_grainJunctions[i].lambda = lambda;
	}

	if (verbose) {
		stringstream filename;
		filename << "ContourLine_" << m_owningGrain->getID() << "_Timestep_"
				<< timestep << ".gnu";
		ofstream file;
		file.open(filename.str());
		for (unsigned int i = 0; i < m_grainBoundary.size(); i++)
			file << m_grainBoundary[i].x << " " << m_grainBoundary[i].y << "\n";
		file.close();

		filename.str("");
		filename << "Junctions_" << m_owningGrain->getID() << "_Timestep_"
				<< timestep << ".gnu";
		file.open(filename.str());
		for (unsigned int i = 0; i < m_grainJunctions.size(); i++)
			file << m_grainJunctions[i].coordinates.x << " "
					<< m_grainJunctions[i].coordinates.y << endl;
		file.close();
	}
	return true;
}

void ExplicitGrainBoundary::buildDirectNeighbours(
		DimensionalBufferReal& distanceBuffer,
		DimensionalBufferIDLocal& idLocal, int timestep, bool verbose) {
	double line_length;
	double h = m_owningGrain->get_h();
	double thetaMis = 0;
	vector<characteristics>::iterator it;

	m_directNeighbourhoodTimestep = timestep;
	m_directNeighbourhood.clear();

	for (unsigned int i = 0; i < m_grainBoundary.size() - 1; i++) {
		double px = m_grainBoundary[i].x;
		double py = m_grainBoundary[i].y;
		int pxGrid = int(px);
		int pyGrid = int(py);
		LSbox* candidate = NULL;
		// evaluate the ID at a outer point -> there must be one, we are at the boundary
		try {
			if (distanceBuffer.getValueAt(pyGrid, pxGrid) > 0) {
				pxGrid = int(px + 1);
				if (distanceBuffer.getValueAt(pyGrid, pxGrid) > 0) {
					pyGrid = int(py + 1);
					if (distanceBuffer.getValueAt(pyGrid, pxGrid) > 0) {
						pxGrid = int(px);
					}
				}
			}

			candidate = m_owningGrain->getNeighbourAt(pyGrid, pxGrid);
			if (candidate == NULL)
				throw std::out_of_range("candidate neighbour is null");
		} catch (const std::out_of_range& oor) {
			m_owningGrain->markAsInvalidMotion();
			return;
		}
		//		if (Settings::IsIsotropicNetwork == 0 && (Settings::ResearchMode != 1)) {
		//			thetaMis = m_owningGrain->computeMisorientation(candidate);
		//			m_grainBoundary[i].mob = m_owningGrain->GBmobilityModel(thetaMis,
		//					candidate);
		//			if (thetaMis <= theta_ref) {
		//				if (thetaMis < 1 * PI / 180)
		//					thetaMis = 1 * PI / 180;
		//			}
		//			m_grainBoundary[i].energy = m_owningGrain->GBEnergyReadShockley(
		//					thetaMis);
		//			if (Settings::IdentifyTwins)
		//				if (m_owningGrain->MisoriToTwinBoundary(candidate) < 8.66025
		//						* PI / 180.0) {
		//					//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
		//					m_grainBoundary[i].energy
		//							= m_owningGrain->GBEnergyReadShockley(1 * PI / 180);
		//				}
		//			//check for nan:
		//			if (m_grainBoundary[i].energy != m_grainBoundary[i].energy) {
		//				cout << m_owningGrain->getID() << "nan in Energy computation "
		//						<< endl;
		//				m_grainBoundary[i].energy = 0.01;
		//			}
		//
		//		} else if (Settings::ResearchMode == 1) {
		//			if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		//				m_grainBoundary[i].energy
		//						= m_owningGrain->getWeigthFromHandler(
		//								m_owningGrain->getID(), candidate->getID());
		//				m_grainBoundary[i].mob = 1.0;
		//			}
		//
		//			else if (Settings::ResearchProject != 0) {
		//				m_grainBoundary[i].energy = 1.0;
		//				m_grainBoundary[i].mob = 1.0;
		//			} else {
		//				thetaMis = m_owningGrain->computeMisorientation(candidate);
		//				theta_ref = 42 * PI / 180;
		//				if (thetaMis <= theta_ref)
		//					m_grainBoundary[i].energy = 0.3;
		//				else
		//					m_grainBoundary[i].energy = gamma_hagb;
		//				m_grainBoundary[i].mob = 1.0;
		//			}
		//		} else if (Settings::IsIsotropicNetwork) {
		//			m_grainBoundary[i].energy = 1.0;
		//			m_grainBoundary[i].mob = 1.0;
		//		}
		thetaMis = m_owningGrain->computeMisorientation(candidate);
		double StoredElasticEnergy = candidate->get_StoredElasticEnergy();
		double magneticEnergy = candidate->get_magneticEnergy();
		m_grainBoundary[i].mob = m_owningGrain->GBmobilityModel(thetaMis,
				candidate);
		m_grainBoundary[i].energy = m_owningGrain->GBEnergyReadShockley(
				thetaMis, candidate);
		line_length = (m_grainBoundary[i] - m_grainBoundary[i + 1]).len();
		line_length *= h;

		//  save length in GrainCharaczeristics

		for (it = m_directNeighbourhood.begin();
				it != m_directNeighbourhood.end(); it++) {
			if (candidate == (*it).directNeighbour)
				break;
		}
		if (it == m_directNeighbourhood.end()) {
			m_directNeighbourhood.push_back(
					characteristics(candidate, 0, m_grainBoundary[i].energy,
							thetaMis, m_grainBoundary[i].mob,
							StoredElasticEnergy, magneticEnergy));
			it = m_directNeighbourhood.end();
			it--;
		}
		it->length += line_length;
	}
	m_grainBoundary[m_grainBoundary.size() - 1].energy =
			m_grainBoundary[0].energy;
	m_grainBoundary[m_grainBoundary.size() - 1].mob = m_grainBoundary[0].mob;
}

double ExplicitGrainBoundary::computeVolume(int timestep, bool verbose) {
	double current_volume = 0;
	double h = m_owningGrain->get_h();

	for (unsigned int i = 0; i < m_grainBoundary.size() - 1; i++) {
		// Gaussian Trapez Formula:
		current_volume += (m_grainBoundary[i].y + m_grainBoundary[i + 1].y)
				* (m_grainBoundary[i].x - m_grainBoundary[i + 1].x);
	}

	current_volume = abs(current_volume) * h * h * 0.5;
	return current_volume;
}

double ExplicitGrainBoundary::computeEnergy(int timestep, bool verbose) {
	double total_energy = 0;
	for (unsigned int i = 0; i < m_grainBoundary.size() - 1; i++) {
		total_energy += (m_grainBoundary[i] - m_grainBoundary[i + 1]).len()
				* m_grainBoundary[i].energy * m_owningGrain->get_h();
	}
	return total_energy;
}

double ExplicitGrainBoundary::computePerimeter(int timestep, bool verbose) {
	double perimeter = 0;
	for (unsigned int i = 0; i < m_grainBoundary.size() - 1; i++) {
		perimeter += (m_grainBoundary[i] - m_grainBoundary[i + 1]).len()
				* m_owningGrain->get_h();
	}
	return perimeter;
}

void ExplicitGrainBoundary::buildBoundarySectors(
		DimensionalBufferIDLocal& idLocal, int timestep, bool verbose) {
	m_constantSectors.clear();
	m_interpolatingSectors.clear();

	if (m_grainJunctions.size() == 0)
		return;

	m_constantSectors.push_back(ContourSector(&m_grainJunctions[0]));
	setPointsOnBoundary(m_constantSectors[0]);
	int currentlyBuilt = 0;
	for (unsigned int i = 1; i < m_grainJunctions.size(); i++) {
		ContourSector currentSector(&m_grainJunctions[i]);
		setPointsOnBoundary(currentSector);
		if (m_constantSectors[currentlyBuilt].isSegmentWithinSector(
				m_grainBoundary, currentSector.getRightContourPoint())) {
			m_constantSectors[currentlyBuilt].mergeWith(&m_grainJunctions[i]);
			m_constantSectors[currentlyBuilt].setLeftContourPoint(
					currentSector.getLeftContourPoint());
			if (i == m_grainJunctions.size() - 1) {
				if (m_constantSectors[0].isSegmentWithinSector(m_grainBoundary,
						m_constantSectors[currentlyBuilt].getLeftContourPoint())) {
					if (currentlyBuilt == 0) {
						m_constantSectors[0].setLeftContourPoint(-1);
						m_constantSectors[0].setRightContourPoint(-1);
					} else {
						for (unsigned int j = 0;
								j
										< m_constantSectors[currentlyBuilt].getJunctions().size();
								j++) {
							m_constantSectors[0].mergeWith(
									m_constantSectors[currentlyBuilt].getJunctions()[j]);
						}
						m_constantSectors[0].setRightContourPoint(
								m_constantSectors[currentlyBuilt].getRightContourPoint());
					}
				}
			}
		} else {
			if (i == m_grainJunctions.size() - 1) {
				if (m_constantSectors[0].isSegmentWithinSector(m_grainBoundary,
						currentSector.getLeftContourPoint())) {
					m_constantSectors[0].setRightContourPoint(
							currentSector.getRightContourPoint());
					m_constantSectors[0].mergeWith(&m_grainJunctions[i]);
				} else {
					m_constantSectors.push_back(currentSector);
					currentlyBuilt++;
				}
			} else {
				m_constantSectors.push_back(currentSector);
				currentlyBuilt++;
			}
		}
	}

	int sec_len = m_constantSectors.size();
	m_interpolatingSectors.clear();
	//If there is only one sector then grain is completely inside one constant region
	//There is no need for interpolating sectors
	if (m_constantSectors.size() > 1
			|| m_constantSectors[0].getLeftContourPoint()
					!= m_constantSectors[0].getRightContourPoint()) {
		for (unsigned int i = 0;
				i < m_constantSectors.size() && m_constantSectors.size() > 1;
				i++) {
			ContourSector& secR = m_constantSectors[i];
			ContourSector& secL = m_constantSectors[(i + 1 + sec_len) % sec_len];

			ContourSector new_right = secR.generateLeftInterpolatingContour(
					m_grainBoundary, m_owningGrain->get_h());
			ContourSector new_left = secL.generateRightInterpolatingContour(
					m_grainBoundary, m_owningGrain->get_h());
			bool merged = false;

			for (int j = new_right.getRightContourPoint();
					j != PERIODIC(new_right.getLeftContourPoint() + 1); j =
							PERIODIC(j + 1)) {
				if (j == new_left.getLeftContourPoint()
						|| j == new_left.getRightContourPoint()) {
					merged = true;
				}
			}
			if (merged) {
				m_interpolatingSectors.push_back(
						secL.generateInterpolatingSector(secR, m_grainBoundary,
								m_owningGrain->get_h()));
			} else {
				m_interpolatingSectors.push_back(new_right);
				m_interpolatingSectors.push_back(new_left);
			}
		}
	}

	if (verbose) {
		ofstream outputFile;
		stringstream outputFilename;
		outputFilename << "Sectors_" << m_owningGrain->getID() << "_Timestep_"
				<< timestep << ".gnu";
		outputFile.open(outputFilename.str());
		for (unsigned int i = 0; i < m_constantSectors.size(); i++) {
			m_constantSectors[i].debugPrintSector(m_grainBoundary, outputFile);
		}
		for (unsigned int i = 0; i < m_interpolatingSectors.size(); i++) {
			m_interpolatingSectors[i].debugPrintSector(m_grainBoundary,
					outputFile);
		}
		outputFile.close();
	}
}

double ExplicitGrainBoundary::getWeight(int i, int j,
		DimensionalBufferIDLocal& idLocal, int timestep, bool verbose) {
	double lambda;
	int segment = projectToGrainBondary(SPoint(j, i, 0, 0), lambda);

	for (unsigned int q = 0; q < m_constantSectors.size(); q++)
		if (m_constantSectors[q].isSegmentWithinSector(m_grainBoundary,
				segment)) {
			return m_constantSectors[q].getWeight(m_owningGrain,
					&m_grainBoundary, segment, lambda);
		}

	for (unsigned int q = 0; q < m_interpolatingSectors.size(); q++)
		if (m_interpolatingSectors[q].isSegmentWithinSector(m_grainBoundary,
				segment)) {
			return m_interpolatingSectors[q].getWeight(m_owningGrain,
					&m_grainBoundary, segment, lambda);
		}
	return m_owningGrain->getGBEnergyTimesGBMobility(i, j);
}

int ExplicitGrainBoundary::projectToGrainBondary(SPoint point,
		double& out_lambda) const {
	double minDist = 1000000.0;
	int minimal_segment = 0;
	// first search on vertices:
	double dist = 0;
	int asize = m_grainBoundary.size();
	int nearestPoint = -1;
	for (int k = 0; k < asize - 1; k++) { //check p_1 .. p_N, skip last entry which is again p_1
		SPoint u = m_grainBoundary[k];
		dist = (point - u).lenSqr();
		if (dist < minDist) {
			nearestPoint = k;
			minDist = dist;
			minimal_segment = k;
		}
	}
	// start local search on adjacent line segments - here the order plays an important role, as the vertex is the next point on the inface to point
	minDist = 1000000.0;
	int next = nearestPoint + 1;
	SPoint u = m_grainBoundary[next] - m_grainBoundary[nearestPoint];
	double lambda = (point - m_grainBoundary[nearestPoint]).dot(u);
	lambda /= u.dot(u);
	if (lambda < 0)
		dist = (point - m_grainBoundary[nearestPoint]).lenSqr();
	else if (lambda > 1)
		dist = (m_grainBoundary[next] - point).lenSqr();
	else
		dist = (point - (m_grainBoundary[nearestPoint] + u * lambda)).lenSqr();
	if (dist < minDist) {
		minimal_segment = nearestPoint;
		minDist = dist;
		out_lambda = lambda;
	}
	// m_grainboundary {p_1,p_2,..,p_N, p_1}
	if (nearestPoint == 0)
		nearestPoint = asize - 1;

	int previous = nearestPoint - 1;
	u = m_grainBoundary[nearestPoint] - m_grainBoundary[previous];
	lambda = (point - m_grainBoundary[previous]).dot(u);
	lambda /= u.dot(u);
	if (lambda < 0)
		dist = (point - m_grainBoundary[previous]).lenSqr();
	else if (lambda > 1)
		dist = (m_grainBoundary[nearestPoint] - point).lenSqr();
	else
		dist = (point - (m_grainBoundary[previous] + u * lambda)).lenSqr();
	if (dist < minDist) {
		minimal_segment = previous;
		minDist = dist;
		out_lambda = lambda;
	}

	if (out_lambda > 1)
		out_lambda = 1;
	if (out_lambda < 0)
		out_lambda = 0;

	// original algorithm:
	//	for (unsigned int k = 1, l = 0; k < m_grainBoundary.size(); k++, l++) {
	//		SPoint u = m_grainBoundary[k] - m_grainBoundary[l];
	//		double lambda = (point - m_grainBoundary[l]).dot(u);
	//		lambda /= u.dot(u);
	//		if (lambda < 0)
	//			dist = (point - m_grainBoundary[l]).lenSqr();
	//		else if (lambda > 1)
	//			dist = (m_grainBoundary[k] - point).lenSqr();
	//		else
	//			dist = (point - (m_grainBoundary[l] + u * lambda)).lenSqr();
	//
	//		if (dist < minDist) {
	//			minimal_segment = l;
	//			minDist = dist;
	//			out_lambda = lambda;
	//			if (out_lambda > 1)
	//				out_lambda = 1;
	//			if (out_lambda < 0)
	//				out_lambda = 0;
	//		}
	//	}
	return minimal_segment;
}

double ExplicitGrainBoundary::DistanceToGrainBondary(SPoint point) const {
	double minDist = 1000000.0;
	int minimal_segment = 0;
	// first search on vertices:
	double dist = 0;
	int asize = m_grainBoundary.size();
	int nearestPoint = -1;
	for (int k = 0; k < asize - 1; k++) { //check p_1 .. p_N, skip last entry which is again p_1
		SPoint u = m_grainBoundary[k];
		dist = (point - u).lenSqr();
		if (dist < minDist) {
			nearestPoint = k;
			minDist = dist;
			minimal_segment = k;
		}
	}
	// start local search on adjacent line segments - here the order plays an important role, as the vertex is the next point on the inface to point
	minDist = 1000000.0;
	int next = nearestPoint + 1;
	SPoint u = m_grainBoundary[next] - m_grainBoundary[nearestPoint];
	double lambda = (point - m_grainBoundary[nearestPoint]).dot(u);
	lambda /= u.dot(u);
	if (lambda < 0)
		dist = (point - m_grainBoundary[nearestPoint]).lenSqr();
	else if (lambda > 1)
		dist = (m_grainBoundary[next] - point).lenSqr();
	else
		dist = (point - (m_grainBoundary[nearestPoint] + u * lambda)).lenSqr();
	if (dist < minDist) {
		minimal_segment = nearestPoint;
		minDist = dist;
	}
	// m_grainboundary {p_1,p_2,..,p_N, p_1}
	if (nearestPoint == 0)
		nearestPoint = asize - 1;

	int previous = nearestPoint - 1;
	u = m_grainBoundary[nearestPoint] - m_grainBoundary[previous];
	lambda = (point - m_grainBoundary[previous]).dot(u);
	lambda /= u.dot(u);
	if (lambda < 0)
		dist = (point - m_grainBoundary[previous]).lenSqr();
	else if (lambda > 1)
		dist = (m_grainBoundary[nearestPoint] - point).lenSqr();
	else
		dist = (point - (m_grainBoundary[previous] + u * lambda)).lenSqr();
	if (dist < minDist) {
		minimal_segment = previous;
		minDist = dist;
	}
	// original algorithm:
	//	for (unsigned int k = 1, l = 0; k < m_grainBoundary.size(); k++, l++) {
	//		SPoint u = m_grainBoundary[k] - m_grainBoundary[l];
	//		double lambda = (point - m_grainBoundary[l]).dot(u);
	//		lambda /= u.dot(u);
	//		if (lambda < 0)
	//			dist = (point - m_grainBoundary[l]).lenSqr();
	//		else if (lambda > 1)
	//			dist = (m_grainBoundary[k] - point).lenSqr();
	//		else
	//			dist = (point - (m_grainBoundary[l] + u * lambda)).lenSqr();
	//
	//		if (dist < minDist) {
	//			minimal_segment = l;
	//			minDist = dist;
	//			out_lambda = lambda;
	//			if (out_lambda > 1)
	//				out_lambda = 1;
	//			if (out_lambda < 0)
	//				out_lambda = 0;
	//		}
	//	}
	return minDist;
}

void ExplicitGrainBoundary::setPointsOnBoundary(ContourSector& sector) {
	double dist = 0;
	//Construct left side
	int l_segment = sector.getJunctions()[0]->contourSegment;
	//dist += ((m_grainBoundary[PERIODIC(l_segment + 1)]	- m_grainBoundary[l_segment]) * (1	- sector.getJunctions()[0]->lambda)).len();
	dist += ((m_grainBoundary[PERIODIC(l_segment + 1)]
			- m_grainBoundary[l_segment])
			* (1 - sector.getJunctions()[0]->lambda)).len()
			* m_owningGrain->get_h();
	l_segment = PERIODIC(l_segment + 1);

	while (dist < Settings::ConstantSectorRadius) {
		//dist += (m_grainBoundary[PERIODIC(l_segment + 1)] - m_grainBoundary[l_segment]).len();
		dist += (m_grainBoundary[PERIODIC(l_segment + 1)]
				- m_grainBoundary[l_segment]).len() * m_owningGrain->get_h();
		l_segment = PERIODIC(l_segment + 1);
	}
	dist = 0;
	//Construct right side
	int r_segment = sector.getJunctions()[0]->contourSegment;
	//dist += ((m_grainBoundary[PERIODIC(r_segment + 1)]	- m_grainBoundary[r_segment]) * (sector.getJunctions()[0]->lambda)).len();
	dist +=
			((m_grainBoundary[PERIODIC(r_segment + 1)]
					- m_grainBoundary[r_segment])
					* (sector.getJunctions()[0]->lambda)).len()
					* m_owningGrain->get_h();
	r_segment = PERIODIC(r_segment);
	while (dist < Settings::ConstantSectorRadius) {
		//dist += (m_grainBoundary[PERIODIC(r_segment - 1)] - m_grainBoundary[r_segment]).len();
		dist += (m_grainBoundary[PERIODIC(r_segment - 1)]
				- m_grainBoundary[r_segment]).len() * m_owningGrain->get_h();
		r_segment = PERIODIC(r_segment - 1);
	}
	sector.setRightContourPoint(r_segment);
	sector.setLeftContourPoint(l_segment);
}

characteristics& ExplicitGrainBoundary::getDirectNeighbourCaracteristic(
		LSbox* neighbour) {
	static characteristics defaultCharacteristic(NULL, 0, 1.0, 1.0);
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		if (m_directNeighbourhood[i].directNeighbour == neighbour) {
			return m_directNeighbourhood[i];
		}
	}

	return defaultCharacteristic;
}

map<int, double>& ExplicitGrainBoundary::getlocalMODF() {
	calculateDiscreteEnergyDistribution();
	return m_localMODF;
}

void ExplicitGrainBoundary::calculateDiscreteEnergyDistribution() {
	m_localMODF.clear();
	double dh;
	if (Settings::LatticeType == 0)
		dh = 62.8 / (double) Settings::DiscreteSamplingRate;
	else
		dh = 92. / (double) Settings::DiscreteSamplingRate;

	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		int key = (int) ((m_directNeighbourhood[i].mis_ori * 180 / PI) / dh
				- 0.5);
		double length = m_directNeighbourhood[i].length;
		map<int, double>::iterator it = m_localMODF.find(key);
		if (it == m_localMODF.end()) {
			m_localMODF.insert(make_pair(key, length));
		} else {
			(*it).second += length;

		}
	}
}

bool ExplicitGrainBoundary::isBoxDirectNeighbour(LSbox* neighbour) {
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		if (m_directNeighbourhood[i].directNeighbour == neighbour)
			return true;
	}
	return false;
}
double ExplicitGrainBoundary::get_f_StoredElasticEnergy(int neighbor) {
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		if (m_directNeighbourhood[i].directNeighbour->getID() == neighbor)
			return (m_owningGrain->get_StoredElasticEnergy()
					- m_directNeighbourhood[i].StoredElasticEnergy)
					* m_directNeighbourhood[i].mobility;
	}
	return -1;
}
double ExplicitGrainBoundary::getMobility(int neighbor) {
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		if (m_directNeighbourhood[i].directNeighbour->getID() == neighbor)
			return m_directNeighbourhood[i].mobility;
	}
	return -1;
}
double ExplicitGrainBoundary::get_f_magneticEnergy(int neighbor) {
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		if (m_directNeighbourhood[i].directNeighbour->getID() == neighbor)
			return ((m_owningGrain->get_magneticEnergy()
					- m_directNeighbourhood[i].magneticEnergy)
					* m_directNeighbourhood[i].mobility);
	}
	return -1;
}

SPoint ExplicitGrainBoundary::calculateCentroid() {

	int nVertices = m_grainJunctions.size();
	SPoint cent;
	cent.x = 0.0;
	cent.y = 0.0;
	double areaSigned = 0.0;
	double x0 = 0.0; // First vertex' x value
	double y0 = 0.0; // Second vertex' Y value
	double x1 = 0.0; // Next vertex' x value
	double y1 = 0.0; // Next vertex' y value
	double loopArea = 0.0; // Intermediate area

	int i, j = 0;
	for (i = 0, j = nVertices - 1; i < nVertices; j = i++) {

		x0 = m_grainJunctions[i].coordinates.x;
		y0 = m_grainJunctions[i].coordinates.y;
		x1 = m_grainJunctions[j].coordinates.x;
		y1 = m_grainJunctions[j].coordinates.y;
		loopArea = x0 * y1 - x1 * y0;
		areaSigned += loopArea;
		cent.x += (x0 + x1) * loopArea;
		cent.y += (y0 + y1) * loopArea;
	}

	areaSigned *= 0.5;
	cent.x /= (6.0 * areaSigned);
	cent.y /= (6.0 * areaSigned);

	return cent;
}

double ExplicitGrainBoundary::getMeanM() const {
	double avg_m = 0;
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		avg_m +=
				m_directNeighbourhood[i].directNeighbour->getDirectNeighbourCount();
	}
	avg_m /= m_directNeighbourhood.size();
	return avg_m;

}

double ExplicitGrainBoundary::getMeanA() const {
	double avg_a = 0;
	for (unsigned int i = 0; i < m_grainJunctions.size(); i++) {
		if (i == m_grainJunctions.size() - 1) {
			avg_a += (m_grainJunctions[0].coordinates
					- m_grainJunctions[i].coordinates).len();
		} else {
			avg_a += (m_grainJunctions[i].coordinates
					- m_grainJunctions[i + 1].coordinates).len();
		}
	}
	avg_a *= m_owningGrain->get_h();
	avg_a /= m_grainJunctions.size();
	return avg_a;
}

vector<int> ExplicitGrainBoundary::getDirectNeighbours() {
	vector<int> neighbors;
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		neighbors.push_back(m_directNeighbourhood[i].directNeighbour->getID());
	}
	return neighbors;
}

vector<double> ExplicitGrainBoundary::getGbSegmentLength() {
	vector<double> SegmentLengths;
	for (unsigned int i = 0; i < m_directNeighbourhood.size(); i++) {
		SegmentLengths.push_back(m_directNeighbourhood[i].length);
	}
	return SegmentLengths;
}
vector<Face>* ExplicitGrainBoundary::getFaces() {
	vector<Face>* my_faces = new vector<Face>;
	for (vector<characteristics>::iterator it = m_directNeighbourhood.begin();
			it != m_directNeighbourhood.end(); it++) {
		if (it->directNeighbour->getID() < m_owningGrain->getID())
			my_faces->push_back(
					Face(it->length, m_owningGrain->getID(),
							it->directNeighbour->getID()));
	}
	return my_faces;
}

