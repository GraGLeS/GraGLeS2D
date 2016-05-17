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
#include "junction.h"
#include "box.h"

grainhdl* GrainJunction::handler = NULL;

double GrainJunction::getWeight(LSbox* me) {
	unsigned int neighbours[3] = { 0xFFFFFF, 0xFFFFFF, 0xFFFFFF };
	for (int i = 0, j = 0; i < junction_type; i++) {
		if (me->getID() != grains[i]) {
			neighbours[j] = grains[i];
			j++;
		}
	}
	if (!handler->project == E_TRIPLE_JUNCTION_DRAG_SINGLE) {

		if (neighbours[0] == 0xFFFFFF || neighbours[1] == 0xFFFFFF) {
			return 1;
		}
	}
	double averageMobility = 0;
	double sigma;
	double gamma[3] = { 0.0, 0.0, 0.0 };
	double gamma_hagb = 1.0;
	double theta_ref = 15.0 * PI / 180;
	double theta_mis;

	if (Settings::MicrostructureGenMode == E_GENERATE_WITH_VORONOY
			|| Settings::MicrostructureGenMode == E_READ_FROM_FILE
			|| Settings::MicrostructureGenMode
					== E_READ_VOXELIZED_MICROSTRUCTURE) {

		if (Settings::ResearchMode == 0) {

			theta_mis = me->computeMisorientation(neighbours[0]);
			//characteristics& myNeighbor = me->m_grainBoundary.getDirectNeighbourCaracteristic(handler->getGrainByID(neighbours[0]));
			//averageMobility += myNeighbor.mobility;
			averageMobility += me->GBmobilityModel(theta_mis,
					handler->getGrainByID(neighbours[0]));
			gamma[0] = me->GBEnergyReadShockley(theta_mis,
					handler->getGrainByID(neighbours[0]));
			if (Settings::IdentifyTwins)
				if (me->MisoriToTwinBoundary(
						handler->getGrainByID(neighbours[0])) < 8.66025 * PI
						/ 180.0) {
					//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
					gamma[0] = me->GBEnergyReadShockley(1 * PI / 180,
							handler->getGrainByID(neighbours[0]));
				}

			theta_mis
					= handler->getGrainByID(neighbours[0])->computeMisorientation(
							neighbours[1]);
			characteristics
					& myNeighbor =
							handler->getGrainByID(neighbours[0])->m_grainBoundary.getDirectNeighbourCaracteristic(
									handler->getGrainByID(neighbours[1]));
			averageMobility += myNeighbor.mobility;
			gamma[1]
					= handler->getGrainByID(neighbours[0])->GBEnergyReadShockley(
							theta_mis, handler->getGrainByID(neighbours[1]));
			if (Settings::IdentifyTwins)
				if (handler->getGrainByID(neighbours[0])->MisoriToTwinBoundary(
						handler->getGrainByID(neighbours[1])) < 8.66025 * PI
						/ 180.0) {
					//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
					gamma[1]
							= handler->getGrainByID(neighbours[0])->GBEnergyReadShockley(
									1 * PI / 180,
									handler->getGrainByID(neighbours[1]));
				}

			//myNeighbor = me->m_grainBoundary.getDirectNeighbourCaracteristic(handler->getGrainByID(neighbours[1]));
			//averageMobility += myNeighbor.mobility;

			theta_mis = me->computeMisorientation(neighbours[1]);
			averageMobility += me->GBmobilityModel(theta_mis,
					handler->getGrainByID(neighbours[1]));
			gamma[2] = me->GBEnergyReadShockley(theta_mis,
					handler->getGrainByID(neighbours[1]));
			if (Settings::IdentifyTwins)
				if (me->MisoriToTwinBoundary(
						handler->getGrainByID(neighbours[1])) < 8.66025 * PI
						/ 180.0) {
					//if the GB is TWIN Boundary set the Energy to the lowest one allowed : theta_mis = 1 degree
					gamma[2] = me->GBEnergyReadShockley(1 * PI / 180,
							handler->getGrainByID(neighbours[1]));
				}

		} else if (Settings::ResearchMode == 1) {
			if (handler->project == E_TRIPLE_JUNCTION_DRAG_SINGLE) {
				if (coordinates.x < handler->delta + 2
						+ handler->get_grid_blowup() || coordinates.x
						> handler->get_ngridpoints()
								- handler->get_grid_blowup() - 2
						|| coordinates.y < handler->delta + 2
								+ handler->get_grid_blowup() || coordinates.y
						> handler->get_ngridpoints()
								- handler->get_grid_blowup() - 2) {
					//return 1.0;
					//cout << "Corner in 4er" << endl;

					return 0.0;

				}
			} else {
				gamma[0] = 1.0;
				gamma[1] = 1.0;
				gamma[2] = 1.0;
			}
		} else if (Settings::ResearchMode == 2) {
			{

				//! Generates for one half of the GBs a high angle
				//! and for the other half a low angle GB energy
				theta_ref = 42 * PI / 180;
				theta_mis = me->computeMisorientation(neighbours[0]);
				if (theta_mis <= theta_ref)
					gamma[0] = 0.3;
				else
					gamma[0] = gamma_hagb;

				theta_mis
						= handler->getGrainByID(neighbours[0])->computeMisorientation(
								neighbours[1]);
				if (theta_mis <= theta_ref)
					gamma[1] = 0.3;
				else
					gamma[1] = gamma_hagb;

				theta_mis = me->computeMisorientation(neighbours[1]);
				if (theta_mis <= theta_ref)
					gamma[2] = 0.3;
				else
					gamma[2] = gamma_hagb;

			}
		}

	}
	//! Distributes weights on edges for E_GENERATE_TESTCASE
	else if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		gamma[0] = handler->weightsMatrix[me->getID()][neighbours[0]];
		gamma[1] = handler->weightsMatrix[neighbours[0]][neighbours[1]];
		gamma[2] = handler->weightsMatrix[me->getID()][neighbours[1]];

		//		cout << gamma[0] << " ++ " << gamma[1] << " ++ " << gamma[2] << endl;
	} else {
		cout << "error in junction.cpp - no rule found" << endl;
		exit(2);
	}

	// find the asociated weight
	sigma = gamma[0] - gamma[1] + gamma[2];
	if (Settings::IsIsotropicNetwork == 0) {
		if (Settings::UseMobilityModel > 0 && Settings::TriplePointDrag > 0) {
			averageMobility /= 3;
			double ds = handler->get_ds();
			double drag = 1 / ((1 / (ds * Settings::TriplePointDrag)) + 1
					/ averageMobility);
			sigma *= drag;
		} else if (handler->project == E_TRIPLE_JUNCTION_DRAG_SINGLE
				&& handler->loadCurvature) {
			if (handler->loop > handler->loadCurvatureLoop)
				sigma *= Settings::TriplePointDrag;
			//else "no drag factor mapped"
		}
	}

	if (sigma < 0.01) {
		//cout << "negative sigma " << endl;
		sigma = 0.01;
	}
	return sigma;

}
