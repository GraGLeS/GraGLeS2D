/*
 * SquaresGrainScheduler.cpp
 *
 *  Created on: Dec 15, 2015
 *      Author: cm654063
 */

#include "SquaresGrainScheduler.h"
#include "IterativeGrainScheduler.h"
#include "Settings.h"
#include <omp.h>
#include <iostream>

SquaresGrainScheduler::SquaresGrainScheduler(int numberOfThreads,
		int totalNumberOfGrains, int startIndex) :
	m_totalNumberOfThreads(numberOfThreads),
			m_totalNumberOfGrains(totalNumberOfGrains),
			m_startIndex(startIndex) {
	m_threadWorkVectors.resize(m_totalNumberOfThreads);
}
SquaresGrainScheduler::~SquaresGrainScheduler() {

}
void SquaresGrainScheduler::buildGrainWorkloads(
		vector<vector<SPoint>>& contours, int n_gridpoints) {
	if (Settings::ExecuteInParallel == 0 || omp_get_max_threads() == 1) {
		for (int i = 1; i < contours.size(); i++)
			m_threadWorkVectors.at(0).push_back(i);
	} else if (omp_get_max_threads() < 4) {
		for (int i = 1; i < contours.size(); i++)
			m_threadWorkVectors.at(0).push_back(i);
		std::cout<< "method not implemented for less than 4 threads, use sequential one. "<< endl;
	} else {
		vector<vector<unsigned int>> list;
		list.resize(4);
		for (unsigned int i = 1; i < contours.size(); i++) {
			int pos_x, pos_y;
			SPoint center = find_center(contours[i], n_gridpoints);
			//check for coordinate system
			if (center.x > 2) { // coordinates are in gridpoint ids
				pos_x = int(center.x / n_gridpoints + 0.5);
				pos_y = int(center.y / n_gridpoints + 0.5);
			} else { // coordinates are in between [0,1]
				pos_x = int(center.x + 0.5);
				pos_y = int(center.y + 0.5);
			}
			if (pos_x == 0 && pos_y == 0)
				list[0].push_back(i);
			else if (pos_x == 1 && pos_y == 0)
				list[1].push_back(i);
			else if (pos_x == 0 && pos_y == 1)
				list[2].push_back(i);
			else
				list[3].push_back(i);

		}
		//determine number of threads per square:
		int rest = omp_get_max_threads() % 4;
		int threads = omp_get_max_threads() / 4;
		int NumberOfThreads[4] = { 0, 0, 0, 0 };
		if (rest != 0) {
			for (int i = 0; i < 4; i++) {
				NumberOfThreads[i] = threads;
				if (rest > 0) {
					NumberOfThreads[i]++;
					rest--;
				}
			}
		} else {
			NumberOfThreads[0] = threads;
			NumberOfThreads[1] = threads;
			NumberOfThreads[2] = threads;
			NumberOfThreads[3] = threads;
		}

		/*
		 THERE ARE 4 NUMA NODES ON THE FIRST LEVEL
		 2--------3
		 |        |
		 |        |
		 |        |
		 0--------1
		 EACH OF WHICH HAS 32 CORES WHICH HOST THE FOLLOWING THREADS:
		 NODE 0 = THREADS { 0,...,31}
		 NODE 1 = THREADS {32,...,63}
		 NODE 2 = THREADS {64,...,95}
		 NODE 3 = THREADS {96,...,127}
		 THIS DISTRIBUTION IS GARAUNTEED BY A RUNTIME BINDIND OF THE THREADS WITH NUMA.H
		 */
		int offsetID = 0;
		for (int i = 0; i < list.size(); i++) {
			buildSubWorkloads(offsetID, list[i], NumberOfThreads[i]);
			offsetID += NumberOfThreads[i];
		}
	}
}

SPoint SquaresGrainScheduler::find_center(vector<SPoint>& contour,
		int n_gridpoints) {
	double x, y;
	double xmax = 0;
	double xmin = n_gridpoints;
	double ymax = 0;
	double ymin = xmin;
	for (int k = 0; k < (contour).size(); k++) {
		y = contour[k].y;
		x = contour[k].x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
	}
	SPoint center((xmax - xmin) / 2 + xmin, (ymax - ymin) / 2 + ymin, 0, 0);
	return center;
}

vector<unsigned int>& SquaresGrainScheduler::getThreadWorkload(int threadID) {
	return m_threadWorkVectors.at(threadID);
}
void SquaresGrainScheduler::buildSubWorkloads(int offsetThreadID,
		vector<unsigned int>& grainlist, int NumberOfThreads) {
	int NumberOfGrains = grainlist.size();
	for (int i = 0; i < NumberOfThreads; i++) {
		int threadID = offsetThreadID + i;
		m_threadWorkVectors.at(threadID).reserve(
				NumberOfGrains / NumberOfThreads + 1);
		for (int j = 0; j < NumberOfGrains / NumberOfThreads + 1; j++) {
			int pos = j * NumberOfThreads + i;
			if (pos >= NumberOfGrains)
				break;
			else {
				int grainID = grainlist[pos];
				m_threadWorkVectors.at(threadID).push_back(grainID);
			}

		}
	}
}

