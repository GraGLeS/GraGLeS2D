/*
 * SquaresGrainScheduler.h
 *
 *  Created on: Dec 14, 2015
 *      Author: cm654063
 */

#ifndef SQUARESGRAINSCHEDULER_H_
#define SQUARESGRAINSCHEDULER_H_

#include "IGrainScheduler.h"
using namespace std;

struct SPoint;

class SquaresGrainScheduler: public IGrainScheduler {
public:
	SquaresGrainScheduler(int numberOfThreads, int totalNumberOfGrains, int startIndex = 1);
	~SquaresGrainScheduler();

	void buildGrainWorkloads(vector<vector<SPoint>>& contours, int n_gridpoints);
	void buildSubWorkloads(int offsetThreadID, vector<unsigned int>& grainlist,
			int NumberOfThreads);
	SPoint find_center(vector<SPoint>& contour, int n_gridpoints);
	vector<unsigned int>& getThreadWorkload(int threadID);

private:
	int m_totalNumberOfThreads;
	int m_totalNumberOfGrains;
	int m_startIndex;
	vector<vector<unsigned int> > m_threadWorkVectors;
};

#endif /* SQUARESGRAINSCHEDULER_H_ */
