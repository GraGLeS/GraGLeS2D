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
#ifndef		__ITERATIVE_GRAIN_SCHEDULER__
#define		__ITERATIVE_GRAIN_SCHEDULER__
#include "IGrainScheduler.h"
using namespace std;

class IterativeGrainScheduler : public IGrainScheduler
{
public:
	IterativeGrainScheduler(int numberOfThreads, int totalNumberOfGrains, int startIndex = 1);
	~IterativeGrainScheduler();

	void buildGrainWorkloads(vector<vector<SPoint>>& contours, int n_gridpoints);
	vector<unsigned int>&	getThreadWorkload(int threadID);

private:
	int								m_totalNumberOfThreads;
	int 							m_totalNumberOfGrains;
	int 							m_startIndex;
	vector<vector<unsigned int> >	m_threadWorkVectors;
};

#endif		//__ITERATIVE_GRAIN_SCHEDULER__




