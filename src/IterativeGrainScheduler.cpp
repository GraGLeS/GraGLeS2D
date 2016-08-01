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
#include "IterativeGrainScheduler.h"


IterativeGrainScheduler::IterativeGrainScheduler(int numberOfThreads,
		int totalNumberOfGrains, int startIndex) :
	m_totalNumberOfThreads(numberOfThreads),
			m_totalNumberOfGrains(totalNumberOfGrains),
			m_startIndex(startIndex) {
	m_threadWorkVectors.resize(m_totalNumberOfThreads);
}
IterativeGrainScheduler::~IterativeGrainScheduler(){
}
void IterativeGrainScheduler::buildGrainWorkloads(vector<vector<SPoint>>& contours, int n_gridpoints) {
	for (int i = 0; i < m_totalNumberOfThreads; i++) {
		m_threadWorkVectors.at(i).reserve(
				m_totalNumberOfGrains / m_totalNumberOfThreads + 1);
		for (int j = 0; j < m_totalNumberOfGrains / m_totalNumberOfThreads + 1; j++) {
			int id = j * m_totalNumberOfThreads + m_startIndex + i;
			if (id > m_totalNumberOfGrains)
				break;
			else
				m_threadWorkVectors.at(i).push_back(id);

		}
	}
}

vector<unsigned int>& IterativeGrainScheduler::getThreadWorkload(int threadID) {
return m_threadWorkVectors.at(threadID);
}
