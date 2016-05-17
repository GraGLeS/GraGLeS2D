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




