#ifndef		__IGRAIN_SCHEDULER__
#define		__IGRAIN_SCHEDULER__

#include <vector>
#include "spoint.h"


using namespace std;

struct SPoint;
class IGrainScheduler
{
public:
	IGrainScheduler(){}
	virtual ~IGrainScheduler(){}
	virtual void buildGrainWorkloads(vector<vector<SPoint>>&, int) = 0;
	virtual std::vector<unsigned int>&	getThreadWorkload(int threadID) = 0;
};

#endif		//__IGRAIN_SCHEDULER__
