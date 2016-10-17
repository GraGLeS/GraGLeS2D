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

#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "ggLS.h"
#include "ExpandingVector.h"
#include "spoint.h"
#include "Settings.h"
#include <omp.h>
#include "misorientation.h"
#include "IGrainScheduler.h"
#include "dimensionalBuffer.h"

#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
// #define min(x,y) (x<y?x:y)
// #define max(x,y) (x>y?x:y)

using namespace voro;
using namespace std;

class LSbox;
class mathMethods;

class Quaternion;

/*!
 * \class grainhdl
 * \brief Class that manages the grain growth simulation.
 */
class grainhdl {
protected:
	int ngrains;
	double dt;
	double h;
	double m_Energy_deltaMAX;
	int realDomainSize;
	int ngridpoints;
	int grid_blowup;

	int Mode;
	IGrainScheduler* m_grainScheduler;

public:
	unsigned int currentNrGrains;
	mathMethods* mymath;
	unsigned int loop;
	//! control variable for research mode
	bool loadCurvature;
	unsigned int loadCurvatureLoop;
	bool convolutionCorrection;
	bool calcCentroid;
	bool calcRegression;
	E_RESEARCH_PROJECT project;
	bool constantE;
	vector<double> time;
	vector<double> totalenergy;
	vector<int> nr_grains;
	vector<double> MODF;

	//! A 2D vector which stores weights.
	vector<vector<double> > weightsMatrix;

	double ds;
	double *ST;
	double *part_pos;
	double delta;
	double *bunge;
	double deviation;
	double BoundaryGrainTube;
	double Realtime;
	double TimeSlope;
	double maxVol;
	MisorientationHdl* m_misOriHdl;
	Quaternion *TwinBoundary;
	DimensionalBuffer<int>* IDField;

	vector<LSbox*> grains;
	LSbox* boundary;

	grainhdl();
	~grainhdl();

	void setResearchAdjustments();
	void setResearchAdjustments(E_RESEARCH_PROJECT project);
	void setSimulationParameter();
	void read_HeaderCPG();

	void VOROMicrostructure();
	void readMicrostructureFromVertex();
	void readMicrostructure();
	void save_Full_Microstructure_for_Restart();
	void read_voxelized_microstructure();

	void createParamsForSim(const char* param_filename,
			const char* vertex_dum_filename = NULL);

	void find_neighbors();

	void distanceInitialisation();
	void convolution(double& plan_overhead);
	void createConvolutionPlans();
	void destroyConvolutionPlans();
	void save_conv_step();
	void comparison_box();
	void tweakIDLocal();

	void updateSecondOrderNeighbors();
	void level_set();
	void redistancing();

	virtual void run_sim();
	void save_NrGrainsStats();
	void clear_mem();
	void save_Texture();
	void save_id();
	void save_regLine();
	void save_Memory_Print();
	void plot_contour();
	void gridCoarsement();

	void saveAllContourLines();
	void save_NetworkPlot();
	void switchDistancebuffer();

	void saveSpecialContourEnergies(int id);
	void removeGrain(int id);
	// 	wrapper functions:

	void set_h(double hn);
	void set_realDomainSize(int realDomainSizen);
	//! Used if points are set manually
	int read_ScenarioPoints();
	void get_biggestGrainVol();
	void find_correctTimestepSize();
	inline LSbox* getGrainByID(unsigned int ID) {
		if (ID == 0)
			return boundary;
		else if (ID > 0 && ID < grains.size())
			return grains[ID];
		else
			return NULL;
	}

	inline long get_ngrains() {
		return ngrains;
	}
	inline int get_realDomainSize() {
		return realDomainSize;
	}
	inline int get_ngridpoints() {
		return ngridpoints;
	}
	inline double get_h() {
		return h;
	}
	inline int get_grid_blowup() {
		return grid_blowup;
	}
	inline int get_loop() {
		return loop;
	}
	inline double get_dt() {
		return dt;
	}
	inline double getBoundaryGrainTube() {
		return BoundaryGrainTube;
	}
	inline double get_ds() {
		return ds;
	}
	inline double get_maxVol() {
		return maxVol;
	}

protected:
	void initEnvironment();
	void initNUMABindings();
	void buildBoxVectors(vector<vector<SPoint>>& contours);
	void buildBoxVectors(vector<vector<SPoint>>& contours, vector<double>& q1,
			vector<double>& q2, vector<double>& q3, vector<double>& q4);
	void buildBoxVectors(int* ID, vector<vector<SPoint>>& contours,
			Quaternion* Quaternionen, double* StoredElasticEnergy);
	int m_ThreadPoolCount;
	vector<ExpandingVector<char> > m_ThreadMemPool;
};
#endif
