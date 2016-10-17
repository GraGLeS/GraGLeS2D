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
#ifndef BOX_H
#define BOX_H

#include "ggLS.h"
#include "dimensionalBufferIDLocal.h"
#include "dimensionalBufferReal.h"
#include "junction.h"
#include "dimensionalBuffer.h"
#include "pooledDimensionalBufferDouble.h"
#include "spoint.h"
#include "contourSector.h"
#include "minimalisticBoundary.h"
#include "grainBoundary.h"
#include "charasteristic.h"
#include "RTree.h"
#include "Quaternion.h"

#ifdef USE_MKL
#include "mkl_dfti.h"
#endif

using namespace std;

class MisorientationHdl;
class grainhdl;
class DimensionalBufferReal;
class MarchingSquaresAlgorithm;
class Settings;
class Quaternion;
class ExplicitGrainBoundary;
struct SPoint;
/*!
 * \struct VolEvolution
 * \brief This structure stores the evolution of area and the corresponding number of vertices.
 */

/*!
 * \class LSbox
 * \brief Class encapsulating a Level Set Box.
 *
 * LSbox class contains the coordinates of the box in the actual grid. <br>
 * For each point that the LSbox covers, it stores: <br>
 * - Distances to the actual grain boundary. <br>
 * - The ID of the closest grain.
 *
 * The class also stores the coordinates of the LSbox, Euler angles that represent
 * the orientation, the volume of the grain, <br>
 * the energy of the grain and a pointer to the \b grainhdl object.
 *
 */
class LSbox {
	friend class GrainJunction;
private:
	unsigned int m_ID;
	bool m_exists;
	grainhdl* m_grainHandler;
	bool m_isMotionRegular;
	bool m_intersectsBoundaryGrain;
	DimensionalBufferIDLocal m_IDLocal;
	double m_meanDa;
	Quaternion* m_orientationQuat;
	double m_volume;
	double m_energy;
	double m_perimeter;
	int m_newXMin;
	int m_newXMax;
	int m_newYMin;
	int m_newYMax;
	double m_StoredElasticEnergy;
	double m_magneticEnergy;
	SPoint m_centroid;
	vector<SPoint> m_regressionPoints;
	vector<SPoint> m_triangleCetroid;
	MinimalisticBoundary m_minimalBoundary;
	DimensionalBufferReal* m_inputDistance;
	DimensionalBufferReal* m_outputDistance;

	vector<LSbox*> m_comparisonList;
	vector<LSbox*> m_secondOrderNeighbours;

#ifdef USE_FFTW
	fftwp_plan m_backwardsPlan;
	fftwp_plan m_forwardPlan;
#elif USE_MKL
	DFTI_DESCRIPTOR_HANDLE m_handle;
	DFTI_DESCRIPTOR_HANDLE m_b_handle;
	MKL_LONG m_dimensions[2];
	MKL_LONG m_f_input_strides[3];
	MKL_LONG m_f_output_strides[3];
	MKL_LONG m_b_input_strides[3];
	MKL_LONG m_b_output_strides[3];
#endif
protected:
	ExplicitGrainBoundary m_grainBoundary;
public:
	//Constructors to document
	LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner);
	LSbox(int id, const vector<SPoint>& vertices, double q1, double q2,
			double q3, double q4, grainhdl* owner);
	LSbox(int aID, vector<SPoint>& contour, grainhdl* owner);
	LSbox(int id, int nvertex, double* vertices, double phi1, double PHI,
			double phi2, grainhdl* owner);
	LSbox(int id, const vector<SPoint>& vertices, Quaternion ori,
			double StoredElasticEnergy, grainhdl* owner);
	//Dtors
	~LSbox();
	void calculateDistanceFunction();
	void executeRedistancing();
	void executeExactRedist();
	void extractContour();
	void markAsInvalidMotion();
	inline bool isMotionRegular() const {
		return m_isMotionRegular;
	}
	void executeComparison();
	double getDistanceFromInputBuff(int i, int j);
	void executeSetComparison();
	void computeSecondOrderNeighbours();
	void computeDirectNeighbours(
			const RTree<unsigned int, int, 2, float>& rtree);
	double computeVolume();
	bool isMotionRegular(int old_contour, int new_contour, double old_volume,
			double new_volume);
	void clearContourGrainArea();
	void computeVolumeAndEnergy();
	double getGBEnergyTimesGBMobility(int i, int j);
	double getGBEnergyTimesGBMobility(LSbox* neighbour);
	double getGBEnergy(LSbox* neighbour);
	double GBmobilityModel(double thetaMis, LSbox* candidate);
	double getWeigthFromHandler(int i, int j);
	void constructBoundarySectors();
	double getWeight(int i, int j, bool minimal = false);

	void marchingSquares(DimensionalBufferReal* which);
	vector<int> getDirectNeighbourIDs();
	vector<double> getGBLengths();
	map<int, double>& getlocalMODF() {
		return m_grainBoundary.getlocalMODF();
	}
	bool checkIntersection(LSbox* box2);
	void executeConvolution(ExpandingVector<char>& mem_pool);
	void reizeIDLocalToDistanceBuffer();
	void recalculateIDLocal();
	void setIDLocal(int ID);

	//Debug printing functions
	void plot_box_contour(int timestep = -1, bool plot_energy = false,
			ofstream* dest_file = NULL, bool absCoordinates = false,
			int threadID = 0);
	void plot_full_grain(int timestep = -1, bool plot_energy = false,
			ofstream* dest_file = NULL, bool absCoordinates = false);
	void plot_box_parameters(ofstream* dest_file = NULL);
	void plot_grain_junctions(int timestep = -1, ofstream* dest_file = NULL,
			bool absCoordinates = false);
	void plot_box(bool distanceplot, int select, string simstep, bool local);

	double computeMisorientation(LSbox* grain_2);
	double computeMisorientation(unsigned int grainID);
	void resizeGrid(int newSize);

	void preallocateMemory(ExpandingVector<char>& memory_dump);

#ifdef USE_FFTW
	void makeFFTPlans(double *in, double* out,fftw_complex *fftTemp, fftw_plan *fftplan1, fftw_plan *fftplan2);
	void makeFFTPlans(float *in, float* out,fftwf_complex *fftTemp, fftwf_plan *fftplan1, fftwf_plan *fftplan2);
	void makeFFTPlans(ExpandingVector<char>& memory_dump);
	void convolutionGeneratorFFTW(fftwp_complex *fftTemp, fftwp_plan fftplan1, fftwp_plan fftplan2);
	void executeFFTW(fftw_plan fftplan);
	void executeFFTW(fftwf_plan fftplan);
	void destroyFFTWs();
#elif defined USE_MKL
	void convolutionGeneratorMKL(MKL_Complex16* fftTemp);
#endif
	void switchInNOut();
	//todo: refactor with a proper name
	void boundaryCondition();
	inline bool intersectsBoundaryGrain() const {
		return m_intersectsBoundaryGrain;
	}
	void updateFirstOrderNeigbors();

	bool isNeighbour(LSbox* candidate);
	bool BoundaryIntersection();
	void calculateMagneticEnergy();

	//todo: Analyze if function is required
	void measureAngles(vector<double>& turningAngles,
			vector<GrainJunction> junctions);
	vector<double> linearRegression(vector<SPoint>& points2D);
	void calculateTriangleCentroid(vector<SPoint>& triangleCentroid,
			vector<SPoint> triangle);
	void calculateCentroid(SPoint& centroid, vector<GrainJunction> junctions);

	double MisoriToTwinBoundary(LSbox* candidate);
	double GBEnergyReadShockley(double theta, LSbox* candidate);
	double get_h();

	void outputMemoryUsage(ofstream& output);
	TextureData collectTextureData();

	inline vector<Face>* get_Faces(){
		return m_grainBoundary.getFaces();
	}
	inline int getDirectNeighbourCount() {
		return m_grainBoundary.getDirectNeighboursCount();
	}
	inline bool grainExists() const {
		return m_exists;
	}
	inline double getVolume() const {
		return m_volume;
	}
	inline double getEnergy() const {
		return m_energy;
	}
	inline double getPerimeter() const {
		return m_perimeter;
	}
	inline unsigned int getID() const {
		return m_ID;
	}
	inline int getMinX() const {
		return m_outputDistance->getMinX();
	}
	inline int getMaxX() const {
		return m_outputDistance->getMaxX();
	}
	inline int getMinY() const {
		return m_outputDistance->getMinY();
	}
	inline int getMaxY() const {
		return m_outputDistance->getMaxY();
	}
	inline SPoint getCentroid() const {
		return m_centroid;
	}
	inline double getMeanDa() const {
		return m_meanDa;
	}
	inline double getMeanM() const {
		return m_grainBoundary.getMeanM();
	}
	inline double getMeanA() const {
		return m_grainBoundary.getMeanA();
	}
	inline const Quaternion* getOrientationQuat() {
		return m_orientationQuat;
	}
	inline const vector<SPoint>& getRegressionPoints() const {
		return m_regressionPoints;
	}
	inline double get_StoredElasticEnergy() {
		return m_StoredElasticEnergy;
	}
	inline double get_magneticEnergy() {
		return m_magneticEnergy;
	}
	LSbox* getNeighbourAt(int i, int j);
};
#endif
