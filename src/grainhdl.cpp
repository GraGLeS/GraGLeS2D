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
#include "grainhdl.h"
#include "Settings.h"
#include "grahamScan.h"
#include "spoint.h"
#include "RTree.h"
#include "utilities.h"
#include <sys/time.h>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

using namespace rapidxml;

grainhdl::grainhdl() :
	m_ThreadPoolCount(0) {
}
grainhdl::~grainhdl() {
	delete mymath;
}

void grainhdl::setSimulationParameter() {
	initEnvironment();
	mymath = new mathMethods();
	// 	readInit();

	//! The NumberOfParticles passed via parameters.xml is altered
	//! for MicrostructureGenMode 4
	if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		Settings::NumberOfParticles = read_ScenarioPoints();
	}
	if (Settings::MicrostructureGenMode == E_READ_VERTEX) {
		FILE * getNumber;
		getNumber = fopen(Settings::ReadFromFilename.c_str(), "r");
		fscanf(getNumber, "%d\n", &Settings::NumberOfParticles);
	}
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	realDomainSize = sqrt(Settings::NumberOfParticles) * Settings::NumberOfPointsPerGrain; // half open container of VORO++
	discreteEnergyDistribution.resize(Settings::DiscreteSamplingRate);
	fill(discreteEnergyDistribution.begin(), discreteEnergyDistribution.end(),
			0);

	switch (Settings::ConvolutionMode) {
	case E_LAPLACE: {
		dt = 0.8 / double(realDomainSize * realDomainSize);
		break;
	}
	case E_LAPLACE_RITCHARDSON: {
		dt = 0.8 / double(realDomainSize * realDomainSize);
		break;
	}
	case E_GAUSSIAN: {
		dt = 0.75 / double(realDomainSize * realDomainSize);
		break;
	}
	default:
		break;
	}
	h = 1.0 / double(realDomainSize);

	//Recalculate the setting parameters for the sector radiuses
	Settings::ConstantSectorRadius *= h;
	Settings::InterpolatingSectorRadius *= h;

	delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
	tubeRadius = sqrt(2) * 1.5 * h + 0.001;
	grid_blowup = Settings::DomainBorderSize;
	BoundaryGrainTube = grid_blowup;
	ngridpoints = realDomainSize + (2 * grid_blowup);
	boundary = new LSbox(0, 0, 0, 0, this);
	// 	(*boundary).plot_box(false,2,"no.gnu");
	//!grains.resize(Settings::NumberOfParticles + 1);
	grains.resize(Settings::NumberOfParticles + 1);

	KernelNormalizationFactor = 2 * (Settings::NumberOfPointsPerGrain + (2
			* grid_blowup)) * (Settings::NumberOfPointsPerGrain + (2
			* grid_blowup));

	switch (Settings::MicrostructureGenMode) {
	case E_GENERATE_WITH_VORONOY: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = 15 * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;

		VOROMicrostructure();
		// 			generateRandomEnergy();
		break;
	}
	case E_READ_VERTEX: {
		bunge = NULL;
		deviation = 0;
		ST = new double[ngrains * ngrains];
		std::fill_n(ST, ngrains * ngrains, 0);
		readMicrostructureFromVertex();

		break;
	}
	case E_READ_FROM_FILE: {
		cout << "Reading file..." << endl;
		bunge = NULL;
		deviation = 0;

		ST = NULL;
		readMicrostructure();
		break;
	}
		//!
		//! This case handles the processing of
		//! grain construction by means of a file input
		//! with 2D point information
		//!
	case E_GENERATE_TESTCASE: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = 15 * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;
		VOROMicrostructure();
		// 			generateRandomEnergy();

		//! Read weights from a file
		if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE && !Settings::IsIsotropicNetwork) {

			//! Read from file in vector

			ifstream mapStream("surfaceTension.dat");

			if (!mapStream) {
				cerr << "File can't be opened." << endl;
				exit(-1);
			}

			while (mapStream) {
				string line;
				getline(mapStream, line);
				istringstream iss(line);
				vector<double> values;
				copy(istream_iterator<double> (iss),
						istream_iterator<double> (), back_insert_iterator<
								vector<double> > (values));
				if (values.size() > 0)
					weightsMatrix.push_back(values);
			}

			cout << "Line: " << weightsMatrix.size() << endl;
			cout << "Column (in the first line): " << weightsMatrix[0].size()
					<< endl;
		}

		break;
	}
	}

	GrainJunction::handler = this;

	// 	construct_boundary();
	//program options:
	cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
	cout << "Number of Grains: " << ngrains << endl;
	cout << "simulated Timesteps: " << Settings::NumberOfTimesteps << endl;
	cout << "DELTA TUBE: " << delta << endl;
	cout << "Timestepwidth " << dt << endl;
	cout << "Number of Gridpoints: " << ngridpoints << endl << endl;

	cout << endl << "******* start simulation: *******" << endl << endl;
}

void grainhdl::VOROMicrostructure() {

	stringstream filename, plotfiles;
	int cell_id;
	double x, y, z, rx, ry, rz;

	grains.resize(Settings::NumberOfParticles + 1);

	vector<LSbox*> local_grains;

	std::vector<LSbox*>::iterator itg;

	// stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
	part_pos = new double[3 * Settings::NumberOfParticles];

	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	if (randbedingung == false)
		realDomainSize -= 1;

	voronoicell_neighbor c;
	container con(0, 1, 0, 1, 0, 1, 5, 5, 5, randbedingung, randbedingung,
			randbedingung, 2);
	c_loop_all vl(con);

	//!
	//! Particles are added deliberately in the container according to the input file data.
	//!
	if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		FILE* pointSketch;
		pointSketch = fopen(Settings::ReadFromFilename.c_str(), "r");
		if (pointSketch == nullptr) {

			cout << "There does not exist a file";
			exit(1);
		}

		double pointX, pointY;

		for (int k = 0; k < ngrains; k++) {

			fscanf(pointSketch, "%lf\t%lf\n", &pointX, &pointY);
			con.put(k,pointX,pointY,0);
		}

		fclose(pointSketch);
	}

	/**********************************************************/
	// Randomly add particles into the container
	if (Settings::MicrostructureGenMode != E_GENERATE_TESTCASE) {
		for (int i = 0; i < ngrains; i++) {
			x = rnd();
			y = rnd();
			z = 0;
			con.put(i, x, y, z);
		}
	}
	/**********************************************************/

	for (int i = 0; i < realDomainSize; i++)
		for (int j = 0; j < realDomainSize; j++) {
			x = double(i * h);
			y = double(j * h); // only point within the domain
			if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
				cell_id++;
				part_pos[3 * (cell_id - 1)] = rx;
				part_pos[3 * (cell_id - 1) + 1] = ry;
				part_pos[3 * (cell_id - 1) + 2] = rz;
			}

			else
				fprintf(stderr, "# find_voronoi_cell error for %g %g 0\n", x, y);

		}

	vector<vector<SPoint>> m_initialContours;
	m_initialContours.resize(Settings::NumberOfParticles+1);

	if (vl.start())
		do {
			con.compute_cell(c, vl);
			int box_id = vl.pid() + 1;
			GrahamScan scanner(c, box_id, part_pos);
			scanner.generateCovnexHull(m_initialContours[box_id]);


		} while (vl.inc());

	buildBoxVectors(m_initialContours);

	delete[] part_pos;
}

void grainhdl::readMicrostructure() {
	FILE * levelset;
	levelset = fopen(Settings::ReadFromFilename.c_str(), "r");
	if (levelset == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}
	int id;
	double q1, q2, q3, q4,xl, yl;

	int nvertices;
	//double* vertices = new double[1000];
	vector<double> vertices;
	vertices.resize(2048);

	for (int nn = 1; nn <= ngrains; nn++) {

		fscanf(levelset, "%d\t %d\t %lf\t %lf\t%lf\t%lf\n", &id, &nvertices,
				&q1, &q2, &q3, &q4);
		vertices.resize(nvertices*2);
		for (int j = 0; j < nvertices; j++) {
			fscanf(levelset, "%lf\t %lf\n", &xl, &yl);
			vertices[2 * j] = xl;
			vertices[(2 * j) + 1] = yl;
		}
		fscanf(levelset, "\n");
		LSbox* newBox =
				new LSbox(id, nvertices, &vertices[0], q1, q2, q3, q4, this);
		grains[nn] = newBox;
	}
	fclose(levelset);
}

void grainhdl::readMicrostructureFromVertex() {
	FILE * levelset;
	levelset = fopen(Settings::ReadFromFilename.c_str(), "r");

	long id;
	int nedges;
	double phi1, PHI, phi2, xr, yr, xl, yl;
	double* edges;

	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;

	grains.resize(ngrains + 1);

	for (int nn = 0; nn < ngrains; nn++) {

		fscanf(levelset, "%ld\t %d\t %lf\t %lf\t%lf\n", &id, &nedges, &phi1,
				&PHI, &phi2);
		edges = new double[nedges * 4];
		cout << id << " || " << nedges << " || " << phi1 << " || " << PHI
				<< " || " << phi2 << endl;

		for (int j = 0; j < nedges; j++) {
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);
			cout << xl << " ||\t " << yl << " ||\t " << xr << " ||\t " << yr
					<< " ||\t " << endl;
			int k = 4 * j;
			edges[k] = xl;
			edges[k + 1] = yl;
			edges[k + 2] = xr;
			edges[k + 3] = yr;
		}

		LSbox* newBox = new LSbox(id, nedges, edges, phi1, PHI, phi2, this);
		grains[id] = newBox;

		delete[] edges;
	}

	ST = new double[ngrains * ngrains]; //Create ST array and fill with zeros
	std::fill_n(ST, ngrains * ngrains, 0);

	for (int i = 0; i < ngrains; i++) {
		double buffer;
		fscanf(levelset, "%lf\t", &buffer);
		for (int j = 0; j < ngrains; j++) {
			while (j < i) {
				fscanf(levelset, "%lf\t", &buffer);
				j++;
			}
			fscanf(levelset, "%lf\t", &buffer);
			ST[j + (ngrains * i)] = (double) buffer;
			ST[i + (ngrains * j)] = ST[j + (ngrains * i)];
		}
		fscanf(levelset, "\n");
	}
	fclose(levelset);

	for (int i = 0; i < ngrains; i++) {
		for (int j = 0; j < ngrains; j++) {
			cout << ST[i + (ngrains * j)] << "  \t";
		}
		cout << endl;
	}
}

void grainhdl::distanceInitialisation() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if(grains[id] != NULL)
					grains[id]->calculateDistanceFunction();
			}
		}
}

void grainhdl::convolution(double& plan_overhead) {
	double timer;
	timeval t;

	plan_overhead = 0;
	gettimeofday(&t, NULL);
	timer = t.tv_sec + t.tv_usec/1000000;
	createConvolutionPlans();
	gettimeofday(&t, NULL);
	plan_overhead += t.tv_sec + t.tv_usec/1000000 - timer;
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if(grains[id] != NULL)
					grains[id]->executeConvolution(m_ThreadMemPool[omp_get_thread_num()]);
			}
		}
	gettimeofday(&t, NULL);
	timer = t.tv_sec + t.tv_usec/1000000;
	destroyConvolutionPlans();
	gettimeofday(&t, NULL);
	plan_overhead += t.tv_sec + t.tv_usec/1000000 - timer;
}

void grainhdl::createConvolutionPlans() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if(grains[id] != NULL)
					grains[id]->preallocateMemory(m_ThreadMemPool[omp_get_thread_num()]);
			}
		}
	for(unsigned int i=1; i<grains.size(); i++)
	{
		if(grains[i] != NULL)
			grains[i]->makeFFTPlans(m_ThreadMemPool[(i - 1)%Settings::MaximumNumberOfThreads]);
	}
}

void grainhdl::destroyConvolutionPlans() {
	for(unsigned int i=1; i<grains.size(); i++)
	{
		if(grains[i] != NULL)
			grains[i]->destroyFFTWs();
	}
}

void grainhdl::comparison_box() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if(grains[id] != NULL)
				{
					grains[id]->executeComparison();
					grains[id]->executeSetComparison();
				}
			}
		}
}

void grainhdl::level_set() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				if (grains[id]->grainExists() == false) {
					delete grains[id];
					grains[id] = NULL;
				} else
					grains[id]->extractContour();
			}
		}
}

void grainhdl::redistancing() {
	currentNrGrains = 0;
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
#pragma omp atomic
				currentNrGrains += 1;
				grains[id]->computeVolumeAndEnergy();
				grains[id]->executeRedistancing();
			}
		}
}

void grainhdl::save_texture() {
	ofstream myfile;
	FILE* enLenDis;
	stringstream filename;

	int numberGrains = 0;
	double totalLength = 0;
	double total_energy = 0.0;

	filename << "Texture" << "_" << loop << ".ori";
	myfile.open(filename.str());

	filename.str("");
	filename << "EnergyLengthDistribution_" << loop << ".txt";
	enLenDis = fopen(filename.str().c_str(), "w");

	double dh = Settings::HAGB / (double) Settings::DiscreteSamplingRate;
	double euler[3];
	vector<characteristics>::iterator it2;
	vector<LSbox*>::iterator it;

	std::fill(discreteEnergyDistribution.begin(),
			discreteEnergyDistribution.end(), 0.0);

	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it != NULL && (*it)->grainExists() == true && (*it)->isMotionRegular() == true) {
			numberGrains++;
			total_energy += (*it)->getEnergy();

			(*mymath).quaternion2Euler((*it)->getOrientationQuat(), euler);

			totalLength += (*it)->getPerimeter() * 0.5;
			map<int, double>& currentEnergyDistribution = (*it)->getDiscreteEnergyDistribution();
			for(const auto& iterator : currentEnergyDistribution)
			{
				if(iterator.first < discreteEnergyDistribution.size())
				{
					discreteEnergyDistribution[iterator.first] += iterator.second / 2;
				}
			}
			//! If ResearchMode is activated additional data (e.g. the mean (euclidean) triple-junction-distance) is stored in the Texture files
			if(Settings::ResearchMode) {
				double mean_a = (*it)->getMeanA();
				//! Determine the mean number "m" of the direct neighbors' number of faces

				double mean_m = (*it)->getMeanM();

				myfile << (*it)->getID() << '\t' << (*it)->getDirectNeighbourCount() << '\t'
					   << (*it)->intersectsBoundaryGrain() << '\t' << (*it)->getVolume() << '\t'
					   << (*it)->getMeanDa() / Settings::AnalysisTimestep << '\t'
					   << (*it)->getPerimeter() << '\t' << (*it)->getEnergy() << '\t'
					   << mean_a << '\t' << mean_m << '\t'
					   << euler[0] << '\t' << euler[1] << '\t' << euler[2] ;

			}
			else {
				myfile << (*it)->getID() << '\t' << (*it)->getDirectNeighbourCount() << '\t'
					   << (*it)->intersectsBoundaryGrain() << '\t' << (*it)->getVolume() << '\t'
					   << (*it)->getMeanDa() / Settings::AnalysisTimestep << '\t'
					   << (*it)->getPerimeter() << '\t' << (*it)->getEnergy() << '\t'
					   << euler[0] << '\t' << euler[1] << '\t' << euler[2] ;
			}

			if(Settings::NeighbourTracking){
				vector<int> ids =(*it)->getDirectNeighbourIDs();
				vector<double> lengths = (*it)->getGBLengths();
				if (ids.size() != lengths.size()) cout << "length differences in ID and Length array! " << endl ;
				for(int i=0; i < (*it)->getDirectNeighbourCount(); i++){
				myfile << '\t' << ids[i] << '\t' << lengths[i] << '\t';
				}
			}

			myfile << endl;
			double sum = 0.0;
			for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
				sum += discreteEnergyDistribution[i];
			}
			for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
				if (discreteEnergyDistribution[i] > 0)
					discreteEnergyDistribution[i]
							= discreteEnergyDistribution[i] / sum;
			}

		}
	}

	if (!Settings::IsIsotropicNetwork) {
		for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
			fprintf(enLenDis, "%lf\t%lf\n", (float) (dh * (i + 1)),
					(float) discreteEnergyDistribution[i]);
			printf("%lf\t%lf\n", (float) (dh * (i + 1)),
					(float) discreteEnergyDistribution[i]);
		}
	}
	totalenergy.push_back(0.5 * total_energy);
	nr_grains.push_back(numberGrains);
	time.push_back(simulationTime);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << numberGrains
			<< endl;
	cout << "Amount of free Energy in the Network :" << 0.5 * total_energy
			<< endl;
	cout << "Total GB Length in Network :" << totalLength << endl << endl
			<< endl;
	myfile.close();
	fclose(enLenDis);

}

double convo_time=0;
double comparison_time=0;
double levelset_time=0;
double redistancing_time=0;
double plan_overhead = 0;
void grainhdl::run_sim() {
	timeval time;
	double timer;
	double overhead;
	distanceInitialisation();
	simulationTime = 0;
	find_neighbors();
	for (loop = Settings::StartTime; loop <= Settings::StartTime
			+ Settings::NumberOfTimesteps; loop++) {
		gridCoarsement();

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec/1000000.0;
		convolution(overhead);
		plan_overhead += overhead;
		gettimeofday(&time, NULL);
		convo_time += time.tv_sec + time.tv_usec/1000000.0 - timer;

		switchDistancebuffer();
		updateSecondOrderNeighbors();
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec/1000000.0;
		comparison_box();
		gettimeofday(&time, NULL);
		comparison_time += time.tv_sec + time.tv_usec/1000000.0 - timer;

		switchDistancebuffer();
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec/1000000.0;
		level_set();
		gettimeofday(&time, NULL);
		levelset_time += time.tv_sec + time.tv_usec/1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec/1000000.0;
		redistancing();
		gettimeofday(&time, NULL);
		redistancing_time += time.tv_sec + time.tv_usec/1000000.0 - timer;

		if (((loop - Settings::StartTime) % int(Settings::AnalysisTimestep))
				== 0 || loop == Settings::NumberOfTimesteps) {
			saveAllContourEnergies();
			save_texture();
			save_sim();
			//! With activated ResearchMode the id to centroid assignment is plotted
			if(Settings::ResearchMode && calcCentroid)
				save_id();
			if (loop == Settings::NumberOfTimesteps)
				saveMicrostructure();

//			stringstream name;
//			name << "MemoryPrint_T"<<loop<<".txt";
//			ofstream file;
//			file.open(name.str());
//			for(int i=1; i< Settings::NumberOfParticles; i++){
//				if(grains[i] != NULL)
//				{
//					grains[i]->outputMemoryUsage(file);
//				}
//			}
//			file.close();
		}
		simulationTime += dt;
		if ((double)currentNrGrains/ (double)Settings::NumberOfParticles < 0.001 ) {
			cout	<< "Network has coarsed to less than 3% of the population. "
					<< "Remaining Grains: " << currentNrGrains
					<<" Break and save."
					<< endl;
			break;
		}
	}
	// 	utils::CreateMakeGif();
	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << simulationTime << endl;
	cout << "Detailed timings: "<< endl;
	cout << "Convolution time: "<< convo_time << endl;
	cout << "     Of which plan overhead is: " << plan_overhead << endl;
	cout << "Comparison time: "<< comparison_time << endl;
	cout << "Redistancing time: "<< redistancing_time << endl;
	cout << "Levelset time: "<< levelset_time << endl;
}

void grainhdl::saveMicrostructure() {
	stringstream param_xml_name;
	param_xml_name << "PARAMETERS_SIM_NONAME_TIMESTEP_" << loop << "_GRAINS_"
			<< currentNrGrains << ".xml";
	stringstream vertex_dump_name;
	vertex_dump_name << "NETWORK_NONAME_TIMESTEP_" << loop << "_GRAINS_"
			<< currentNrGrains << ".dat";

	createParamsForSim(param_xml_name.str().c_str(),
			vertex_dump_name.str().c_str());

	ofstream output;
	output.open(vertex_dump_name.str());
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL || (*it)->grainExists() == false)
			continue;
		(*it)->plot_full_grain(loop, false, &output, true);
	}
	output.close();

}
void grainhdl::createParamsForSim(const char* param_filename,
		const char* vertex_dump_filename) {
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(doc_tree.allocate_attribute("encoding",
			"utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(Settings::generateXMLParametersNode(&doc_tree,
			vertex_dump_filename, loop, currentNrGrains));
	ofstream output;
	output.open(param_filename);
	output << doc_tree;
	output.close();

}

void grainhdl::save_sim() {
	// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);
	ofstream myfile;
	myfile.open("NrGrains&EnergyStatistics.txt");
	for (unsigned int i = 0; i < nr_grains.size(); i++) {
		myfile << time[i] << "\t";
		myfile << nr_grains[i] << "\t";
		myfile << totalenergy[i] << endl;
	}
	myfile.close();

	// 	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
	//cout << "number of distanzmatrices: "<< domains.size() << endl;

}

void grainhdl::updateSecondOrderNeighbors() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->computeSecondOrderNeighbours();
			}
		}
}

void grainhdl::find_neighbors() {
	RTree<unsigned int, int, 2, float> tree;
	int min[2], max[2];
	for(unsigned int i=1; i<=Settings::NumberOfParticles; i++)
	{
		if(grains[i] == NULL)
			continue;
		min[0] = grains[i]->getMinX(); min[1] = grains[i]->getMinY();
		max[0] = grains[i]->getMaxX(); max[1] = grains[i]->getMaxY();
		tree.Insert(min, max, i);
	}

	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->computeDirectNeighbours(tree);
			}
		}
}



void grainhdl::saveSpecialContourEnergies(int id) {
	if (grains[id] == NULL)
		return;
	grains[id]->plot_box_contour(loop, true);
}

void grainhdl::saveAllContourEnergies() {
	ofstream output;
	stringstream filename;
	filename << "Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		(*it)->plot_box_contour(loop, true, &output, true);
	}
	output.close();
}

void grainhdl::save_id() {
	ofstream output;
	stringstream filename;
	filename << "Id_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());


	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		SPoint centroid = (*it)->getCentroid();
		output << centroid.x << "\t" << centroid.y << "\t" << (*it)->getID() << endl;
	}
	output.close();

//! Plot junction points
	stringstream filename2;
	filename2 << "Junctions_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename2.str());

	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		else
			(*it)->plot_grain_junctions(loop, &output);
	}
	output.close();

}

void grainhdl::save_regLine() {
	ofstream output;
	stringstream filename;
	filename << "RegLine_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		for(unsigned int j = 0; j < (*it)->getRegressionPoints().size(); j++) {

			output << (*it)->getRegressionPoints()[j].x << "\t" <<(*it)->getRegressionPoints()[j].y << "\t" << endl;
			if ((j-2)%3 == 0)
				output << "\n" << endl;
		}


		//!output << "\n" << endl;
	}
	output.close();
}

void grainhdl::saveAllContourLines() {
	stringstream filename;
	filename << "NetworkAtTime_" << loop << ".gnu";
	ofstream dateiname;
	dateiname.open(filename.str());
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		(*it)->plot_box_contour(loop, false);
	}
	dateiname.close();
}
void grainhdl::removeGrain(int id) {
	grains[id] = NULL;
}

void grainhdl::switchDistancebuffer() {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->switchInNOut();
			}
		}
}

void grainhdl::gridCoarsement() {
	if ((double)currentNrGrains/(double)ngrains < Settings::GridCoarsementGradient
			&& loop != 0 && Settings::GridCoarsement) {
		int newSize = sqrt(currentNrGrains) * Settings::NumberOfPointsPerGrain;
		cout << "coarsing the current grid in Timestep: "<< loop <<endl;
		cout<< "newSize :" << newSize << endl << endl;
		unsigned int j;
	#pragma omp parallel for private(j)
		for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
			for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
			{
				if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
				{
					int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
					if (grains[id] == NULL)
						continue;
					grains[id]->resizeGrid(newSize);
				}

			}
		realDomainSize = newSize;
		delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
		ngridpoints = realDomainSize + 2 * grid_blowup;
		h = 1.0 / realDomainSize;
		//! DISCREPANCY: Compare to the application of dt in the convolution, time decreasing factor 0.8
		dt = 1.0 / double(realDomainSize * realDomainSize);
		ngrains= currentNrGrains;
	#pragma omp parallel for private(j)
		for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
			for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
			{
				if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
				{
					int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
					if (grains[id] == NULL)
						continue;
					grains[id]->recalculateIDLocal();
				}
			}
	#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < grains.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				if (grains[id] == NULL)
					continue;
				grains[id]->extractContour();
			}
		}

	} else {
		switchDistancebuffer();
	}
}

void grainhdl::clear_mem() {
	if (ST != NULL) {
		delete[] ST;
	}
}

void grainhdl::initEnvironment() {
	//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

	m_ThreadPoolCount = Settings::MaximumNumberOfThreads;
	m_ThreadMemPool.resize(m_ThreadPoolCount);
	for (unsigned int i = 0; i < m_ThreadMemPool.size(); i++) {
		double max_size = Settings::NumberOfPointsPerGrain
				* Settings::NumberOfPointsPerGrain * 50;
		int power_of_two = 1 << (int)(ceil(log2(max_size))+0.5);
		//!int power_of_two = 1 << (int) (ceil(log2(2<<20)) + 0.5); //!27
		m_ThreadMemPool[i].resize(power_of_two);
	}
}

void grainhdl::buildBoxVectors(vector<vector<SPoint>>& contours) {
	unsigned int j;
#pragma omp parallel for private(j)
	for(unsigned int i=0; i< Settings::MaximumNumberOfThreads; i++)
		for(j=0; j < Settings::NumberOfParticles/Settings::MaximumNumberOfThreads + 1; j++)
		{
			if(j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num() < contours.size())
			{
				int id = j*Settings::MaximumNumberOfThreads + 1 + omp_get_thread_num();
				LSbox* grain = new LSbox(id, contours[id], this);
				grains[id] = grain;
			}
		}
}

void grainhdl::set_h(double hn) {
	h = hn;
}
void grainhdl::set_realDomainSize(int realDomainSizen) {
	realDomainSize = realDomainSizen;
	ngridpoints = realDomainSize + 2 * grid_blowup;
}
/**
 * This function analyzes the input file for MicrostructureGenMode 4.
 * The amount of lines of the input file is determined. This number
 * indicates the number of points specified in the file, i.e the number of
 * grains. This means one pair of x-y-coordinates each in every line. A
 * tabulator is used as the separator.
 *
 * @return the amount of points in the input file
 */
int grainhdl::read_ScenarioPoints() {

	int counter = 0;
	string line;
	ifstream reader(Settings::ReadFromFilename.c_str());
	while (std::getline(reader, line)) {
		counter++;
	}
	return counter;
}

/**
 * This piece of code cares for the proper program behavior as
 * intended by the current research project.
 */
void grainhdl::setResearchAdjustments(E_RESEARCH_PROJECT pro) {

	//! Some minor control variables to enhance speed if necessary,
	//! It is important to check their properties before starting a
	//! large simulation project.
	//! Activate calibration factor in convolution corrector step
	convolutionCorrection = false;
	//! Calculate regression lines?
	calcRegression = false;
	//! Calculate centroids?
	calcCentroid = false;
	//! Is there a need to load curvature in the first timesteps?
	loadCurvature = false;
	//! As the reference the observed values of the four-sided-grain case
	//! is used.
	loadCurvatureLoop = ((Settings::NumberOfPointsPerGrain/40.0)*(Settings::NumberOfPointsPerGrain/40.0)*200);
	//! See the ENUM declaration to get the matching of integer and project
	if (pro == 0)
		project = pro;

	//! This project is characterized by constant energies for every
	//! grain boundary in the network. The few exceptions are the
	//! junctions sections on the domain boundaries which need to be
	//! fixed to their position (weight = 0) because of small deviation
	//! from symmetry in the algorithm.
	//!
	//! Application: Investigation of triple junction and grain boundary
	//! regime by means of Lambda.
	//!
	//! Changes are made in: junction.cpp
	//!
	//! Works with: MicrostructureGenMode == E_READ_VERTEX
	if (pro == 1)
		project = pro;

	//! This project is characterized by constant energies for every
	//! grain boundary in the network.
	//!
	//! Application: Investigation of triple junction and grain boundary
	//! regime by means of Lambda.
	//!
	//! Changes are made in: junction.cpp
	//!
	//! Works with: MicrostructureGenMode == E_GENERATE_WITH_VORONOY,
	//! 			MicrostructureGenMode == E_READ_FROM_FILE
	if (pro == 2)
		project = pro;
}

void grainhdl::setResearchAdjustments() {
	return;
}

