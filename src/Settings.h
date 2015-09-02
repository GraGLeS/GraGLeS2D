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
#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>
#include "rapidxml.hpp"

/*!
 * \enum E_MICROSTRUCTURE_GEN_MODE
 * \brief Enumeration used to control how the microstructure in the simulation will be
 * generated.
 */
enum E_MICROSTRUCTURE_GEN_MODE
{
	E_READ_FROM_FILE,
	E_GENERATE_WITH_VORONOY,
	E_READ_VERTEX,
	E_GENERATE_TESTCASE,
	E_INVALID_VAL
};
/*!
 * \enum E_CONVOLUTION_MODE
 * \brief Enumeration used to control the convolution kernel.
 */
enum E_CONVOLUTION_MODE
{
	E_LAPLACE,
	E_LAPLACE_RITCHARDSON,
	E_GAUSSIAN,
	E_INVALID_VALUE
};

enum E_LATTICE_TYPE
{
	E_CUBIC,
	E_HEXAGONAL,
	E_INVALID_LATTICE
};

/*!
 * \enum E_RESEARCH_PROJECT
 * \brief Enumeration used to control the research project execution.
 */
enum E_RESEARCH_PROJECT
{
	E_NO_PROJECT,
	E_TRIPLE_JUNCTION_DRAG_SINGLE,
	E_TRIPLE_JUNCTION_DRAG_NETWORK,
	E_DEFAULT
};
/*!
 * \class Settings
 * \brief Class that holds all global settings that are simulation specific.
 */
class Settings
{
	public:
	static unsigned long StartTime;
	static unsigned int NumberOfParticles;
	static unsigned long NumberOfPointsPerGrain;
	static unsigned long NumberOfTimesteps;
	static unsigned long AnalysisTimestep;
	static unsigned long DiscreteSamplingRate;
	static unsigned long DomainBorderSize;
	static E_MICROSTRUCTURE_GEN_MODE MicrostructureGenMode;
	static E_RESEARCH_PROJECT ResearchProject;
	static std::string ReadFromFilename;
	static unsigned long LatticeType;
	static double HAGB;
	static double TriplePointDrag;
	static unsigned long UseMobilityModel;
	static bool IsIsotropicNetwork;
	static bool UseTexture;
	static double MaxMisOrientation;
	static bool ExecuteInParallel;
	static bool GridCoarsement;
	static bool ResearchMode;
	static double GridCoarsementGradient;
	static unsigned long MaximumNumberOfThreads;
	static E_CONVOLUTION_MODE ConvolutionMode;
	static double ConstantSectorRadius;
	static double InterpolatingSectorRadius;
	static unsigned long NeighbourTracking;


	static void initializeParameters(std::string filename = "");
	static rapidxml::xml_node<>* generateXMLParametersNode(rapidxml::xml_document<>* root, const char* filename, int loop, int grains);
};

#endif	//__SETTINGS_H__
