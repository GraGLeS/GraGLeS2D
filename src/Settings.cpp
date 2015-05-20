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
#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;
using namespace rapidxml;

//Initializing the static setting variables
unsigned long Settings::StartTime = 0;
unsigned int Settings::NumberOfParticles = 0;
unsigned long Settings::NumberOfPointsPerGrain = 0;
unsigned long Settings::NumberOfTimesteps = 0;
unsigned long Settings::AnalysisTimestep = 0;
unsigned long Settings::DiscreteSamplingRate = 0;
unsigned long Settings::DomainBorderSize = 0;
unsigned long Settings::MaximumNumberOfThreads = 0;
double Settings::GridCoarsementGradient = 0.95;
E_CONVOLUTION_MODE Settings::ConvolutionMode = E_INVALID_VALUE;
E_MICROSTRUCTURE_GEN_MODE Settings::MicrostructureGenMode = E_INVALID_VAL;
E_RESEARCH_PROJECT Settings::ResearchProject = E_DEFAULT;
string Settings::ReadFromFilename;
double Settings::HAGB = 0.0;
double Settings::TriplePointDrag = 0.0;
bool Settings::UseMobilityFactor = false;
bool Settings::IsIsotropicNetwork = false;
bool Settings::UseTexture = false;
bool Settings::ExecuteInParallel = false;
bool Settings::GridCoarsement = false;
//! Enables the generation and capture of a wider spectrum of analysing data
bool Settings::ResearchMode = false;
double Settings::ConstantSectorRadius = 0.0;
double Settings::InterpolatingSectorRadius = 0.0;
unsigned long Settings::NeighbourTracking = 0;
unsigned long Settings::DislocationEnergy = 0;

void Settings::initializeParameters(string filename) {
	if (0 == filename.compare(""))
		filename = string("parameters.xml");
	ifstream file(filename);
	if (file.fail()) {
		cout << "Unable to locate simulations parameters. Will now halt !"
				<< endl;
		exit(2);
	}
	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	try {
		tree.parse<0> (&xmlDocument[0]);
	} catch (parse_error& Error) {
		cout << filename.c_str() << "is not a valid XML!" << endl;
		cout << "Exception is " << Error.what() << endl;
	}
	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		cout << "Root node is not 'Parameters'! Will now halt!" << endl;
		exit(2);
	}

	//Now read all parameters if present
	if (0 != rootNode->first_node("StartTime")) {
		StartTime = std::stoul(rootNode->first_node("StartTime")->value());
	}
	if (0 != rootNode->first_node("NumberOfParticles")) {
		NumberOfParticles = std::stoul(
				rootNode->first_node("NumberOfParticles")->value());
	}
	if (0 != rootNode->first_node("NumberOfPointsPerGrain")) {
		NumberOfPointsPerGrain = std::stoul(rootNode->first_node(
				"NumberOfPointsPerGrain")->value());
	}
	if (0 != rootNode->first_node("AnalysisTimestep")) {
		AnalysisTimestep = std::stoul(
				rootNode->first_node("AnalysisTimestep")->value());
	}
	if (0 != rootNode->first_node("NumberOfTimesteps")) {
		NumberOfTimesteps = std::stoul(
				rootNode->first_node("NumberOfTimesteps")->value());
	}
	if (0 != rootNode->first_node("DiscreteSamplingRate")) {
		DiscreteSamplingRate = std::stoul(rootNode->first_node(
				"DiscreteSamplingRate")->value());
	}
	if (0 != rootNode->first_node("DomainBorderSize")) {
		DomainBorderSize = std::stoul(
				rootNode->first_node("DomainBorderSize")->value());
	}
	if (0 != rootNode->first_node("MicrostructureGenMode")) {
		MicrostructureGenMode = (E_MICROSTRUCTURE_GEN_MODE) std::stoi(
				rootNode->first_node("MicrostructureGenMode")->value());
		if (MicrostructureGenMode >= E_INVALID_VAL)
			MicrostructureGenMode = E_INVALID_VAL;
	}
	if (0 != rootNode->first_node("ReadFromFilename")) {
		ReadFromFilename = rootNode->first_node("ReadFromFilename")->value();
	}
	if (0 != rootNode->first_node("HAGB")) {
		HAGB = std::stod(rootNode->first_node("HAGB")->value());
	}
	if (0 != rootNode->first_node("TriplePointDrag")) {
		TriplePointDrag = std::stod(
				rootNode->first_node("TriplePointDrag")->value());
	}
	if (0 != rootNode->first_node("UseMobilityFactor")) {
		UseMobilityFactor = (bool) std::stoul(rootNode->first_node(
				"UseMobilityFactor")->value());
	}
	if (0 != rootNode->first_node("IsIsotropicNetwork")) {
		IsIsotropicNetwork = (bool) std::stoul(rootNode->first_node(
				"IsIsotropicNetwork")->value());
	}
	if (0 != rootNode->first_node("UseTexture")) {
		UseTexture = (bool) std::stoul(
				rootNode->first_node("UseTexture")->value());
	}
	if (0 != rootNode->first_node("ExecuteInParallel")) {
		ExecuteInParallel = (bool) std::stoul(rootNode->first_node(
				"ExecuteInParallel")->value());
	}
	if (0 != rootNode->first_node("MaximumNumberOfThreads")) {
		MaximumNumberOfThreads = std::stoul(rootNode->first_node(
				"MaximumNumberOfThreads")->value());

	}
	if (0 != rootNode->first_node("GridCoarsement")) {
		GridCoarsement = std::stoul(
				rootNode->first_node("GridCoarsement")->value());
	}
	if (0 != rootNode->first_node("GridCoarsementGradient")) {
		GridCoarsementGradient = std::stod(rootNode->first_node(
				"GridCoarsementGradient")->value());
	}
	if (0 != rootNode->first_node("ConvolutionMode")) {
		ConvolutionMode = (E_CONVOLUTION_MODE) std::stoi(rootNode->first_node(
				"ConvolutionMode")->value());
		if (ConvolutionMode >= E_INVALID_VALUE)
			ConvolutionMode = E_INVALID_VALUE;
	}
	//!
	if (0 != rootNode->first_node("ResearchMode")) {
		ResearchMode
				= std::stoul(rootNode->first_node("ResearchMode")->value());
	}
	//! Research modification
	if (0 != rootNode->first_node("ResearchProject")) {
		ResearchProject = (E_RESEARCH_PROJECT) std::stoi(rootNode->first_node(
				"ResearchProject")->value());
		if (ResearchProject >= E_DEFAULT)
			ResearchProject = E_DEFAULT;
	}
	if (0 != rootNode->first_node("ConstantSectorRadius")) {
		ConstantSectorRadius = std::stod(
					rootNode->first_node("ConstantSectorRadius")->value());
	}
	if (0 != rootNode->first_node("InterpolatingSectorRadius")) {
		InterpolatingSectorRadius = std::stod(
						rootNode->first_node("InterpolatingSectorRadius")->value());
		}
	if (0 != rootNode->first_node("NeighbourTracking")) {
		NeighbourTracking = std::stoul(rootNode->first_node("NeighbourTracking")->value());
		}

	if (0 != rootNode->first_node("DislocationEnergy")) {
		DislocationEnergy = std::stoul(rootNode->first_node("DislocationEnergy")->value());
			}
	file.close();
}

#define PUSH_PARAM(param_name) 	\
		temp_string.str("");	\
		temp_string << param_name ;	\
		params->append_node(root->allocate_node(node_element,	\
				root->allocate_string(#param_name),	\
				root->allocate_string(temp_string.str().c_str()) ));

#define PUSH_VALUE(param_name, value) 	\
		temp_string.str("");	\
		temp_string << value ;	\
		params->append_node(root->allocate_node(node_element,	\
				root->allocate_string(#param_name),	\
				root->allocate_string(temp_string.str().c_str()) ));

xml_node<>* Settings::generateXMLParametersNode(xml_document<>* root,
		const char* filename, int loop, int grains) {
	xml_node<>* params = root->allocate_node(node_element, "Parameters", "");
	stringstream temp_string;

	PUSH_VALUE(StartTime, loop);
	PUSH_VALUE(NumberOfParticles, grains);
	PUSH_PARAM(NumberOfPointsPerGrain);
	PUSH_PARAM(AnalysisTimestep);
	PUSH_PARAM(NumberOfTimesteps);
	PUSH_PARAM(DiscreteSamplingRate);
	PUSH_PARAM(DomainBorderSize);
	PUSH_VALUE(MicrostructureGenMode, 0);
	PUSH_VALUE(ReadFromFilename, filename);
	PUSH_PARAM(HAGB);
	PUSH_PARAM(TriplePointDrag);
	PUSH_PARAM(UseMobilityFactor);
	PUSH_PARAM(IsIsotropicNetwork);
	PUSH_PARAM(UseTexture);
	PUSH_PARAM(ExecuteInParallel);
	PUSH_PARAM(MaximumNumberOfThreads);
	PUSH_PARAM(GridCoarsement);
	PUSH_PARAM(GridCoarsementGradient);
	PUSH_PARAM(ConvolutionMode);
	PUSH_PARAM(ResearchMode);
	PUSH_PARAM(ResearchProject);
	PUSH_PARAM(ConstantSectorRadius);
	PUSH_PARAM(InterpolatingSectorRadius);
	PUSH_PARAM(NeighbourTracking);
	PUSH_PARAM(DislocationEnergy);
	return params;
}
#undef PUSH_PARAM
#undef PUSH_VALUE
