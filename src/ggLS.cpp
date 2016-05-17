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
#include "ggLS.h"
#include "Settings.h"
#include <sys/time.h>
using namespace voro;
using namespace std;

int main(int argc, char *argv[]) {

	if (argc > 1)
		Settings::initializeParameters(argv[1]);
	else
		Settings::initializeParameters();

	grainhdl* my_sim = NULL;

	my_sim = new grainhdl();

	my_sim->setResearchAdjustments(Settings::ResearchProject);
	timeval time_start, time_end;
	gettimeofday(&time_start, NULL);
	my_sim->setSimulationParameter();
	gettimeofday(&time_end, NULL);
	double elapsed_secs = (time_end.tv_sec - time_start.tv_sec)
			+ (time_end.tv_usec - time_start.tv_usec) / 1000000.0;
	cout << "elapsed secs for Initializing network (Read):" << elapsed_secs << endl << endl;

	if (Settings::MicrostructureGenMode == E_GENERATE_WITH_VORONOY)
		my_sim->save_Full_Microstructure_for_Restart();

	gettimeofday(&time_start, NULL);
	my_sim->run_sim();

	gettimeofday(&time_end, NULL);
	elapsed_secs = (time_end.tv_sec - time_start.tv_sec) + (time_end.tv_usec
			- time_start.tv_usec) / 1000000.0;
	cout << "elapsed secs for main loop:" << elapsed_secs << endl;

	my_sim->save_NrGrainsStats();

	my_sim->clear_mem();

	delete my_sim;

}
