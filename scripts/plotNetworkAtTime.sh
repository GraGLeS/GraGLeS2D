#
#	GraGLeS 2D A grain growth simulation utilizing level set approaches
#   Copyright (C) 2015  Christian Miessen, Nikola Velinov
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
FIRST=0
if [ $# -lt 1 ]; then
	echo "Timestep not specified!"
	echo "Usage: plotNetworkAtTime.sh <timestep>"
	exit 2
fi

COMMAND="gnuplot -persist -e \"set palette rgbformulae 33,13,10; set view map; set cbrange[0:1]; set cbtics 0.1; set yrange [:] reverse; plot './Network_Timestep_$1.gnu' w l palette title 'Timestep $1'\""
echo $COMMAND
eval $COMMAND
