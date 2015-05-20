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
if [ $# -lt 1 ]; then
	echo "Specify timestep size!"
	echo "Usage: produceAnimation_aux.sh <timestep_size>"
	exit 2
fi
#Dunno how to negate this ...
if [[ $1 =~ ^-?[0-9]+$ ]]; then
	echo
else
	echo "Specify correct timestep size!"
	echo "Usage: produceAnimation_aux.sh  <timestep_size>"
	exit 2
fi
PATTERN=./Network_Timestep_
TIMESTEP=0
while true; do
	CURRENT_FILENAME=$PATTERN`printf "%06d" $TIMESTEP`.png
	echo -e -n "Processing timestemp $TIMESTEP\r"
	CURRENT_PATTERN=$PATTERN$TIMESTEP
	CURRENT_PATTERN=$CURRENT_PATTERN".gnu"
	FIRST=""
	if [ -f  $CURRENT_PATTERN ]; then			
		GNUPLOT_STRING="gnuplot -e \"set palette rgbformulae 33,13,10; set cbrange[0:1.0]; set cbtics 0.1; set term png giant size 1024,768; set output '$CURRENT_FILENAME'; set yrange [:] reverse; plot"
		GNUPLOT_STRING="$GNUPLOT_STRING '$CURRENT_PATTERN' w l palette title 'Timestep $TIMESTEP' \""
		FIRST="done"
	eval $GNUPLOT_STRING
	else
		break;
	fi
	TIMESTEP=`expr $TIMESTEP + $1`
done
echo ""
echo "Done creating intermidiate files."
echo "Creating gif with resolution 1024 , 768"
convert -delay 40 -loop 0 Network_Timestep_*.png NetworkEvolution.gif
echo "Done creating gif. Cleaning up.."
rm *.png
echo "Done."

