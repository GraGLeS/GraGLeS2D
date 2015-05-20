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
#!/bin/sh


### --------------------------------------------------------------- ###
### Usage text
USAGE='
commandUnit.sh
Generates a area variation gif at a certain timestep
--------------------------------

Options
-bg	Include bouandary grains in analysis
-file	Produce output files, i.e. a png file as well as a gnu-plt file.

-h	Show usage text

Arguments: 
<t>	Choose the timestep to be analysed.
<ats>	Give the used analyse timestep to coincide with the requested timestep.
comment=<text>	Pass comment phrase.
'

### --------------------------------------------------------------- ###
### Takes care of the parameter input

# Options
if ([ $# -ne 2 ] && [ $# -ne 3 ] && [ $# -ne 4 ] && [ $# -ne 5 ]); then
	echo "Two arguments are required. The first argument is the time and the second is the employed timestep in the analysis."
	exit 1
fi

while [[ "$1" == -* ]]
do
	[[ "$1" == -bg  ]] && opt_bg=yes && shift && continue
	[[ "$1" == -file  ]] && opt_file=yes && shift && continue

	# Additional options
	[[ "$1" == *h* ]] && echo "$USAGE" && exit
	
	shift

done

# Arguments
 for arg in $*
 do
	[[ "$arg" == *=* ]] && eval $arg

 done


### --------------------------------------------------------------- ###
### main 

# Area variation for one timestep and all grains

# Cleanup (always run without bouandary grains before with boundary grains if you want to save files)
if [ "$opt_bg" ] ; then

	pattern=`echo ./areaVariation_gnuplotFile_*_withBG_*.plt`
else
	pattern=`echo ./areaVariation_gnuplotFile_*.plt`
fi

pattern=${pattern% *}

if ([ -f "$pattern" ] && [ "$opt_bg" ]) ; then
	rm ./areaVariation_gnuplotFile_*_withBG_*.plt
	#ls
fi

if [ -f "$pattern" ] ; then
	rm ./areaVariation_gnuplotFile_*.plt
	#ls
fi

# Pass the arguments and test if the parameterization suffices the minium timestep
atTime=$1
analysisTimestep=$2

if [ $atTime -lt $analysisTimestep ] ; then
	echo "Time 0 does not exist."
	exit 1	
fi

currTimestep=`expr $atTime / $analysisTimestep`
echo "The analysed timestep is: $currTimestep"	
for datei in `find ./ -type f -name "Area*"`
do
	sed -n "${currTimestep}p" $datei >> gnuplotFileTemp_${atTime}.plt
done	
if [ "$opt_bg" ] ; then

	awk 'BEGIN {counter=1}
		{$1=counter ; print $1,$2,$3,$4 ; counter++}' gnuplotFileTemp_${atTime}.plt >> areaVariation_gnuplotFile_${atTime}_withBG_${comment}.plt
else

	awk 'BEGIN {counter=1} 
		$4<1 {$1=counter ; print $1,$2,$3,$4 ; counter++}' gnuplotFileTemp_${atTime}.plt >> areaVariation_gnuplotFile_${atTime}_${comment}.plt
fi
rm gnuplotFileTemp_${atTime}.plt

# Switch to plot-file with boundary grains
if [ "$opt_bg" ] ; then
	plotString=areaVariation_gnuplotFile_${atTime}_withBG_${comment}.plt
	titleString="Area variation at t=${atTime} with boundary grains; $comment"
else
	plotString=areaVariation_gnuplotFile_${atTime}_${comment}.plt
	titleString="Area variation at t=${atTime}; $comment"
fi
if [ "$opt_file" ] ; then
	
	if [ "$opt_bg" ] ; then 
		outputFileName=areaVariation_${atTime}_withBG_${comment}.png
	else
		outputFileName=areaVariation_${atTime}_${comment}.png
	fi

	outputString="set term png giant size 1024,768
	set output \"$outputFileName\""
fi

gnuplot -persist << END_FILE
set xrange[0:]
set yrange[-15:15]
set ytics -15,2.5,15	
set offset graph 0.05,0.05,0.05,0.0
# Provides the option to generate an output file
$outputString
plot "$plotString" using 1:2 title "$titleString"
END_FILE


