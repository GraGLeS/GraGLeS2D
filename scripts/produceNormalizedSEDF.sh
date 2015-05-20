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
### functions block

#The original energy files are altered in the second column. The absolute
#length is substituted by the corresponding bin height. For each file 
# a temporary file with the conversions is returned. Requires two arguments:
# filename and discreteSamplingRate.
normalizeLength(){
	if [ $# -eq 0 ] ; then
		echo "Pass a filename and a sampling rate as parameters to this function."
		exit 1
	fi

	nL_fileName=$1
	nL_samplingRate=$2
	nL_timestep=`echo ${nL_fileName//[A-Z,a-z,.,_,\/]}`
	
	#Calculate the sum for normalizing to relative frequencies
	nL_area=`awk 'BEGIN {sum = 0}
	{sum = sum + $2}
	END {print sum}' $nL_fileName`
	
	#Generate height of the rectangle as normalized relative frequency divided by step size , i.d. width of bin.
	#This ensues that the cumulative area of all rectangles equals 1.
	awk -v nL_area="$nL_area" -v nL_samplingRate="$nL_samplingRate" '{$2=($2/(nL_area))/nL_samplingRate ; print $1,$2 }' $nL_fileName >> tempEnergyFile${nL_timestep}.txt 
	return 0

}
# The column with the lengths are scrutinized and the maimal value is returned.
# A parameter which is constituted of the filename as a string is requisite.
returnMaxLength(){
	if [ $# -eq 0 ] ; then
		echo "Pass a filename as a parameter to this function."
		exit 1
	fi

	rml_fileName=$1
	
	rml_value=`awk 'BEGIN {max = 0}
		{max = $2>max ? $2 : max} 
		END {print max}' $rml_fileName`
	echo "$rml_value"
	return 0
}

# The lines of the input file are counted and returned.
# A parameter specifying the filename is requisite.
numberOfEnergyClasses(){
	if [ $# -eq 0 ] ; then
		echo "Pass a filename as a parameter to this function."
		exit 1
	fi

	noec_fileName=$1

	noec_amountLines=`wc -l $noec_fileName`
	# wc -l returns not only a number but also the filename. So we need to truncate the last output.
	noec_amountLinesCutted=${noec_amountLines% *}
	echo "$noec_amountLinesCutted"
	return 0
}


# This function cuts the first word of the last line of the input file to obtain the biggest energy class.
# A parameter specifying the filename is requisite.
maximumEnergyClass(){
	if [ $# -eq 0 ] ; then
		echo "Pass a filename as a parameter to this function."
		exit 1
	fi

	mec_lastLine=`tail -n 1 $1`
	mec_firstWord=${mec_lastLine% *} # There is a tab character in the original SEDF scipt
	echo "$mec_firstWord"
	return 0
}

# The gnuplot config file is created. This function expects four  parameters:
# the number of energy classes, the maximum energy class, the maximum length, i.e. proportion
# and the output file name as well as a counter for outputfiles.

generateGnuplotFile(){
	
	if [ $# -ne 5 ] ; then 
		echo "You should pass five arguments."
		exit 1
	fi
	
	ggf_numberOfClasses=$1
	ggf_maximumEnergyClass=$2
	ggf_maximumLength=$3
	ggf_fileName=${4%.txt}
	ggf_counter=$5

	# The following paragraph generates the gnu commands and stores them in a file.
	echo "# This test file executes gnu commands. In particular 
	# data is read from a file and used to produce a gnu 
	# output
	
	set title \"Energy distribution histogram\"
	
	# In some cases the following command is necessary
	# reset
	
	n = $ggf_numberOfClasses # Depends on the \"DiscreteSamplingRate\"
	min = 0 
	max = $ggf_maximumEnergyClass  #Equals the \"HAGB\"-value
	
	# Max divided through n yields the rectangle width
	width = max / n
	
	# The function hist employs the variable width to assign a value x to an intervall
	# hist(x,width) = width * floor(x/width) + width/2.0 # mybe just adding 1 is appropriate
		
	# Provi des the option to generate an output file
	set term png giant size 1024,768
	set output \"histogram_$ggf_counter.png\"
	
	set yrange[0:$ggf_maximumLength]
	set xrange[min:max]
	
	# The following commands pimp the output apperance
	set offset graph 0.0,0.025,$ggf_maximumLength*0.05,0.0
	set xtics min, 0.1 , max
	set boxwidth width*0.85
	set style fill solid 0.5
	set tics out nomirror

	# set palette rgbformulae 33,13,10 
	# set cbrange[0:0.6] 
	# set cbtics 0.1
	
	set ylabel \"Proportion id est lengths\"
	set xlabel \"Energy intervalls\"
	
	plot \"${ggf_fileName}.txt\" using 1:2  with boxes lc rgb\"blue\" title \"$ggf_counter timestep\"" > ${ggf_fileName}.plt

	return 0	

}

### --------------------------------------------------------------- ###
### main 

# Test if at leat one energy file exists

#testerFile=EnergyLengthDistribution_*
#testerfile=${testerFile%*}
#echo $testerfile
#if [ -f ] ; then 


#fi

sortedEnergyFiles=`find ./ -type f -name "EnergyLengthDistribution_*" | sort -k2 -t_ -n`

for energyFile in $sortedEnergyFiles
do
	
	firstLine=`head -1 $energyFile`
	firstLineT=${firstLine%	*}
	secondLine=`tail +2 $energyFile | head -1`
	secondLineT=${secondLine%	*}	
	difference=$(echo "$secondLineT -  $firstLineT" | bc)
	normalizeLength $energyFile $difference
done

# Determine maximal y-value of all timesteps
	
globalMaxY=0
sortedEnergyFiles=`find ./ -type f -name "tempEnergyFile*" | sort -k2 -t_ -n`

for energyFile in $sortedEnergyFiles 
do
	
	currentMaxY=`returnMaxLength $energyFile`
	#comparator=`expr $currentMaxY '>' $globalMaxY`	
	comparator=$(echo "$currentMaxY > $globalMaxY" | bc)

	if [ $comparator -eq 1 ] ; then
	       globalMaxY=$currentMaxY
	fi	       
quantityEnergyClasses=`numberOfEnergyClasses $energyFile`
topEnergyClass=`maximumEnergyClass $energyFile`
#echo "$topEnergyClass"
done

# echo "$globalMaxY" ; echo "$quantityEnergyClasses" ; echo "$topEnergyClass"

# Generate .pngs with gnuplot
#counter=1
currentTimestep=0
for energyFile in $sortedEnergyFiles 
do
	currentTimestep=`echo ${energyFile//[A-Z,a-z,.,_,\/]}`
	#echo "$quantityEnergyClasses $topEnergyClass $globalMaxY $energyFile $currentTimestep"
	generateGnuplotFile $quantityEnergyClasses $topEnergyClass $globalMaxY $energyFile $currentTimestep
	#echo "$energyFile"
	gnuplot ${energyFile%.txt}.plt
	rm ${energyFile%.txt}.plt

	#counter=$((counter + 1))
done

# Convert pngs to one SEDF.gif

histogramList=`find ./ -type f -name "histogram_*" | sort -k2 -t_ -n`
# echo $histogramList
convert -delay 40 -loop 0 $histogramList SEDFNormalized.gif
rm histogram*
rm tempEnergyFile*
