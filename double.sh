clear
clear
clear
mkdir -p export
mkdir -p plot
mkdir -p build
rm -rf plot/*

# VERIFICATION
rm -rf export/*
export sourceName="mainDoubleSolution"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# PARAMETERS
declare -a gridType=()
declare -a interpScheme=()
declare -a problemsSolved=0
declare -a inputOptions
declare -a medium

gridType+=("staggered")
interpScheme+=("NA")
gridType+=("collocated")
interpScheme+=("CDS")
gridType+=("collocated")
interpScheme+=("I2DPIS")
numRuns=${#gridType[@]}
problemsSolved+=2
medium="gulfMexicoShale";

# SOLVE
mkdir export/lastrun
echo "-- Solving benchmarking problems"
for ((i=0; i<numRuns; i++));
do
	cd build
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} ${medium} ${problemsSolved} "0" "1"
	echo ""
	cd ..
done
mv export/terzaghi_${medium}_MacroPNumeric_* export/lastrun
for ((i=0; i<numRuns; i++));
do
	cd build
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} ${medium} ${problemsSolved} "1" "0"
	echo ""
	cd ..
done
mv export/lastrun/* export

# PLOT
echo "-- Plotting results"
python3 -W ignore ./postpro/doublePlotSolutionCILAMCE2020.py ${medium}
echo ""

# # STABILITY
# rm -rf export/*
# export sourceName="mainStability"

# # COMPILE
# cd build
# cmake ..
# make
# cd ..
# echo ""

# # PARAMETERS
# problemsSolved+=8
# medium="modifiedHardSediment";

# # RUN
# echo "-- Testing method's stability"
# echo ""
# for ((i=0; i<numRuns; i++));
# do
# 	cd build
# 	./$sourceName ${gridType[$i]} ${interpScheme[$i]} ${medium} ${problemsSolved}
# 	echo ""
# 	cd ..
# done

# # PLOT
# echo "-- Plotting results"
# # python3 -W ignore ./postpro/terzaghiPlotStabilityPaper.py ${medium}
# # python3 -W ignore ./postpro/mandelPlotStabilityPaper.py ${medium}
# python3 -W ignore ./postpro/stripfootPlotStabilityPaper.py ${medium}
# echo ""
rm -rf export/*