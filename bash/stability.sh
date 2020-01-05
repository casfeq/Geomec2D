export sourceName="mainStability"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# DEFINE GRIDS AND INTERPOLATION SCHEMES
declare -a gridType=("staggered" "collocated" "collocated" "collocated" "collocated")
declare -a interpScheme=("NA" "CDS" "1DPIS" "I2DPIS" "C2DPIS" )
numRuns=${#gridType[@]}

# RUN
for ((i=0; i<numRuns; i++));
do
	rm -rf export/*
	cd build
	echo "-- Testing method's stability"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "softSediment"
	cd ..
	echo "-- Plotting results"
	python3 -W ignore ./postpro/terzaghiPlotStability.py "softSediment"
	python3 -W ignore ./postpro/mandelPlotStability.py "softSediment"
	echo ""
done
# echo "-- Plotting results"
# python3 -W ignore ./postpro/terzaghiPlotStabilityComparison.py "softSediment"
# python3 -W ignore ./postpro/mandelPlotStabilityComparison.py "softSediment"
rm -rf export/*
nautilus plot/
echo ""