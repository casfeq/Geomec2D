export sourceName="mainSolution"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# DEFINE GRIDS AND INTERPOLATION SCHEMES
declare -a gridType=("staggered" "collocated" "collocated" "collocated" "collocated")
declare -a interpScheme=("NA" "CDS" "1DPIS" "I2DPIS" "C2DPIS")
numRuns=${#gridType[@]}


# RUN
for ((i=0; i<numRuns; i++));
do
	rm -rf export/*
	cd build
	echo "-- Solving benchmarking problems"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "gulfMexicoShale"
	cd ..
	echo "-- Plotting results"
	python3 -W ignore ./postpro/sealedColumnPlotSolution.py "gulfMexicoShale"
	python3 -W ignore ./postpro/terzaghiPlotSolution.py "gulfMexicoShale"
	python3 -W ignore ./postpro/mandelPlotSolution.py "gulfMexicoShale"
	echo ""
done
rm -rf export/*
nautilus plot/
echo ""