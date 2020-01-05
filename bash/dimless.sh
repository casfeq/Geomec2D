export sourceName="mainDimless"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# DEFINE GRIDS AND INTERPOLATION SCHEMES
declare -a gridType=("staggered" "collocated")
declare -a interpScheme=("NA" "CDS")
numRuns=${#gridType[@]}

# RUN
for ((i=0; i<numRuns; i++));
do
	rm -rf export/*
	cd build
	echo "-- Testing method's dimensionless stability"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "alundum"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "ohioSandstone"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "danianChalk"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "coarseSand"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "hardSediment"
	cd ..
	echo "-- Plotting results"
	python3 -W ignore ./postpro/terzaghiPlotDimless.py "hardSediment" "coarseSand" "danianChalk" \
		"ohioSandstone" "alundum"
	echo ""
done
rm -rf export/*
nautilus plot/
echo ""