clear
clear
clear
mkdir -p export
mkdir -p plot
mkdir -p build
rm -rf plot/*

# PARAMETERS
declare -a gridType=()
declare -a interpScheme=()
declare -a problemsSolved=0
declare -a inputOptions
declare -a medium

gridType+=("staggered")
interpScheme+=("NA")
# gridType+=("collocated")
# interpScheme+=("CDS")
# gridType+=("collocated")
# interpScheme+=("I2DPIS")
numRuns=${#gridType[@]}
# problemsSolved+=2
# problemsSolved+=4
problemsSolved+=8

# VERIFICATION
rm -rf export/*
export sourceName="mainDoubleExact"
medium="gulfMexicoShale";

# COMPILE
cd build
cmake ..
make
cd ..
echo ""


# RUN
echo "-- Solving benchmarking problems"
for ((i=0; i<numRuns; i++));
do
	cd build
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} ${medium} ${problemsSolved} "0.97" "0.03"
	echo ""
	cd ..
done

# PLOT
echo "-- Plotting results"
# python3 -W ignore ./postpro/sealedDoublePlotSolution.py ${medium}
# python3 -W ignore ./postpro/storageDoublePlotSolution.py ${medium}
python3 -W ignore ./postpro/leakingDoublePlotSolution.py ${medium}
echo ""