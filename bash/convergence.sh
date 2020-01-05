export sourceName="mainConvergence"

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
    echo "-- Testing method's convergence"
    ./$sourceName ${gridType[$i]} ${interpScheme[$i]} "gulfMexicoShale"
    cd ..
    echo "-- Plotting results"
    python3 -W ignore ./postpro/terzaghiPlotConvergence.py
    echo ""
done
rm -rf export/*
nautilus plot/
echo ""