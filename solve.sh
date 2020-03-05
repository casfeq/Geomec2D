clear
clear
clear
mkdir -p export
mkdir -p plot
mkdir -p build
rm -rf plot/*
rm -rf export/*

export sourceName="mainDouble"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# RUN
rm -rf export/*
cd build
echo "-- Solving benchmarking problems"
./$sourceName "staggered" "NA" "hardSediment"
./$sourceName "collocated" "CDS" "hardSediment"
./$sourceName "collocated" "1DPIS" "hardSediment"
./$sourceName "collocated" "I2DPIS" "hardSediment"
./$sourceName "collocated" "C2DPIS" "hardSediment"

cd ..
echo "-- Plotting results"
python3 -W ignore ./postpro/doublePorosityPlotStability.py "hardSediment"
echo ""
rm -rf export/*