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
./$sourceName "staggered" "NA" "softSediment"
./$sourceName "collocated" "CDS" "softSediment"
./$sourceName "collocated" "1DPIS" "softSediment"
./$sourceName "collocated" "I2DPIS" "softSediment"
./$sourceName "collocated" "C2DPIS" "softSediment"

cd ..
echo "-- Plotting results"
python3 -W ignore ./postpro/doublePorosityPlotStability.py "softSediment"
echo ""
# rm -rf export/*