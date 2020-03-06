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
./$sourceName "staggered" "NA" "doublePorosity"
./$sourceName "collocated" "CDS" "doublePorosity"
./$sourceName "collocated" "I2DPIS" "doublePorosity"

cd ..
echo "-- Plotting results"
python3 -W ignore ./postpro/doublePorosityPlotStability.py "doublePorosity"
echo ""
rm -rf export/*