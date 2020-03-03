mkdir -p export
mkdir -p plot
mkdir -p build
rm -rf plot/*
rm -rf export/*

export sourceName="mainSolution"

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
./$sourceName "staggered" "NA" "gulfMexicoShale" 8
cd ..
# echo "-- Plotting results"
# python3 -W ignore ./postpro/terzaghiPlotSolution.py "gulfMexicoShale"
# echo ""
# rm -rf export/*