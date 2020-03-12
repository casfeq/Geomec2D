export sourceName="mainDoubleHStability"

# COMPILE
cd ..
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
# ./$sourceName "collocated" "CDS" "doublePorosity"
# ./$sourceName "collocated" "I2DPIS" "doublePorosity"

# # PLOT
# cd ..
# echo "-- Plotting results"
# python3 -W ignore ./postpro/doublePorosityPlotHStability.py "doublePorosity"
# echo ""
# rm -rf export/*

# export sourceName="mainDoubleTStability"

# # COMPILE
# cd build
# cmake ..
# make
# cd ..
# echo ""

# # RUN
# rm -rf export/*
# cd build
# echo "-- Solving benchmarking problems"
# ./$sourceName "staggered" "NA" "doublePorosity"
# ./$sourceName "collocated" "CDS" "doublePorosity"
# ./$sourceName "collocated" "I2DPIS" "doublePorosity"

# # PLOT
# cd ..
# echo "-- Plotting results"
# python3 -W ignore ./postpro/doublePorosityPlotTStability.py "doublePorosity"
# echo ""
# rm -rf export/*