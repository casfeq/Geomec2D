export sourceName="mainStability"

# COMPILE
cd build
cmake ..
make
cd ..
echo ""

# GUI
declare -a gridType=()
declare -a interpScheme=()
declare -a problemsSolved=0
declare -a inputOptions
declare -a medium
ans=$(zenity --list --text "Choose grid type and interpolation schemes" --checklist \
	--column "Pick" --column "Formulation" FALSE "Staggered_Grid" FALSE "Collocated_Grid_with_CDS" \
	FALSE "Collocated_Grid_with_1D_PIS" FALSE "Collocated_Grid_with_I_2D_PIS" \
	FALSE "Collocated_Grid_with_C_2D_PIS" --width=300 --height=300);

IFS="|"
for word in $ans
do
	if [ $word == "Staggered_Grid" ]; then
		gridType+=("staggered")
		interpScheme+=("NA")
	elif [ $word == "Collocated_Grid_with_CDS" ]; then
		gridType+=("collocated")
		interpScheme+=("CDS")
	elif [ $word == "Collocated_Grid_with_1D_PIS" ]; then
		gridType+=("collocated")
		interpScheme+=("1DPIS")
	elif [ $word == "Collocated_Grid_with_I_2D_PIS" ]; then
		gridType+=("collocated")
		interpScheme+=("I2DPIS")
	elif [ $word == "Collocated_Grid_with_C_2D_PIS" ]; then
		gridType+=("collocated")
		interpScheme+=("C2DPIS")
	fi
done
IFS=""
numRuns=${#gridType[@]}

ans=$(zenity --list --text "Choose the benchmarking problems to be solved" --checklist \
	--column "Pick" --column "Problem" FALSE "Terzaghi" FALSE "Mandel" --width=300 --height=300);

IFS="|"
for word in $ans
do
	if [ $word == "Terzaghi" ]; then
		problemsSolved+=2
	fi
	if [ $word == "Mandel" ]; then
		problemsSolved+=4
	fi
done
IFS=""

# RUN
for ((i=0; i<numRuns; i++));
do
	rm -rf export/*
	cd build
	echo "-- Testing method's stability"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} "softSediment"  ${problemsSolved}
	cd ..
	echo "-- Plotting results"
	IFS="|"
	for word in $ans
	do
		if [ $word == "Terzaghi" ]; then
			python3 -W ignore ./postpro/terzaghiPlotStability.py "softSediment"
		fi
		if [ $word == "Mandel" ]; then
			python3 -W ignore ./postpro/mandelPlotStability.py "softSediment"
		fi
	done
	IFS=""
	echo ""
done
# echo "-- Plotting results"
# python3 -W ignore ./postpro/terzaghiPlotStabilityComparison.py "softSediment"
# python3 -W ignore ./postpro/mandelPlotStabilityComparison.py "softSediment"
rm -rf export/*
