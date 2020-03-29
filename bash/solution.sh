export sourceName="mainSolution"

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
	--column "Pick" --column "Problem" FALSE "Sealed_Column" FALSE "Terzaghi" FALSE "Mandel" \
	--width=300 --height=300);

IFS="|"
for word in $ans
do
	if [ $word == "Sealed_Column" ]; then
		problemsSolved+=1
	fi
	if [ $word == "Terzaghi" ]; then
		problemsSolved+=2
	fi
	if [ $word == "Mandel" ]; then
		problemsSolved+=4
	fi
done
IFS=""

cd input
inputOptions=(*.txt)
inputOptions=("${inputOptions[@]%.*}")
cd ..

# medium=$(zenity --list --text "Choose the medium to be simulated" --radiolist --column "Pick" \
# 	--column "Medium" FALSE "${inputOptions[@]}" --width=300 --height=400);

medium="gulfMexicoShale";

# RUN
for ((i=0; i<numRuns; i++));
do
	rm -rf export/*
	cd build
	echo "-- Solving benchmarking problems"
	./$sourceName ${gridType[$i]} ${interpScheme[$i]} ${medium} ${problemsSolved}
	cd ..
	echo "-- Plotting results"
	IFS="|"
	for word in $ans
	do
		if [ $word == "Sealed_Column" ]; then
			python3 -W ignore ./postpro/sealedColumnPlotSolution.py ${medium}
		fi
		if [ $word == "Terzaghi" ]; then
			python3 -W ignore ./postpro/terzaghiPlotSolution.py ${medium}
		fi
		if [ $word == "Mandel" ]; then
			python3 -W ignore ./postpro/mandelPlotSolution.py ${medium}
		fi
	done
	IFS=""
	echo ""
done
rm -rf export/*