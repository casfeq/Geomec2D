mkdir -p export
mkdir -p plot
mkdir -p build
rm -rf plot/*
rm -rf export/*

# GUI
ans=$(zenity --list --text "Choose an application" --radiolist --column "Pick" --column \
	"Application" FALSE "Verification" FALSE "h_Convergence_Analysis" FALSE \
	"t_Convergence_Analysis" FALSE "Stability_Test" FALSE "Dimensionless_Stability_Test" \
	--width=300 --height=300); echo $ans
if [ $ans == "Verification" ]; then
	clear
	./bash/solution.sh
elif [ $ans == "h_Convergence_Analysis" ]; then
	clear
	./bash/convergence.sh
elif [ $ans == "t_Convergence_Analysis" ]; then
	clear
	./bash/tconvergence.sh
elif [ $ans == "Stability_Test" ]; then
	clear
	./bash/stability.sh
elif [ $ans == "Dimensionless_Stability_Test" ]; then
	clear
	./bash/dimless.sh
else
	clear
	zenity --error --text "No application selected "
fi