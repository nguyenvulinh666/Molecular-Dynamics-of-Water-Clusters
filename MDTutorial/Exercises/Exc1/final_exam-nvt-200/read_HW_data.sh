#!/bin/bash

# Input file (replace with your file name)
input_file="HISTORY"

# Output file
output_file="timestep_data.txt"

> $output_file
# Time step multiplier
#multiplier=100

# Time step increment (from the last number)
multiplier=0.001000

# Initialize time variable
time=0

#if ! test -d "/path/to/directory"; then
mkdir ./HW-HW
#fi

for i in {1..58..3}
do
	> "HW-HW/H$((i+1)).txt"	
	> "HW-HW/H$((i+2)).txt"	
done

for i in {1..58..3}
do
	#awk -v var="$i" '/^OW[[:space:]]+'${i}' .*/ {getline; print}' $input_file >> "oxygen_data/O$i.txt"
	awk -v var="$i" '/^HW[[:space:]]+'$((i+1))' .*/ {getline; print}' $input_file >>  "./HW-HW/H$((i+1)).txt"
	awk -v var="$i" '/^HW[[:space:]]+'$((i+2))' .*/ {getline; print}' $input_file >>  "./HW-HW/H$((i+2)).txt"
done

> displacement-HW.txt

paste -d " " ./HW-HW/*.txt >> displacement-HW.txt

#paste -d " " O1.txt O4.txt O7.txt O10.txt O13.txt O16.txt >> combineO1O16.txt

#awk '/^OW[[:space:]]+4 .*/ {getline; print}' $input_file >> "O2.txt"
echo "Time step data has been written to $output_file"
