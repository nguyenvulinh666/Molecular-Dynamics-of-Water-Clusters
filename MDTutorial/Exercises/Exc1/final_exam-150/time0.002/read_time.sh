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

# Process the input file
grep "timestep" "$input_file" | while read line; do
  # Extract the timestep value (the second number)
  timestep=$(echo "$line" | awk '{print $2}')
  
  #echo $timestep
  #echo $multiplier  

  # Calculate the time
  time=$(echo "$timestep * $multiplier" | bc)

  #append to output
  echo "$time" >> "$output_file"

  #Increment the time variable by the increment.
 # time=$(echo "$time + $increment" | bc)

done
#if ! test -d "/path/to/directory"; then
#	mkdir oxygen_data
#fi

for i in {1..373..3}
do
	> "oxygen_data/O$i.txt"
done

for i in {1..373..3}
do
	awk -v var="$i" '/^OW[[:space:]]+'${i}' .*/ {getline; print}' $input_file >> "oxygen_data/O$i.txt"
done

> combine.txt

paste -d " " ./oxygen_data/*.txt >> combine.txt

#paste -d " " O1.txt O4.txt O7.txt O10.txt O13.txt O16.txt >> combineO1O16.txt

#awk '/^OW[[:space:]]+4 .*/ {getline; print}' $input_file >> "O2.txt"
echo "Time step data has been written to $output_file"
