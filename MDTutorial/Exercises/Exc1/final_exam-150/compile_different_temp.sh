timesteps=(0.001 0.002 0.0005 0.0008)
time=20

for ts in "${timesteps[@]}"; do
	step=$(echo "scale=0; $time / $ts" | bc) # scale is the decimal point

	mkdir "time$ts"
	
	# Changing different time step inside the CONTROL file
	sed -i "s/^\(timestep[[:space:]]*\).*/\1$ts/" CONTROL
	sed -i "s/^\(steps[[:space:]]*\).*/\1i $step/" CONTROL

	# Copy from the mail folder to other folder the necessity file for simulation
	cp CONTROL time$ts/CONTROL
	cp CONFIG time$ts/CONFIG
	cp FIELD time$ts/FIELD
	cp read_time.sh time$ts/read_time.sh
	cp compile_data.sh time$ts/compile_data.sh
	cd time$ts
	
	# Run DLPOLY.X module
	./../../../../execute/DLPOLY.X 
	
	# Compile to txt file
	./compile_data.sh
	cd ..
done
