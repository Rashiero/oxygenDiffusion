#!/bin/bash

#####
# This file runs all the routines required to process a microscopy image of vasculature and simulate oxygen perfusion
#####

# filename="P8_Segmented"
# filename="P9_Segmented"
# filename="P10_Segmented"
# filename="P11_Segmented"
# filename="P12_Segmented"

mkdir "${filename}_Run"

echo Processing Images

# Process image from /Data/$filename
python3 image_processing.py $filename
# Outputs $filename_skeleton.tif

echo Processing done - Generating network from skeleton

# Processes image from $filename_skeleton.tif
python3 network_generation.py $filename
# Outputs $filename_primary_network.json

echo Processing network topology

# Process $filename_primary_network.json
python3 network_cleanup.py $filename
# Outputs $filename_final_network.json

# Process $filename_final_network.json
python3 calculate_diameter.py $filename
# Outputs $filename_final_network.json updated with diameters

#Modify parameters file
sed -i "/${filename}/s/^# //g" parameters.py

echo Running Simulation

# Process $filename_final_network.json
python3 main.py $filename
# Outputs $filename_Network.obj and $filename_Tissue_Oxygen_Pressure.vti

sed -i "/${filename}/s/^Nx/# Nx/g" parameters.py
echo Restored parameter file

mv "${filename}_skeleton.tiff" "${filename}_primary_network.json" "${filename}_final_network.json" "${filename}_Network.obj" "${filename}_Tissue_Oxygen_Pressure.vti" "${filename}_Run"
echo Moved all files to ${filename}_Run