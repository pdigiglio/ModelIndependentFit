# Simple bash script to run the fit program after some check on its arguments

#!/bin/bash

fit_executable="runmiFit"
directory="output/"
if [ ! -d "$directory" ]
then
    echo "$directory is not a directory."
    exit 1
fi

model_name="f0_rho0"
root_file="D3PI_${model_name}_RESONANCE_mcmc.root"
if [ ! -f "$directory$root_file" ]
then
    echo "$root_file is not a regular file."
    exit 1
fi

# Call the executable.
./${fit_executable} "${directory}" "$root_file" "$model_name"
