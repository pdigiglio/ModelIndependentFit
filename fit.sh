# Simple bash script to run the fit program after some check on its arguments

#!/bin/bash

fit_executable="runmiFit"

# This is the .root file name to open
if [ -n "$1" -a -n "$2" ]
then
    if [ -f $1 ]
    then
        ./${fit_executable} $1 $2
    else
        echo "$1 is not a regular file."
    fi
else
    # Run it without arguments to get the
    # help output
    ./${fit_executable}
fi

