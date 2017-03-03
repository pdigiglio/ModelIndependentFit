# Simple bash script to run the fit program after some check on its arguments

#!/bin/bash

# Backups the output directory, if needed.
function backup_output() {
    # Select the name of the most recent backup directory
    local last_backup="$(ls -1d * | grep "${directory}_" | tail -n1)"
    local last_hash=""

    # If there's already a backup.
    if [ -n "${last_backup}" ]
    then
        last_hash="$(cd ${last_backup} && sha1sum * | sha1sum && cd ..)"
    fi

    # The hash of the current output directory.
    local current_hash="$(cd ${directory} && sha1sum * | sha1sum && cd ..)"

    # Backup the files only if the hashes differ.
    if [ "$last_hash" != "$current_hash" ]
    then
        # Make a backup of the output directory
        local backup_directory="${directory}_$(date "+%G-%m-%d_%H:%M:%S")"
        echo " > Copying ${directory} -> ${backup_directory}"
        cp -ar "${directory}" "${backup_directory}"
    else
        echo " > Skip backup as there's nothing new since last backup (${last_backup})."
    fi
}

directory="output"
if [ ! -d "$directory" ]
then
    echo "$directory is not a directory."
    exit 1
fi

model_name="f0_f2"
root_file="${model_name}_mcmc.root"
#root_file="model_mcmc.root"
if [ ! -f "$directory/$root_file" ]
then
    echo "$root_file is not a regular file."
    exit 1
fi

# Call the executable.
fit_executable="runmiFit"
if [ ! -x "$fit_executable" ]
then
    echo "$fit_executable is not an executable file."
fi

# Check if the debug option is active
if [ -z "$1" ]
then
    ./${fit_executable} "${directory}" "${root_file}" "${model_name}"
elif [ "$1" == "--debug" ]
then
    gdb --args ./${fit_executable} "${directory}" "${root_file}" "${model_name}"
elif [ "$1" == "--valgrind" ]
then
    valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all ./${fit_executable} "${directory}" "${root_file}" "${model_name}"
else
    echo " > Unknown option '$1'"
fi

backup_output
