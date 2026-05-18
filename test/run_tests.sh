#!/usr/bin/bash -l 

# Capture name of the directory this script is located in and add it to PATH
export script_dir=$(dirname $(realpath $0))
export dir_above=$(dirname $script_dir)
export PATH="$script_dir:$dir_above:$PATH"

# Make R language environment available 
module load R

# Add test functions
source test_functions.sh

# Make copy of scripts with command line input statements uncommented and chooser script commented
sed '16,17 s/^#//' $(which kymograph.R) | sed '20 s/^/#/' > test_kymograph.R

# Make copy of median detection script with command line input statements uncommented
# (and comment the chooser method)
sed '23,24 s/^#//' $(which median_detection.R) | sed '27 s/^/#/' > test_median_detection.R

# make test scripts executable
chmod +x test_kymograph.R test_median_detection.R

# Run and test scripts
test_kymograph.R input/Example.xlsx &> /dev/null # ignore stdout/stderr by sending to /dev/null
test_pdf Example.kymograph.pdf expected/Example.kymograph.pdf

test_median_detection.R input/Example.xlsx &> /dev/null # ignore stdout/stderr by sending to /dev/null
test_pdf change_of_phases.pdf expected/change_of_phases.pdf


# Remove test scripts and test outputs
rm test_kymograph.R test_median_detection.R
rm Example.kymograph.pdf change_of_phases.pdf

result=$( msd.R input/msd_example.csv X Y Z)

test_value "$result" "22.5" "calculated mean square deviation (assuming equal delta times)"

# Indicate testing is finished
finished_testing
