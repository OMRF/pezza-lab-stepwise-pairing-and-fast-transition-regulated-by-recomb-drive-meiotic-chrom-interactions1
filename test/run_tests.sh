#!/usr/bin/bash -l 
DEBUG=${DEBUG:=FALSE}

# Capture name of the directory this script is located in and add it to PATH
export script_dir=$(dirname $(realpath $0))
export dir_above=$(dirname $script_dir)
export PATH="$script_dir:$dir_above:$PATH"

# Make R and Python available
module load R      # If needed change this to whatever makes R available to you
module load python # If needed change to whatever makes Python available to you

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

result=$(msd.R input/msd_example.csv X Y Z 2> /dev/null)

test_value "$result" "22.5" "calculated mean square deviation (assuming equal delta times)"

dresser_unspin_algo.py input/unspin/control_fuzzed.csv --center 0 0 0 --radius 10 &> /dev/null
approx_csv input/unspun/control_fuzzed-gpa.csv expected/unspun/control_fuzzed_unspun.csv

dresser_unspin_algo.py input/unspin/control_real.csv --center 0 0 0 --radius 10 &> /dev/null
approx_csv input/unspun/control_real-gpa.csv expected/unspun/control_real_unspun.csv

# Unless DEBUG is on, remove test scripts and test outputs
if [[ "$DEBUG" == "FALSE" ]]; then
    rm test_kymograph.R \
       test_median_detection.R \
       Example.kymograph.pdf \
       change_of_phases.pdf \
       input/unspun/control_fuzzed-gpa.csv \
       input/unspun/control_real-gpa.csv
fi

# Indicate testing is finished
finished_testing
