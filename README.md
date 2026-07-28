# Stepwise pairing and fast transition regulated by recombination drive meiotic chromosome interactions

# Scripts

## kymograph.R

### USAGE:
    kymograph.R "<folder>/<file name>.xlsx"

### Eample:
    kymograph.R test/input/Example.xlsx

Input XLSX file must contain a sheet named "Distance To Nearest Neighbour".
Column headers must be in the second row.
One column header must be "Distance To Nearest Neighbour".
One column header must be "Time".
Other columns are ignored.
Rows with duplicate Time values are ignored (we expect them to correspond to completely duplicate rows).

## median_detection.R test/input/Example.xlsx

Input file (e.g., `Example.xlsx`) must contain a sheet named "Distance To Nearest Neighbour" with column headers in the second row.
One column header must be "Distance To Nearest Neighbour" and the other "Time".
Rows with duplicate Time values are ignored (we expect them to correspond to completely duplicate rows).

## msd.R

### USAGE:
    msd.R "<folder>/<file name>.csv"

### Example:
    msd.R test/input/msd_example.csv 

The CSV file contains the columns x, y, and z, separated by commas. Other columns are ignored..  

Calculate mean square distance from each data point to the first data point.  

See https://en.wikipedia.org/wiki/Mean_squared_displacement#Derivation_for_n_dimensions.

## dresser_unspin_algo.py

### USAGE:
    dresser_unspin_algo.py "<folder>/<file name>.csv" --center CX CY CZ --radius R

### Example:
    dresser_unspin_algo.py test/input/unspin/control_fuzzed.csv  --center 0 0 0 --radius 10

The CSV file contains the columns "labels,x,y,z,t" representing the the spot's label, x,y,z coordinates and timepoint, respectively. Spot labels end in "<spot number>H" or "<spot number>G". The nuclear center and radius are measured at timepoint 0 (TP0) and given as command-line arguments.  The Hoechst (H) spots used for registration are either specified by their spot numbers or all H spots are used when the cell has exactly three.  

If you're using the uv environment for this script, you can also proceed those commands with `uv run`.  

### Command line flags

`--center CX CY CZ`  (required) measured nuclear center at timepoint 0.  
`--radius R`  (required) measured nuclear radius (assumed stable over the movie).  
`--out-dir DIR`  output directory (default: a sister "unspun" folder beside the input, for example "movie_name/unspun/" where the input data are in "movie_name/raw/").  
`--time-step-seconds S`  seconds per time-point, used for "avg_speed" (default 60.0; may be given as "60").  
`--proj-to-periph`  also project the non-registration (G) spots onto the sphere.  Is off by default so that interior G spots keep their real radial position.  
`--H-spots [A,B,C]`  the three H spot numbers to use as registration spots.  Can omit when the cell has exactly three H spots.  

# Getting started

## Clone repo and then use uv to install Python dependencies (including cellpose 2 and jupyterlab)
    # Assuming R installed
    # Assuming uv installed (see https://docs.astral.sh/uv/getting-started/installation/)

    git clone https://github.com/OMRF/pezza-lab-stepwise-pairing-and-fast-transition-regulated-by-recomb-drive-meiotic-chrom-interactions1.git stepwise-pairing-fast-transition
    cd stepwise-pairing-fast-transition
    uv sync 

## Testing JupyterLab notebooks:

Please create the directory "sample_img" and download the file CS4_g1.czi to it.

Start with the notebook "SpotSelection.ipynb" and afterwards run the notebook "Nucleus_Proofreading_and_Classification.ipynb".

