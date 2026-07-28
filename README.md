# Stepwise pairing and fast transition regulated by recombination drive meiotic chromosome interactions

# Scripts

## kymograph.R test/input/Example.xlsx

Input file (e.g., `Example.xlsx`) must contain a sheet named "Distance To Nearest Neighbour" with column headers in the second row.
One column header must be "Distance To Nearest Neighbour" and the other "Time".
Rows with duplicate Time values are ignored (we expect them to correspond to completely duplicate rows).

## median_detection.R test/input/Example.xlsx

Input file (e.g., `Example.xlsx`) must contain a sheet named "Distance To Nearest Neighbour" with column headers in the second row.
One column header must be "Distance To Nearest Neighbour" and the other "Time".
Rows with duplicate Time values are ignored (we expect them to correspond to completely duplicate rows).

## msd.R distance_file.csv

CSV file `distance_file.csv` contains the columns x, y, and z, separated by
commas. Other columns are ignored..  

Calculate mean square distance from each data point to the first data point.  

See https://en.wikipedia.org/wiki/Mean_squared_displacement#Derivation_for_n_dimensions.

# Getting started

## Clone repo and then use uv to install Python dependencies (including cellpose 2 and jupyterlab)
    # Assuming R installed
    # Assuming uv installed (see https://docs.astral.sh/uv/getting-started/installation/)

    git clone https://github.com/OMRF/pezza-lab-stepwise-pairing-and-fast-transition-regulated-by-recomb-drive-meiotic-chrom-interactions1.git stepwise-pairing-fast-transition
    cd stepwise-pairing-fast-transition
    uv sync 

## Testing JupyterLab notebooks:

Please create the directory "sample_img" and download the file CS4_g1.czi to it.
