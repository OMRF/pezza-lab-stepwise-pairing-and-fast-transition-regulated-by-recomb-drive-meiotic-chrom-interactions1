# Stepwise pairing and fast transition regulated by recombination drive meiotic chromosome interactions

# Scripts

## kymograph.R file.xls

`file.xls` is an XLS file containing the sheet "Distance To Nearest Neighbour"
with same-named column ("Distance To Nearest Neighbour") and a column named
"Time". Rows with duplicate Time values are ignored (we expect them to
correspond to completely duplicate rows).

## median_detection.R changeofphases.csv

file `changeofphases.csv` contains a list of numbers, one number per line.

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
