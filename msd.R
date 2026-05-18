#!/bin/env -S Rscript --vanilla

# Mean Square Displacemnt

args <- commandArgs(trailingOnly = TRUE)

calc_msd <- function(df, x="x", y="y", z="z") {

    first_x <- df[1, x]
    first_y <- df[1, y]
    first_z <- df[1, z]

    # x,y,z positions excluding reference (first) position
    x_positions <- df[-1, x]
    y_positions <- df[-1, y]
    z_positions <- df[-1, z]

    # Calculate displacements from reference
    delta_x <- x_positions - first_x
    delta_y <- y_positions - first_y
    delta_z <- z_positions - first_z

    # Calculate means of squares for each dimension
    means_of_squares <- mean((delta_x + delta_y + delta_z)**2)

    # Combine all three dimensions
    msd <- sum(means_of_squares)

    return(msd)
}

# Run function on input CSV file
if (length(args) == 1) {
    input_df <- read.csv(args[1])
    message('Expecting column headers of "x", "y", and "z"')

    result <- calc_msd(input_df)
    cat(result)
} else if (length(args) == 4) {
    input_df <- read.csv(args[1])
    x <- args[2]
    y <- args[3]
    z <- args[4]
    message(paste0("Found x, y, z column headers: ", x, ", ", y,  ", ", z))

    result <- calc_msd(input_df, x, y, z)
    cat(result)
}
