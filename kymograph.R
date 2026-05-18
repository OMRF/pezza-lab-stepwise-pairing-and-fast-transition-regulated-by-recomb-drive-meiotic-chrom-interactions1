#!/bin/env -S Rscript --vanilla

## Turn on detailed feedback for errors
#options(error = function() traceback(3))

library(ggplot2) # plotting library
library(readxl)  # for read_excel

# 3 different ways to provide a file name
# Uncomment your preferred way

# Harcoded variable.
#filename <- "test/input/demo.xls"

# Input Excel file name from the command line
#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]

# Use interactive window to choose a file
filename <- file.choose()


# Read in Excel file, 
#  selecting the sheet with the distance data (i.e., "Distance To Nearest Neighbour"),
#  and skipping the first row (i.e., get column headers from the second row instead).
Distances_raw <- read_excel(filename,
                            sheet = "Distance To Nearest Neighbour",
                            skip = 1) # Skip first row

# Make sure the distances are treated as numbers
Distances_raw$`Distance To Nearest Neighbour` <-
  as.numeric(Distances_raw$`Distance To Nearest Neighbour`)

# Create TRUE/FALSE vector indicating rows with the first instance of each time point
#   (so far, duplicate time points correspond to duplicate rows, so just keeping one works)
keepers <-  ! duplicated(Distances_raw$Time)

# Keep just the desired rows 
#  restrictions for row selection is to the left of the comma 
#                                |
#                                | empty to the right of comma means no column restrictions
#                                | |  (i.e., keep all the columns)
Distances <- Distances_raw[keepers,] 

# Since we are plotting the distance between points, centered on 0
# Plot the first point half the distance to the left (- distance /2 )
# and the other point half the distance to the right (+ distance / 2)

Distances$leftDistance  <- - Distances$`Distance To Nearest Neighbour`/2
Distances$rightDistance <-   Distances$`Distance To Nearest Neighbour`/2

# Make sure time is treated as numeric
Distances$Time <- as.numeric(Distances$Time)

# Build up a ggplot2 object, one step at a time, notice the "+" symbols between
# each step

# Store this plot in the variable called "kymograph". We'll print it later.
kymograph <-
  
  # Create new ggplot2 object 
  ggplot(
    data = Distances, # Use Distances as the data for the plot
    aes(y = Time,     # Select Time for the Y axis
        col = Time)   # Use Time to assign color
  ) +
    
  # Add left-side points
  geom_point( 
    aes(x = leftDistance), # Plot x at leftDistance
    size = 1               # Make points size of 1
  ) +
    
  # Add right-side points
  geom_point(
    aes(x = rightDistance), # Plot x at rightDistance
    size = 1
  ) +
  
  # Start Y at the top
  scale_y_reverse() +
    
  # Limit the display to show a distance of at most 10
  xlim(-5,5) +
    
  # Remove the legend
  theme(legend.position = "none") +
  
  # And finally, label the X axis as micrometers (Unicode 00B5 is the micro symbol)
  xlab("distance (\U00B5m)")  # No more plus symbols (last ggplot layer)

# Now that we've created our plot object, let's print it to a file. At this 
# point, if you are using RStudio interactively, you can also type "kymograph" on
# the command line to see it in an RStudio window.

# Create output filename
out_filename <- sub(".xls(x)?$", ".kymograph.pdf", filename)
out_filename <- basename(out_filename)

pdf(out_filename) # Open new PDF file
print(kymograph)  # Print plot to it
dev.off()         # finalize PDF (close the file)
