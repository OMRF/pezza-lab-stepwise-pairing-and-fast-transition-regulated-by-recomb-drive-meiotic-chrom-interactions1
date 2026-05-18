#!/bin/env -S Rscript --vanilla

######
##
## Change point analysis option for Dr. Roberto Pezza
## We will make use of 2 approaches: changepoint and "chow" method
##
## AUTHOR(s): Nathan Pezant (with a few modifications by Christopher Bottoms)
######
library(changepoint) # needed for change point method
library(strucchange) # needed for "chow" method.
library(readxl)      # for read_excel

pdf("change_of_phases.pdf")

# 3 different ways to provide a file name
# Uncomment your preferred way

# Harcoded variable.
#filename <- "test/input/demo.csv"

# Input CSV file name from the command line
#args=commandArgs(trailingOnly = TRUE)
#filename=args[1]

# OR Use interactive window to choose a file
filename=file.choose()

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

## plot the data to get an idea of what to expect
print(plot(Distances$`Distance To Nearest Neighbour`))

######
##
## section 1: changepoint
##
######

### changepoint has 3 options: cpt.mean, cpt.var, cpt.meanvar
### We will use all 3 and see which we like best

res.mean=cpt.mean(Distances$`Distance To Nearest Neighbour`) # based on region mean
res.var=cpt.var(Distances$`Distance To Nearest Neighbour`) # based on region variance
res.meanvar=cpt.meanvar(Distances$`Distance To Nearest Neighbour`) # based on region mean and variance

### to view the actual cut point from each analysis we can do the following
res.mean@cpts
res.var@cpts
res.meanvar@cpts

### note: the last point is always listed as a "change point".

n.mean=length(res.mean@cpts) ### retrieve the number of sections determined by the cutpoints.
n.var=length(res.var@cpts)
n.meanvar=length(res.meanvar@cpts)

### now we need to calculate the average (mean) for each region for each method

region_avg.mean=mean(Distances$`Distance To Nearest Neighbour`[1:res.mean@cpts[1]]) # calculate the mean for the first section
for(i in 2:n.mean){
  temp=mean(Distances$`Distance To Nearest Neighbour`[res.mean@cpts[i-1]:res.mean@cpts[i]]) # calculate the mean for each section
  region_avg.mean=c(region_avg.mean,temp)
  rm(temp)
}


region_avg.var=mean(Distances$`Distance To Nearest Neighbour`[1:res.var@cpts[1]]) # calculate the mean for the first section
for(i in 2:n.var){
  temp=mean(Distances$`Distance To Nearest Neighbour`[res.var@cpts[i-1]:res.var@cpts[i]]) # calculate the mean for each section
  region_avg.var=c(region_avg.var,temp)
  rm(temp)
}


region_avg.meanvar=mean(Distances$`Distance To Nearest Neighbour`[1:res.meanvar@cpts[1]]) # calculate the mean for the first section
for(i in 2:n.meanvar){
  temp=mean(Distances$`Distance To Nearest Neighbour`[res.meanvar@cpts[i-1]:res.meanvar@cpts[i]]) # calculate the mean for each section
  region_avg.meanvar=c(region_avg.meanvar,temp)
  rm(temp)
}

### create a table to display the "results"
### I need this to include 1) each index, and 2) its corresponding sectional mean

tab.mean=as.data.frame(matrix(nrow=length(Distances$`Distance To Nearest Neighbour`),ncol=2)) # create a table with 85 rows and 2 columns
tab.mean$V1=seq(1,length(Distances$`Distance To Nearest Neighbour`),1) # Fill in the indices

tab.var=as.data.frame(matrix(nrow=length(Distances$`Distance To Nearest Neighbour`),ncol=2)) # create a table with 85 rows and 2 columns
tab.var$V1=seq(1,length(Distances$`Distance To Nearest Neighbour`),1) # Fill in the indices

tab.meanvar=as.data.frame(matrix(nrow=length(Distances$`Distance To Nearest Neighbour`),ncol=2)) # create a table with 85 rows and 2 columns
tab.meanvar$V1=seq(1,length(Distances$`Distance To Nearest Neighbour`),1) # Fill in the indices


### put in the first avg to initialize
tab.mean$V2[1:res.mean@cpts[1]]=region_avg.mean[1]
### Loop through the other sections
for(i in 2:n.mean){
  tab.mean$V2[res.mean@cpts[i-1]:res.mean@cpts[i]]=region_avg.mean[i]
}

### put in the first avg to initialize
tab.var$V2[1:res.var@cpts[1]]=region_avg.var[1]
### Loop through the other sections
for(i in 2:n.var){
  tab.var$V2[res.var@cpts[i-1]:res.var@cpts[i]]=region_avg.var[i]
}

tab.meanvar$V2[1:res.meanvar@cpts[1]]=region_avg.meanvar[1]
### Loop through the other sections
for(i in 2:n.meanvar){
  tab.meanvar$V2[res.meanvar@cpts[i-1]:res.meanvar@cpts[i]]=region_avg.meanvar[i]
}

### Plot the results and save to file
print(plot(Distances$`Distance To Nearest Neighbour`))
points(tab.mean$V1,tab.mean$V2,col="red",pch=20)
points(tab.var$V1,tab.var$V2,col="blue",pch=20)
points(tab.meanvar$V1,tab.meanvar$V2,col="green",pch=20)
#dev.off()


########
##
## section 2: "Chow" method
## This method doesn't identify change points, but rather tests if any given
## point serves as a change point.
## We will test every point and report p-values, then adjust for multiple testing.
##
########

### create a data frame to track each index's p-value for being a break-point
break_point=as.data.frame(matrix(nrow=1,ncol=2))
names(break_point)=c("index","p-value") # we will track the break-point p-value for each point and record them here

## initialize the break-point resutls table with index 1 and p-value 1 (The first nor last points can be break-points)
break_point[1,]=c(1,1) # artificially setting the p-value for index 1 to p-value=1
for(i in 2:(nrow(Distances)-2)){ # loop through the remaining indices and record p-values
  temp=sctest(Distances$`Distance To Nearest Neighbour`~1,type="Chow",point=i) #the actual statistical test
  break_point[i,]=c(i,temp$p.value) # append newest results to the break_point table
  rm(temp)
}

break_point$padj=p.adjust(break_point$`p-value`)
to_add=min(break_point$padj[which(break_point$padj!=0)])/100
break_point$padj=break_point$padj+to_add ### we will be taking the log of these values and we can't log 0
break_point$neglog10p=-1*log10(break_point$padj)
print(plot(break_point$index,break_point$neglog10p,main=paste0("amount added to padj \n before log10 = ",formatC(to_add))))

dev.off()
