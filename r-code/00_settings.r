###################
#### libaries #####
###################

library("ROCR") # prediction()
library(compiler) # fct: cmpfun(...)

source("r-code/00_helper_functions.r")
mir <- "http://cran.us.r-project.org" # define R repository

########################
#### date and time #####
########################

# Date format - forms part of names of created files or graphs
today<-format(Sys.time(), "%Y.%m.%d")

# The IDate class is a simple wrapper around the Date 
# class that tries to keep an integer storage format. The ITime class, the 
# time of day, is stored as the number of seconds in a day.
#source("C:/Programme/R/R-2.11.1/library/IDateTime.R")


# no scientific notation
options(scipen=999)
