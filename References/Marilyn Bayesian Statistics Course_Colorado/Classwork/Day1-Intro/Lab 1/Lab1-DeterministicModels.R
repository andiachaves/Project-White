# Lab 1: Deterministic models
# Exercise 1 - tadpoles data - try to find the correct function to model your data

setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Lectures\\")

data<-read.csv("tadpoles.csv")

x<-data$TBL

##Pick a function to apply to the data (use "Deterministic-functions-cheat-sheet4.pdf")
##- chose a Ricker function:

ricker = function(x, a = 1, b = 1)
{
a * x * exp(-b*x)
}

plot(data$TBL, data$Kill, col="blue", pch=16,xlim=c(0,40), ylim=c(0,10))

curve(ricker(x,a=2.5,b=0.2), add = TRUE)

##Ricker does an ok job of fitting the data
##But Logistic would do better

###Add a c parameter, so that we can shift the function


logistic <- function(x,a=-0.59,b=12,c=1.67)
{
exp(a*(b-x))/(1+exp(c*a*(b-x)))
}

curve(logistic(x), add=TRUE)

##Try a Michaelis-Menten function

michmenton <- function(x,a,b)
{
a*x/b+x
}

curve(michmenton(x,a=0.2,b=1), add=TRUE)


