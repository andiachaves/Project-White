if(!require(compiler))
{
  install.packages("compiler", repos = mir)
  require(compiler)
}


#' ## script for converting .Rmd files to .R scripts.

#' #### Kevin Keenan 2014
#' http://rstudio-pubs-static.s3.amazonaws.com/12734_0a38887f19a34d92b7311a2c9cb15022.html
#'
#' This function will read a standard R markdown source file and convert it to
#' an R script to allow the code to be run using the "source" function.
#'
#' The function is quite simplisting in that it reads a .Rmd file and adds
#' comments to non-r code sections, while leaving R code without comments
#' so that the interpreter can run the commands.
#'
#'
# library "stringi" needs to be loaded
rmd2rscript <- function(infile){
  # read the file
  flIn <- readLines(infile)
  # identify the start of code blocks
  cdStrt <- which(grepl(flIn, pattern = "```{r*", perl = TRUE))
  # identify the end of code blocks
  cdEnd <- sapply(cdStrt, function(x){
    preidx <- which(grepl(flIn[-(1:x)], pattern = "```", perl = TRUE))[1]
    return(preidx + x)
  })
  # define an expansion function
  # strip code block indacators
  flIn[c(cdStrt, cdEnd)] <- ""
  expFun <- function(strt, End){
    strt <- strt+1
    End <- End-1
    return(strt:End)
  }
  idx <- unlist(mapply(FUN = expFun, strt = cdStrt, End = cdEnd,
                       SIMPLIFY = FALSE))
  # add comments to all lines except code blocks
  comIdx <- 1:length(flIn)
  comIdx <- comIdx[-idx]
  for(i in comIdx){
    flIn[i] <- paste("#' ", flIn[i], sep = "")
  }
  # create an output file
  #nm <- strsplit(infile, split = "\\.")[[1]][1]
  nm <- stri_sub(infile,1,-5)
  flOut <- file(paste(nm, "_rmd2r.R", sep = ""), "w")
  for(i in 1:length(flIn)){
    cat(flIn[i], "\n", file = flOut, sep = "\t")
  }
  close(flOut)
}



### POSTGRESQL

dbExistsTable <- function (con, name, ...)
{
  as.logical (
    dim (
      dbGetQuery (con,
                  paste ("select schemaname,tablename from pg_tables where schemaname='",
                         rev(strsplit(name, ".", fixed=TRUE)[[1]])[2],
                         "' and tablename='",
                         rev(strsplit(name, ".", fixed=TRUE)[[1]])[1], "'", sep="")
      )) [1])
}

dbRemoveTable <- function (con, name, ..., cascade=FALSE)
{
  if (dbExistsTable (con, name)) {
    dbGetQuery (con, paste ("drop table ", name, ifelse (cascade, "
                                                         cascade", ""), ";", sep=""))
  }
  }








detachAllPackages <- function() {

  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

  package.list <- setdiff(package.list,basic.packages)

  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}







#******************************  Data exploration ***************************#


# USEFUL DEFINITIONS
unit <- list(
  mgl = expression(paste("(",mg~l^{-1},")")),
  mgkg = expression(paste("(",mg~kg^{-1},")")),
  mugl = expression(paste("(",mu,g~l^{-1},")")),
  mugkg = expression(paste("(",mu,g~kg^{-1},")")),
  mum = expression(paste("(",mu,m,")")),
  temp = expression(paste("T (",degree,"C)")),
  oxy = expression(paste(O[2]," (",mg~l^{-1},")")),
  phyto= expression(paste("Phytoplankton (",mg~l^{-1},")")),
  phos = expression(paste(PO[4]^{3-phantom(0)},"(",mu,g~l^{-1},")"))
)


# Standard error of a mean
se<-function(x) sqrt(var(x)/length(x))


# Error bars for boxplots
#Usage error.bar(means, standard errors and bar labels)
error.bars<-function(yv,z,nn,ylab,ylim,main,cex.names){
  xv<-barplot(yv,ylim=ylim,names=nn,ylab=ylab,main=main,cex.names=cex.names)
  g=(max(xv)-min(xv))/50
  for (i in 1:length(xv)) {
    lines(c(xv[i],xv[i]),c(yv[i]+z[i],yv[i]-z[i]))
    lines(c(xv[i]-g,xv[i]+g),c(yv[i]+z[i], yv[i]+z[i]))
    lines(c(xv[i]-g,xv[i]+g),c(yv[i]-z[i], yv[i]-z[i]))
  }}



# rounding up:
# k=2, 321-> 400, 3221->3300
# k=3  321 -> 1000, 3221 -> 4000
roundout <- function(x,k){(floor( abs(x/10^k) +0.499999999999 ) +1)*10^k}



# 3 significant digits
NumToString <- function(nval, sigdig=3, xpcut=4) {   # Filter out zero as a special case; calculate digits after the decimal
  # from number of significant digits
  if (nval == 0) digs = sigdig - 1 else digs <- sigdig - 1 - floor(log10(abs(nval)))
  # Can't be negative
  if (digs < 0) digs <- 0
  # Switch to scientific notation if digits after the decimal > xp
  if (digs > xpcut) { fchar <- "e"; digs <- sigdig - 1 } else fchar <- "f"
  sprintf(paste("%.", digs, fchar, sep=""), nval) }
#library(broman) # rounding





###########################################
#                                         #
#          Correlation plots              #
#                                         #
###########################################



cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X,method="spearman",use="pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R
}





## put correlations on the panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=1, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))

  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

panel.smooth2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = 1, ...)
}


panel.lines2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                       cex = 1, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}

}




panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}




# function for upper panel in plot using cor.test()
panel.cor.spear <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x, y,method="spearman")$estimate
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

  test <- cor.test(x,y,method="spearman")
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))

  text(0.5, 0.5, txt, cex = cex * abs(r))
  text(.8, .8, Signif, cex=cex, col=2)
}


# function for upper panel in plot using lm()
"panel.lm" <-
  function (x, y, col = par("col"), bg = NA, pch = par("pch"),
            cex = 1, col.lm = "red",  ...)
  {   ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylim = ylim, xlim= xlim)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    abline(lm(y[ok]~ x[ok]),
           col = col.lm, ...)
  }





# usage:
# pairs(data, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(..., nrpoints=0)},
# diag.panel=panel.hist, upper.panel=panel.cor)






#######################################################
#                                                     #
#   Max. combinations of environmental predictors     #
#                                                     #
#######################################################


# Scripts allows to select environmental predictor based on
# subjective judgement using ecological expertise (interactive)
# Script takes correlation into account! (data.frames: predictors.cor -> df -> xy)
# combination is saved in "y"

# creation date: 19.04.2013



PredictorSelection<-function(df){

# SELECTION starts here:

# get rid of factors!, sonst Wanrmeldung beim letzten loop:
# Warnmeldung:
# In `[<-.factor`(`*tmp*`, iseq, value =
# "rs13,rs265,bgl,fk_klasse,nfk_klasse,kak_klasse,kfa_klasse,natveg,kultpfla,fipu_a_ln,fipu_o_ln,fipu_s_ln") :
#  invalid factor level, NAs generated

xy<-data.frame(rownames(df),lapply(df,as.character),stringsAsFactors=FALSE)


# SELECT variables with 0 correlation
y<-as.vector(xy[xy[,2]==0,1])
xy<-xy[-which(xy[,2]==0),] # delete selected variables


gg<-NULL
while(nrow(xy)>0){
cat("\nLeft variables:",nrow(xy),"\n")
print(as.vector(xy[,1]))
#print(xy)

x<-readline("\nSelect an additional environmental variable! ")
while(!(x %in% xy[,1])){
x<-readline("\nWrong spelling - please try again! ")
}


# STEP 1: add variable to selection
y<-c(y,x)

# STEP 2: delete correlated variables from dataframe
x.cor1<-strsplit(as.vector(xy[which(xy[,1]==x),2]),"\\,") # select corr variables
x.cor2<-which(xy[,1] %in% x.cor1[1][[1]]) # which rows?
if(length(x.cor2)!=0){
xy<-xy[-x.cor2,] # delete selected variables (their rows)
}
# STEP 3: delete variable x
xy<-xy[-which(xy[,1] %in% x),]



# Step 4: delete correlated variables from the "correlated variables column" of other variables
if(nrow(xy)>0){
for (i in 1:nrow(xy)){
g<-strsplit(as.vector(xy[i,2]),"\\,")[[1]] # splitted string of correlated variables(column 2) belonging to variable i
xy[i,2]<-paste(g[!(g %in% x.cor1[1][[1]])],collapse=",") # get
if(length(g[!(g %in% x.cor1[1][[1]])])==0){ # if variable has 0 correlated variables -> add it to y
y<-c(y,as.vector(xy[i,1]))
gg<-c(gg,which(xy[,1] %in% xy[i,1])) # record rownumber of variable to be dropped
}

}
if (length(gg)>0){
xy<-xy[-gg,]
gg<-NULL}
}

}



pred.uncor<-y
print(pred.uncor)

return(pred.uncor)
}






#######################################
#                                     #
#              GLM                    #
#                                     #
#######################################



# store warning message if it occurs
keepWarnings <- function(expr) {
  localWarnings <- list()
  value <- withCallingHandlers(expr,
                               warning = function(w) {
                                 localWarnings[[length(localWarnings)+1]] <<- w
                                 invokeRestart("muffleWarning")
                               })
  list(value=value, warnings=localWarnings)
}






#### Function to calculate all possible combinations without duplicates in different columns  ###

# highly correlated variables not in the same model (spearman >0.5)


combs=NULL
allcombs<-function(x) {
  nx<-length(x)

  if(nx>=4){
    #*****
    for(i in 2:4) { # 4 = max number of variables in models

      indices<-combn(1:nx,i)


      #~~~ indices (combinations) check for correlations ~~~#

      for (p in 1:nx){

        #if (predictors.cor[p,1]!=0){
        pcc<-as.character(predictors.cor[p,1])
        pcc.spl<-unlist(strsplit(pcc,","))

        if (i==2){
          if(length(which(indices[1,]==p & indices[2,] %in% which(x %in% pcc.spl)))>0){
            meet.crit <- which(indices[1,]==p & indices[2,] %in% which(x %in% pcc.spl))
            indices<-indices[,-meet.crit]}}

        if (i==3){
          if(length(which(indices[1,]==p &
                          (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl))
          ))>0){
            meet.crit <- which(indices[1,]==p &
                                 (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl)))
            indices<-indices[,-meet.crit]}

          # Are there correlation within model predictors? ,i.e. btw 2nd and 3rd row
          usec<-unique(indices[2,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc2<-as.character(predictors.cor[usec[z],1])
            pcc2.spl<-unlist(strsplit(pcc2,","))
            if(length(which(indices[1,]==p & indices[2,]==usec[z] &
                            (indices[3,] %in% which(x %in% pcc2.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[2,]==usec[z] &
                                   (indices[3,] %in% which(x %in% pcc2.spl)))
              indices<-indices[,-meet.crit]
            }}
        }# closes i=3 loop

        if (i==4){
          if(length(which(indices[1,]==p &
                          (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl) | indices[4,] %in% which(x %in% pcc.spl))
          ))>0){
            meet.crit <- which(indices[1,]==p &
                                 (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl) | indices[4,] %in% which(x %in% pcc.spl)))
            indices<-indices[,-meet.crit]
          }

          # Are there correlation within model predictors? ,i.e. btw 2nd AND 3rd and 4thw
          usec<-unique(indices[2,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc2<-as.character(predictors.cor[usec[z],1])
            pcc2.spl<-unlist(strsplit(pcc2,","))
            if(length(which(indices[1,]==p & indices[2,]==usec[z] &
                            (indices[3,] %in% which(x %in% pcc2.spl)| indices[4,] %in% which(x %in% pcc2.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[2,]==usec[z] &
                                   (indices[3,] %in% which(x %in% pcc2.spl)| indices[4,] %in% which(x %in% pcc2.spl)))
              indices<-indices[,-meet.crit]
            }}



          # Are there correlation within model predictors? ,i.e. btw 3rd and 4th row
          usec<-unique(indices[3,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc3<-as.character(predictors.cor[usec[z],1])
            pcc3.spl<-unlist(strsplit(pcc3,","))
            if(length(which(indices[1,]==p & indices[3,]==usec[z] &
                            (indices[2,] %in% which(x %in% pcc3.spl)| indices[4,] %in% which(x %in% pcc3.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[3,]==usec[z] &
                                   (indices[2,] %in% which(x %in% pcc3.spl)| indices[4,] %in% which(x %in% pcc3.spl)))
              indices<-indices[,-meet.crit]
            }}

        }# closes i=4 loop


      }
      #}

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

      for(j in 1:dim(indices)[2]) {
        xlist<-list(x[indices[,j]])
        combs<-c(combs,xlist)
      }


    }
  } #close if clausen nx >=4


  ################################################################################################




  if(nx==2){
    #*****
    for(i in 2) { # 4 = max number of variables in models

      indices<-combn(1:nx,i)


      #~~~ indices (combinations) check for correlations ~~~#

      for (p in 1:nx){

        #if (predictors.cor[p,1]!=0){
        pcc<-as.character(predictors.cor[p,1])
        pcc.spl<-unlist(strsplit(pcc,","))

        if (i==2){
          if(length(which(indices[1,]==p & indices[2,] %in% which(x %in% pcc.spl)))>0){
            meet.crit <- which(indices[1,]==p & indices[2,] %in% which(x %in% pcc.spl))
            indices<-indices[,-meet.crit]}}

        if (i==3){
          if(length(which(indices[1,]==p &
                          (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl))
          ))>0){
            meet.crit <- which(indices[1,]==p &
                                 (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl)))
            indices<-indices[,-meet.crit]}

          # Are there correlation within model predictors? ,i.e. btw 2nd and 3rd row
          usec<-unique(indices[2,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc2<-as.character(predictors.cor[usec[z],1])
            pcc2.spl<-unlist(strsplit(pcc2,","))
            if(length(which(indices[1,]==p & indices[2,]==usec[z] &
                            (indices[3,] %in% which(x %in% pcc2.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[2,]==usec[z] &
                                   (indices[3,] %in% which(x %in% pcc2.spl)))
              indices<-indices[,-meet.crit]
            }}
        }# closes i=3 loop

        if (i==4){
          if(length(which(indices[1,]==p &
                          (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl) | indices[4,] %in% which(x %in% pcc.spl))
          ))>0){
            meet.crit <- which(indices[1,]==p &
                                 (indices[2,] %in% which(x %in% pcc.spl) | indices[3,] %in% which(x %in% pcc.spl) | indices[4,] %in% which(x %in% pcc.spl)))
            indices<-indices[,-meet.crit]
          }

          # Are there correlation within model predictors? ,i.e. btw 2nd AND 3rd and 4thw
          usec<-unique(indices[2,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc2<-as.character(predictors.cor[usec[z],1])
            pcc2.spl<-unlist(strsplit(pcc2,","))
            if(length(which(indices[1,]==p & indices[2,]==usec[z] &
                            (indices[3,] %in% which(x %in% pcc2.spl)| indices[4,] %in% which(x %in% pcc2.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[2,]==usec[z] &
                                   (indices[3,] %in% which(x %in% pcc2.spl)| indices[4,] %in% which(x %in% pcc2.spl)))
              indices<-indices[,-meet.crit]
            }}



          # Are there correlation within model predictors? ,i.e. btw 3rd and 4th row
          usec<-unique(indices[3,which(indices[1,]==p)]) # unique values 2nd row
          for (z in 1:length(usec)){
            pcc3<-as.character(predictors.cor[usec[z],1])
            pcc3.spl<-unlist(strsplit(pcc3,","))
            if(length(which(indices[1,]==p & indices[3,]==usec[z] &
                            (indices[2,] %in% which(x %in% pcc3.spl)| indices[4,] %in% which(x %in% pcc3.spl))
            ))>0){
              meet.crit <- which(indices[1,]==p & indices[3,]==usec[z] &
                                   (indices[2,] %in% which(x %in% pcc3.spl)| indices[4,] %in% which(x %in% pcc3.spl)))
              indices<-indices[,-meet.crit]
            }}

        }# closes i=4 loop


      }
      #}

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

      for(j in 1:dim(indices)[2]) {
        xlist<-list(x[indices[,j]])
        combs<-c(combs,xlist)
      }


    }
  } #close if clausen nx ==2







  #*****
  return(combs)
}







###################  Function to prepare input data for GLM formula ##########################


pred.function<-function(pred){

  ### Quadratic terms
  # add variables to quadratic terms if necessary
  # until now, combinations in this style: slope.sq + av.sq + FK_KLASSE
  # next loop to add normal variables were appropiate: slope + slope.sq + av + av.sq + FK_KLASSE

  x<-pred

  if (length(x[grepl(".sq", x)])>0){
    pred2<-rep(0,length(pred)+length(x[grepl(".sq", x)]))

    p=0
    for (s in 1:length(pred)){
      p=p+1
      if((pred[s] %in% x[grepl(".sq", x)])==FALSE){
        pred2[p]<-pred[s]}

      if((pred[s] %in% x[grepl(".sq", x)])==TRUE){
        pred2[p]<-gsub(".sq","",pred[s])
        p=p+1
        pred2[p]<-pred[s]}
    }

    pred<-pred2

  }
  return(pred)}


df.function<-function(pred){
  z<-length(pred) # index of variables
  xy.reg<-data.frame(matrix(0,length(hyal_pres),1+length(pred))) # define new dataset for model
  # response
  xy.reg[,1]<-hyal_pres
  colnames(xy.reg)[1]<-"response"
  # predictors
  for (j in 1:z){
    xy.reg[,j+1]<-glm.preds[pred[j]]
    colnames(xy.reg)[j+1]<-pred[j]
  }
  return(xy.reg)
}


df.function.138<-function(pred){
  z<-length(pred) # index of variables
  xy.reg<-data.frame(matrix(0,length(hyal_pres),1+length(pred))) # define new dataset for model
  # response
  xy.reg[,1]<-hyal_pres
  colnames(xy.reg)[1]<-"response"
  # predictors
  for (j in 1:z){
    xy.reg[,j+1]<-glm.preds.138[pred[j]]
    colnames(xy.reg)[j+1]<-pred[j]
  }
  return(xy.reg)
}


df.function.schriesheim<-function(pred){
  z<-length(pred) # index of variables
  xy.reg<-data.frame(matrix(0,length(y),1+length(pred))) # define new dataset for model
  # response
  xy.reg[,1]<-y
  colnames(xy.reg)[1]<-"response"
  # predictors
  for (j in 1:z){
    xy.reg[,j+1]<-glm.preds[pred[j]]
    colnames(xy.reg)[j+1]<-pred[j]
  }
  return(xy.reg)
}


# total run time
# Usage:
# before run enter: start.time<-Sys.time()
# after run enter: end.time<-Sys.time();run.time(end.time)


run.time<-function(end.time){
  run.sec<-trunc(difftime(end.time,start.time,unit="sec"))[[1]]
  run.min<-trunc(difftime(end.time,start.time,unit="mins"))[[1]]
  run.hours<-trunc(difftime(end.time,start.time,unit="hours"))[[1]]

  if(run.hours==0){
    m<-run.min*60
    run.sec2<-run.sec-m
    print(paste("Total run time: ",run.min," minutes and ",run.sec2," seconds",sep=""))}

  if(run.hours>0){
    h<-run.hours*60
    run.min2<-run.min-h
    m<-run.min*60
    run.sec2<-run.sec-m
    print(paste("Total run time: ",run.hours," hours and ", run.min2," minutes and ",run.sec2," seconds",sep=""))}
}








########################################################################
#                                                                      #
#   checking for duplicates, could be more efficient with binomial     #
#                                                                      #
########################################################################


# function expects all environmental variables as a DF, with the response (1/2) as the LAST column

combineBinomial <- function(data,nvariables=1){
  data <- data[do.call(order,data),]
  data$duplicated = duplicated(data)
  numPred <- dim(data)[1]
  data$illsum <- rep(NA, numPred)
  data$Nsum <- rep(NA, numPred)

  lastUnique <- 1
  matdata <- as.matrix(data[,(dim(data)[2]-3):dim(data)[2] ])

  for (row in 1:numPred){
    if (matdata[row,2]==T){
      matdata[lastUnique,3] = matdata[lastUnique,3] + matdata[row,1]
      matdata[lastUnique,4] = matdata[lastUnique,4] + 1
    }
    else{
      lastUnique <- row
      matdata[row,3] = matdata[row,1] # put healthy(0)/ill(1) in
      matdata[row,4] = 1
    }
  }

  data[,(dim(data)[2]-3):dim(data)[2] ] = matdata


# PartII
#> head(tmpfdf,3)
#id ill duplicated illsum Nsum
#1   1   0          0      0   36
#37  2   0          0      0   48
#51  2   1          0      2    2

tmpdf<-data[!data$duplicated,]

tmpdf$duplicated<-duplicated(tmpdf[,1:nvariables]) # checks for more >1 entry for id/site
TrueDuplicates<-which(tmpdf$duplicated==T)
NoDuplicates<-which(tmpdf$duplicated==F)

for(duplicates in TrueDuplicates){

  # 1. step get total ills/ successes, e.g. 2 aus Zeile 3 in Zeile 2
  tmpdf[duplicates-1,(ncol(tmpdf))-1]<-tmpdf[duplicates,(ncol(tmpdf))-1]
  # 2. step sum ill and healthy to get total
  tmpdf[duplicates-1,(ncol(tmpdf))]<-tmpdf[duplicates,(ncol(tmpdf))]+tmpdf[duplicates-1,(ncol(tmpdf))]
}

# 3. step remove duplicates
finaldf<-tmpdf[NoDuplicates,]

  return(finaldf)
}

combineBinomial <- cmpfun(combineBinomial)



# TEST DATA SET
# site <- c(1,1,1,2,2,2,3,3,3)
# elevation <- c(0,0,0,10,20,30,10,10,20)
# slope <- c(0,0,0,0,10,10,10,10,10)
# inf <- c(1,0,0,1,1,1,0,0,0)
# testDf <- data.frame(site,elevation,slope,inf)
# data <- testDf[,-1]


# > testDf[,-1]
#   elevation slope inf
# 1         0     0   1
# 2         0     0   0
# 3         0     0   0
# 4        10     0   1
# 5        20    10   1
# 6        30    10   1
# 7        10    10   0
# 8        10    10   0
# 9        20    10   0

# function should combine 1,2,3 and 7,8 and 5,9

#after part I tmpdf:
# > tmpdf
#   elevation slope inf duplicated illsum Nsum
# 2         0     0   0          0      0    2
# 1         0     0   1          0      1    1
# 4        10     0   1          0      1    1
# 7        10    10   0          0      0    2
# 9        20    10   0          0      0    1
# 5        20    10   1          0      1    1
# 6        30    10   1          0      1    1




# function to disentangle binomial response(success/total observations) into bernoulli (0 or 1)
# first success, second total observations, last: predictions


# illsum, Nsum = vector, predicted = matrix

untangleBinomial<-function(illsum,Nsum,predicted){

  predictionsExpanded <-apply(predicted,2,function(x) rep(x,Nsum))


  ones <-mapply(rep, 1, illsum)
  zeros<-mapply(rep, 0, Nsum-illsum)


  recurse <-function(l1,l2){
    if(is(l1[[1]],"list"))
      mapply(recurse,l1,l2,SIMPLIFY=F)
    else
		  {mapply(c,l1,l2,SIMPLIFY=F)}
	 }

  observed <-unlist(recurse(ones,zeros))

return(data.frame(observed, predictionsExpanded ))
}


untangleBinomial <- cmpfun(untangleBinomial)


# TEST FOR untangleBinomial
# e.g.
# test data
# illsum<-c(2,0,1,4)
# Nsum<-c(4,1,8,4)
# predicted<-cbind(c(0.3,0.6,0.8,0.2),c(0.3,0.6,0.8,0.2))
#untangleBinomial(illsum, Nsum, predicted)




#####################
#                   #
#  AUC calculation  #
#                   #
#####################

#lumpchains combines the chains in mcmc.list
# copied from emdbook package by Ben Bolker
lumpchains <- function(x){
  x2 <- do.call("rbind", x)
  mcpars <- sapply(x, attr, "mcpar")
  class(x2) <- "mcmc"
  if (var(mcpars[1, ]) > 0 || var(mcpars[3, ]) > 0)
    stop("can't combine chains with unequal start/thin")
  attr(x2, "mcpar") <- c(mcpars[1, 1], sum(mcpars[2, ]),
                         mcpars[3,1])
  x2
}



# predicts from a stan model
# matrix ultiplications: d = vector, dd = data.frame
# 1) apply(as.matrix(dd)%*%diag(d),1,sum)
# 2) t(t(dd)*d)

# if random effects are used, optional arguments need to be specified
getPredictions <- function(pars, environmentVariables,numberRandomeffects=0,randomEffects=NULL,nObs=NULL,hyperPriors=NULL){

if(numberRandomeffects==2){
mu<-rep(0,nObs)
for( i in 1:nObs ) {
  mu[i] <- pars[paste(hyperPriors[1],paste("[",randomEffects[,1][i],"]",sep=""),sep="")] + pars[paste(hyperPriors[2],paste("[",randomEffects[,2][i],"]",sep=""),sep="")] # u and v need quotes!
    }
linear<- pars["b0"] +
apply(as.matrix(environmentVariables)%*%diag(pars[grep("par",names(pars))]),1,sum) + mu
}

if(numberRandomeffects==1){
mu<-rep(0,nObs)
for( i in 1:nObs ) {
  mu[i] <- pars[paste(hyperPriors[,1],paste("[",randomEffects[,1][i],"]",sep=""),sep="")] # u and v need quotes!
    }

linear<- pars["b0"] +
apply(as.matrix(environmentVariables)%*%diag(pars[grep("par",names(pars))]),1,sum) + mu
}

if(numberRandomeffects==0){
linear<- pars["b0"] +
apply(as.matrix(environmentVariables)%*%diag(pars[grep("par",names(pars))]),1,sum)
}
return(as.data.frame(inv.logit(linear)))
}






getAUC <- function(illsum,Nsum,bionomialPredictions, makeplot = F){

  expandedpredictions <- untangleBinomial(illsum,Nsum,bionomialPredictions)

  combined<-prediction(expandedpredictions[,2:(ncol(bionomialPredictions)+1)],matrix(rep(expandedpredictions[,1],ncol(bionomialPredictions)),ncol= ncol(bionomialPredictions)))

  if (makeplot==T){
    Perform<-performance(combined,"tpr", "fpr")
    plot(Perform, avg='vertical', spread.estimate='boxplot')
  }
  Perform<-performance(combined,"auc")
  return(unlist(Perform@y.values))
}




################
#              #
#  Stan graph  #
#              #
################

# libraries: rstan, MCMCpack, plyr
stanGraphInput<-function(stan.fit,n.ch=4,q=0.1,parGrep="par|b0"){
# extract different elements used for plotting
# returns 2 objects: datq (used for conditional distribution of parameters) and var.median
#$datq
#   X1 X2            X3           X4
#1   1  1    0.75554996   5.61069092
#2   2  2    1.47901014   4.64270081
#3   3  3   -2.19399495   3.00084311
#4   4  4   -0.09040486   3.53912401
#5   5  5    0.46681738   6.18981024
#6   6  6   -0.84793811   2.00565756
# Number - lowerX,upperX, lowerY,upperY



s<-rstan::extract(stan.fit,permuted=F,inc_warmup=F) # extract samples, library rstan
s.coda <- do.call(mcmc.list, alply(s, 2, mcmc)) # create coda object, libraries: MCMCpack, plyr

dnam<-dimnames(s.coda[[1]])$parameters # get just the names of the parameters
grepdnam<-which(grepl(parGrep,dnam))

# create empty data.frame for medians
var.median<-data.frame(matrix(0,nrow=length(grepdnam)*n.ch,ncol=2)) # nrow*4 -> 4 medians
var.median[,1]<-rep(1:length(grepdnam),each=n.ch)
# create empty data.frame for lower and upper interval
datq<-data.frame(matrix(0,nrow=length(grepdnam),ncol=4)) # ncol=4 -> lower & upper interval + 2 values for x axis
datq[,1:2]<-1:length(grepdnam) # values for x-axis of lines


# fill data.frames
h=1
for(g in 1:length(grepdnam)){ # loop over all specified variables
var<-dimnames(s.coda[[1]])$parameters[grepdnam[g]]
var.median[h:(h+3),2]<-sapply(1:n.ch,function(m) median(s.coda[,var][[m]])) # calculates medians from the 4 chains
h<-h+n.ch
if (n.ch==1){dat<-s.coda[,var][[1]]}
if (n.ch==2){dat<-c(s.coda[,var][[1]],s.coda[,var][[2]])}
if (n.ch==3){dat<-c(s.coda[,var][[1]],s.coda[,var][[2]],s.coda[,var][[3]])}
if (n.ch==4){dat<-c(s.coda[,var][[1]],s.coda[,var][[2]],s.coda[,var][[3]],s.coda[,var][[4]])}
datq[g,3:4]<-quantile(dat,probs=c(q,1-q)) # get upper and lower interval
}
return(list("parameters"=dnam,"var.median" = var.median,"datq"= datq))
}






# Function for data preparation

combineSitePredictions<-function(df1,df2){
  # df1: last column: site
  # df2: first column site, second column predictions

  #alternative: zz <- merge(df1, df2, all = TRUE)


  # dataset examples
  #id<-c("1","2","3","4","5","6")
  #site<-c("1","1","2","3-4","1-4","3-4")
  #df1<-data.frame(id,site)

  #site<-as.factor(c(1,2,3,4))
  #hyal<-c(0.2,1,0.4,0.8)
  #df2<-data.frame(site,hyal)

  # merging
  zz<-join(df1, df2, type="left")

  na.sites<-which(is.na(zz[,ncol(zz)])) # which rows contain overlapping sites and sites not found in df2
  u.na.sites<-unique(zz[na.sites,ncol(zz)-1]) # unique na sites
  u.overlap.sites<-as.character(u.na.sites[grep("-",u.na.sites)]) # only overlapping sites, for other there is not match in df2 anyway


  list.overlap.sites<-strsplit(x=u.overlap.sites, split="-") # convert "-" in list
   # checks if e.g. site "417-418" both numbers are in df2
  for (i in 1:length(u.overlap.sites)){
  if(suppressWarnings(!is.na(mean(df2[df2[,1] %in% list.overlap.sites[[i]],2])))==TRUE){
  zz[zz[,ncol(zz)-1]==u.overlap.sites[i],ncol(zz)]<-mean(df2[which(df2[,1] %in% list.overlap.sites[[i]]),2])
  }}

return(zz)}




# store warning message if it occurs
keepWarnings <- function(expr) {
  localWarnings <- list()
  value <- withCallingHandlers(expr,
  warning = function(w) {
  localWarnings[[length(localWarnings)+1]] <<- w
  invokeRestart("muffleWarning")
  })
  list(value=value, warnings=localWarnings)
}


# Logit's
logit<- function(x) log(x/(1-x))
inv.logit<-function(x) 1/(1+exp(-x))

#inv.logit<-function(x) plogis(x)  # exp(x)/(1+exp(x))
#logit<-function(x) qlogis(x)






################
#              #
#  Diagnosis   #
#              #
################

# function calculates Bayesian p-value

# Was das macht ist dass ich das quantil einfach auf die Mitte setze, d.h wenn die Simulation 50% Nullen vorhersagt, und 0
# beobachtet ist, dann sollte der Wert 0.25 sein.

#   <         <=
# sim obs  sim  obs
#  0   0    0    0
#  0   0    0    0
#  0   0    0    0
#  1   0    1    0
#  1   0    1    0
#  1   0    1    0

#=  0        0.5     -> 0.25

# es sind skalierte residuen, und sie sind auch nicht in Prozent, sondern quantile ... ich würde schreiben über die Plots
# schreiben "Quantile residuals", und über die Skala "residuals as quantiles from simulation"


# # example of current configuration
# simulatedValues<- c(2,4,6,8)
# observedValue <- 2
# low <- sum(simulatedValues < observedValue) / length(simulatedValues)
# up <- sum(simulatedValues <= observedValue) / length(simulatedValues)
# > low
# [1] 0
# > up
# [1] 0.25
# > mean(c(low, up))
# [1] 0.125
#
# observedValue <- 8
# low <- sum(simulatedValues < observedValue) / length(simulatedValues)
# up <- sum(simulatedValues <= observedValue) / length(simulatedValues)
# > low
# [1] 0.75
# > up
# [1] 1
# > mean(c(low, up))
# [1] 0.875


# > observedValue <- 10
# > low <- sum(simulatedValues < observedValue) / length(simulatedValues)
# > up <- sum(simulatedValues <= observedValue) / length(simulatedValues)
# > low
# [1] 1
# > up
# [1] 1
# > mean(c(low, up))
# [1] 1
#
# > observedValue <- 5
# > low <- sum(simulatedValues < observedValue) / length(simulatedValues)
# > up <- sum(simulatedValues <= observedValue) / length(simulatedValues)
# > low
# [1] 0.5
# > up
# [1] 0.5
# > mean(c(low, up))
# [1] 0.5



doUnivariateTest <- function(Site='', dfSim,dfObs,plot = T){
  simulatedValues <-dfSim  # e.g. length= 20.000; length of combined chains
  observedValue <- dfObs   # length = 1

  low <- sum(simulatedValues < observedValue) / length(simulatedValues)
  up <- sum(simulatedValues <= observedValue) / length(simulatedValues)

  pv = mean(c(low, up))

  #pv = lowerP
  #else if (upperP > 0.5) pv = upperP
  #else if (lowerP < 0.5 & upperP > 0.5) pv = 0.5

  if (plot == T){
    hist(simulatedValues, main = paste("Site = ",Site, "\nmean p-value = ",pv))
    abline(v=observedValue, col = "red", lwd = 3)
  }
  return(pv)
}

