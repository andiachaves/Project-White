#' #Model
#'

# Date format - forms part of names of created files or graphs
today <- format(Sys.time(), "%Y.%m.%d")
set.seed(2017)

# load libraries
library(formatR)
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE) # automatically save the compilied Stan program
options(mc.cores = parallel::detectCores()) # execute multiple Markov chains in paral

library(MCMCpack) # rwish function
library(modeest)
library(IDPmisc)
#'
#'
#' ## Load data
#'
rm(list=ls(all=TRUE))
# Load settings
source("r-code/00_settings.r")

load(file="data/Whiteplague_workingdata.Rdata")
head(Whiteplague_workingdata)
dim(Whiteplague_workingdata)
#'
#'
#' ## Data preparation
#'
df <- Whiteplague_workingdata[,c("Ntot","ill","Plague","Chl_sca","Oxy_sca","Sal_sca","Turb_sca","Depth_tran_sca","Temp_sca","DHW_sca")]

df$healthy <- df$Ntot-df$ill
df$Month   <- as.factor(Whiteplague_workingdata$Month_12)
df$Site    <- as.factor(as.numeric(Whiteplague_workingdata$Station_name))
df$Obsid   <- as.factor(1:nrow(df))
df$SiteMonth <- as.factor(as.numeric(as.factor(paste(df$Month,df$Site,sep="-"))))

df$Site_index <- as.integer(df$Site)
df$Month_index=as.integer(df$Month)
df$Obsid_index=as.integer(df$Obsid)
df$SiteMonth_index=as.integer(as.numeric(as.factor(paste(df$Month,df$Site,sep="-"))))

Environment <- df[,c("Chl_sca","Oxy_sca","Sal_sca","Turb_sca","Depth_tran_sca","Temp_sca","DHW_sca")]
#'
#'
#'
#' ## Model definition
#'
dataList = list(
  TotalObservedcolonies = df$Ntot, #Total observed colonies Shallow (2 years) Deep (4 years)
  Diseasedcolonies = df$ill, # Total numer of disease colonies Shallow (2 years) Deep (4 years)
  nWP = nrow(df), #White Plague prevalence total number of rows
  observedWpPrevalence = df$Plague, # Total Prevalence of White Plague Shallow (2 years) Deep (4 years)

  Environment = Environment,
  nEnvi = ncol(Environment),

  SiteMonth = df$SiteMonth_index,
  nSiteMonth =length(unique(df$SiteMonth_index))
)
#'
#'
#'
#' ## Running the model
#'
#'
time <- proc.time()
stan.white <- stan("r-code/Whiteplague_rstanModel.stan",
                   data = dataList,
                   iter = 30000,
                   chains = 4,
                   control = list(adapt_delta = 0.99,max_treedepth = 15),
                   refresh=10)
(proc.time() - time)/60/60
