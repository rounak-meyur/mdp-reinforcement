library(Matrix)
library(spatstat.utils)

#-----------------------------------------------------------------------------------------------------------------------------------------
# Set up directories
setwd("~/GitHub/mdp-reinforcement/src/R")
parent.dir <- dirname(getwd())
data.dir <- file.path(parent.dir,"case-data")

# Get the csv data
branchdat <- read.csv(file.path(data.dir,"branchdat.csv"))
gendat <- read.csv(file.path(data.dir,"gendat.csv"))
busdat <- read.csv(file.path(data.dir,"busdat.csv"))



### Function definitions
get.data <- function(bus,gen,branch)
{
  # Construct the gmin and gmax vectors
  busnum <- bus$bus_id
  genbus <- gen$bus_id
  nb <- nrow(bus)
  ng <- nrow(gen)
  gmax <- rep(0,nb)
  gmin <- rep(0,nb)
  for(i in 1:ng)
  {
    busind <- match(genbus[i],busnum)
    gmin[busind] = gmin[busind] + gen$pmin[i]/100.0
    gmax[busind] = gmax[busind] + gen$pmax[i]/100.0
  }
  
  # Construct the dmin and dmax vectors
  dmin <- rep(0,nb)
  dmax <- bus$pd/100.0
  
  # Construct the line rating vector
  flim <- branch$rateA/100.0
  
  # Return as an object
  return(list("gmin"=gmin,"gmax"=gmax,"dmax"=dmax,
              "dmin"=dmin,"flim"=flim))
}

make.PTDF = function(bus,branch)
{
  nb <- nrow(bus)
  nl <- nrow(branch)
  busnum <- bus$bus_id
  
  # Shivani: This requires further update (Insert code here)
  # For each component find a slack bus and append it
  # to this vector
  # Then iterate through each component and make sure that
  # each component has a slack bus, else identify a bus as
  # a slack bus
  # Change the entry of bus type for these buses to 3
  
  
  # The columns and rows to be chosen
  noslack <- (1:nb)[bus$type!=3]
  
  # Compute the PTDF matrix
  colind <- c(match(branch$fbus,busnum),
              match(branch$tbus,busnum))
  rowind <- c(seq(1,nl),seq(1,nl))
  entries <- c(rep(1,nl),rep(-1,nl))
  A <- sparseMatrix(i=rowind,j=colind,
                    x=entries,dims=c(nl,nb))
  
  # Note that I used the status column of the branchdat file
  # When you remove a line, just change the entry of status
  # column corresponding to the branch to 0
  b <- (branch$status)/(branch$x)
  
  # These are inv(X)*A and t(A)*inv(X)*A respectively
  Bf <- sparseMatrix(i=rowind,j=colind,
                     x=c(b,-b),dims=c(nl,nb))
  B <- t(A)%*%Bf
  
  # Create the PTDF matrix
  S <- matrix(0,nl,nb)
  S[,noslack] = as.matrix(Bf[,noslack]) %*% inv(as.matrix(B[noslack,noslack]))
  
  return(S)
}


data <- get.data(busdat,gendat,branchdat)
S <- make.PTDF(busdat,branchdat)

# Draw random samples
d = c()
for(i in 1:nbus) d = append(d, runif(1,dmin[i],dmax[i]))

while(sum(d) - sum(gmin) < 0) 
{
  d = c()
  for (i in 1:nrow(busdat)) d = append(d, runif(1,dmin[i],dmax[i]))
}

g_ib = ((sum(d)-sum(gmin))/sum(gmax-gmin))*(gmax-gmin)
g = g_ib + gmin

# Compute power injections and flows
p = g - d
f <- S %*% p


# Success/Failure
# fail = 0
# success = 0
# for(l in 1:nbr)
# {
#   if(inside.range(abs(f[l]), c(0,flim[l])))
#   {
#     success = success+1
#   }
#   else
#   {
#     print(l)
#     fail = fail+1
#   }
# }

