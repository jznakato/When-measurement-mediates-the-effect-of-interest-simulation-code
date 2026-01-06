## ==================================================================================================
## generate.cluster: Function to generate data for one cluster. 
## Input
# - "dgp":The "type" of missingness  generating process i.e simple, complex, lb
# - effect ; Whether treatment effect is TRUE or FALSE
# - N : The cluster size i.e number of individuals within each cluster
# - j : jth cluster. j=1,2...J
# - verbose: whether to print summary stats and histograms for propensity scores (TRUE or FALSE)
## output: A dataset containing individual level variables and cluster-level variables for each cluster
## Update from original manuscript :
# Changed get.Y2 function by reducing the E1 in response to reviewer
# Added a simple get.delta.simple function for a simpler missingness mechanism upon reviewer's request
#======================================================================================================
generate.cluster <- function(dgp, effect, N, j, verbose){
  
# N is the number of individuals in cluster j 
  UE1 <- rep(runif(1, 0, 1), N)
  UE2 <- rep(runif(1, 0, 0.5), N)
  
## cluster-level covariates (Ec's)
  E1 <- rep(abs(rnorm(1, UE1, 1)), N) 
  E2 <- rep(abs(rnorm(1,UE2, 1)), N)
  
### Obtain fW1, fW2, fW3
## Individual level covariates
  W1 <-  as.numeric( round(scale((runif(N, 18, 60))),2))  ## individual level covariates (age) scaled  
  W2 <- rbinom(N,1,UE1)
  W3<- rbinom(N, 1, UE2) 
  
  
## Measurement Indicator
  U.Delta <- runif(N, 0,1)
  
  if(effect){
  
  if(dgp=='simple'){
    Delta.1 <- get.delta.simple(A=1, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                                UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    Delta.0 <- get.delta.simple(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                                UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
  } else if(dgp=='complex'){
    Delta.1 <- get.delta.complex(A=1, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                                 UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    Delta.0 <- get.delta.complex(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                 UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
  } else if(dgp=='lb'){
    Delta.1 <- get.delta.lb(A=1, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                            UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    Delta.0 <- get.delta.lb(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                            UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
  }  
  }else{
    
    if(dgp=='simple'){
     
      Delta.0 <- get.delta.simple(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                                  UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    } else if(dgp=='complex'){
     
      Delta.0 <- get.delta.complex(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                   UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    } else if(dgp=='lb'){
   
      Delta.0 <- get.delta.lb(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3, 
                              UE1=UE1, UE2=UE2, U.Delta=U.Delta, verbose=verbose)
    }
    Delta.1<-Delta.0
  }
  

## underlying HIV  status (Y1) 
# **** no effect of A 
  U.Y1<-  runif(N,0, 1)
  
  p.Y1 <- plogis(-0.8 + 0.3*E1 + 0.15*E2 + 2*W1  + 1.5*W2 - 3.0*W3 - 0.55*UE1 + 0.15*UE2) #original manuscript
 
  
  if(verbose) { print(summary( p.Y1)); hist(p.Y1) }
  Y1_star <- as.numeric(U.Y1 >p.Y1 ) ## have a  1-p.Y1 probability of being greater than 

## comment (01/23/2025)
# Y1 is generated as binary indicator of having HIV i.e AT risk is Y1=0 . 
# However, In manuscript Y1=1 is indicator of being at risk. 
  
#  generate outcome Y2 (uptake of PREP/PEP)
  U.Y2 <- runif(N, 0, 1)
  
  Y2.0 <- get.Y2(A=0, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2, UE1=UE1, UE2=UE2,
                 Y1_star=Y1_star, Delta=Delta.0, U.Y2=U.Y2, verbose=verbose)
  if(effect){
    Y2.1<- get.Y2(A=1, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2, UE1=UE1, UE2=UE2,
                  Y1_star=Y1_star, Delta=Delta.1, U.Y2=U.Y2, verbose=verbose)
  }else{
# if under the null, then set the counterfactual outcome Y2.1 to Y2.0
    Y2.1 <- Y2.0
  }
  
  
  dt<-data.frame( cbind(id=j, n=N, alpha=1, U=1, W1 = W1, W2 = W2, W3=W3,
                        E1 = E1, E2 = E2,  Delta.1= Delta.1, Delta.0= Delta.0, 
                        Y1_star=Y1_star, Y2.1=Y2.1, Y2.0=Y2.0))
  dt
}



get.delta.lb <- function(A,E1,E2,W1,W2,W3, UE1, UE2, U.Delta, verbose=F) {
   pscore_delta <- plogis(0.1 + 2*A- 2*W1  + 0.99*W2*A +0.75*W3 - 0.75*E1 + 0.85*E2 - 0.5*UE1 +0.5*UE2 ) 
  if(verbose) { print(summary(pscore_delta)); hist(pscore_delta) }
  Delta <- as.numeric(U.Delta < pscore_delta) 
## is the same way of saying probability of measurement is given by pscore_delta
   Delta
  
}

## For each cluster, get.delta is a simpler complex linear function with main terms and no interactions
get.delta.simple<- function(A, E1, E2, W1, W2, W3, UE1, UE2, U.Delta, verbose=F) {

  pscore_delta <- A*plogis(1 + 0.5*W1  - 1.5*W2  + 0.3*W3 + 2*UE1 ) + (1-A)*plogis(-0.1 + 0.5*W1 +1.05*W2  - 0.3*W3 - 2*UE1 )
  if(verbose) { print(summary(pscore_delta)); hist(pscore_delta) }
  Delta <- as.numeric(U.Delta < pscore_delta) ##is the same way of saying probability of measurement is given by pscore_delta
  
  Delta
}

## get.delta is a more complex non-linear function
get.delta.complex <- function(A,E1,E2,W1,W2,W3,UE1, UE2, U.Delta, verbose=F) {
## generate delta such that A has a strong impact on delta
 
 pscore_delta <- W3*plogis(0.1+ 0.8*A- 2*W1  + 0.8*W2*A +0.75*W3 - 0.75*E1 +0.75*E2 - 0.15*UE1 +0.15*UE2 )+
   (1-W3)*plogis(0.1 + 0.8*A + 2*W1  - 0.8*W2*A - 0.75*W3 + 0.75*E1 - .75*E2 + 0.15*UE1 - 0.15*UE2 ) 


if(verbose) { print(summary(pscore_delta)); hist(pscore_delta) }
  Delta <- as.numeric(U.Delta < pscore_delta) 

  Delta

}



# misisngness mechanism depends on unobserved UEs but we still take care of this using 2 stage tmle
# UEs are confounders of Y1*star and Delta . in a sinlge stage we would not have identification

get.Y2 <- function(A, W1, W2, W3, E1, E2, UE1, UE2, Y1_star, Delta, U.Y2, verbose=F){
  
  
# calculating the probability of starting PrEP once one links (if HIV-)
  
p.Y2 <- plogis(0.2 + 0.08*A + 0.5*W1 + 0.4*W2 -1*E1 + .87*E2 - 0.08*A*W3 - 0*UE1 +0*E2 )  ##reduce the effect of E1 (main review)


  if(verbose) { print(summary( p.Y2)); hist(p.Y2) }
  
  t<- case_when (Y1_star==1 & Delta==1~0, # HIV+ (known) and linked
                 Y1_star==1 & Delta==0~0, # HIV+ (unknown) but didnt link
                 Y1_star==0 & Delta==0~0, # HIV- (unknown) but didnt link
                 Y1_star==0 & Delta==1 ~ as.numeric(U.Y2 < p.Y2))  # HIV- (known) AND linked
  
  t
  
}


## =====================================================================================================
## get.full.data : Generates individual data within each cluster and repeats the process for J clusters
## input : 
# dgp: Data generating process for delta
# J : Number of clusters
# N.mean : Mean cluster size
# N.sd : sd of cluster-sizes
# effect: TRUE or FALSE
## Output: Full dataset including counterfactuals
## =====================================================================================================
get.full.data <- function(dgp, J=64, N.mean=100, N.sd=20, effect=effect, verbose=F){
  
  n.indv<- round(rnorm(J, N.mean, N.sd))
  if(N.mean >=30){
    n.indv[n.indv < 30] <- 30
  } else{
    n.indv[n.indv < 10] <- 10
  }
  
  full.data<- NULL
  
  for(j in 1:J){
    data.comm.j <- generate.cluster(dgp=dgp, effect=effect, N=n.indv[j], j=j, verbose=verbose) 
    full.data <- rbind(full.data, data.comm.j)
  }
  
  
  X<- aggregate(x = full.data, by = list(id = full.data$id), FUN = mean)
  
## randomize observed treatment

  A <- c(rep(1, J/2), rep(0, J/2) )
  X$A <- sample(A)
  
## TRANSLATE TREATMENT TO FULL DATA
  full.data$A = X$A[match(full.data$id, X$id)]
##  observed  delta
  full.data$Delta<- ifelse(full.data$A == 1, full.data$Delta.1, full.data$Delta.0)
## set observed Y1 here 
  full.data$Y1 <- full.data$Y1_star*full.data$Delta
  
## observed Y2
full.data$Y2 = ifelse(full.data$A == 1, full.data$Y2.1, full.data$Y2.0)
  
  
full.data
}


## =====================================================================================================
## get.truth: Generates the true individual level and cluster-level effects
## =====================================================================================================
get.truth <- function(dgp, J, N.mean=200, N.sd=20, effect,verbose=F){
  
full.data<-get.full.data(dgp=dgp, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)

# individual-level effects (not needed right now)
  i1 <- mean(full.data$Y2.1)/mean(1- full.data$Y1_star)
  i0 <- mean(full.data$Y2.0)/mean(1- full.data$Y1_star)
  
# cluster-level effects 
## calculate true counterfactual proportion Yc within each cluster
# aggregate to the data to the cluster-level 
  truth.C <- aggregate(full.data[,c('U', 'A', 'Y2.1', 'Y2.0','Y1_star')], by=list(full.data$id), mean)
  
  Yc1 <- truth.C$Y2.1/(1-truth.C$Y1_star)
  Yc0 <- truth.C$Y2.0/(1-truth.C$Y1_star)
  c1 <- mean(Yc1)
  c0 <- mean(Yc0)
  
  truth<-data.frame( cbind(
    int.ind = i1,
    con.ind = i0, 
    RD.ind = i1 - i0,
    RR.ind= i1/i0,
    int.clust = c1, 
    con.clust = c0,
    RD.clust= c1 - c0,
    RR.clust= c1/c0,
    Psi_den = mean(1-truth.C$Y1_star),
# adding in measures of within/between variability
    k0 = sd(Yc0)/mean(Yc0),
    k1 = sd(Yc1)/mean(Yc1)
  ))
  truth
# return summary data.frame AND the cluster-specific summaries
  return(list(truth=truth, Yc=truth.C))
  
}

## =====================================================================================================
## get.obs: gets observed data including observed Y1 
## =====================================================================================================
get.obs <- function(O){
  
  # set observed outcome Y as (measurement indicator)x(underlying outcome)
  O$Y1<- case_when(O$Delta==0~NA,
                   O$Delta==1~O$Y1_star*O$Delta)
  
  
  O
}


