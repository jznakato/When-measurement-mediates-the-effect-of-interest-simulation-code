# Code by Joy Nakato 
# With review/editing by Laura Balzer


## =====================================================================================================
## get_full_data: Generates the "full" data for the cluster randomized trial (CRT)
# - gets individual data within each cluster and repeats the process for J clusters
# - "full" because counterfactuals are also returned
## Input: 
# dgp: Data generating process for delta
# J: Number of clusters
# N_mean: Mean cluster size
# N_sd: sd of cluster-sizes
## Output: 
# - Full dataset including counterfactuals
## =====================================================================================================
get_full_data <- function(dgp, J=50, N_mean=200, N_sd=10, verbose=F){
  
  # get cluster size 
  N_j <- round(rnorm(J, N_mean, N_sd))

  full_data<- NULL
  
  for(j in 1:J){
    # generate the data for each cluster j in 1:J
    data_j <- generate_cluster(dgp=dgp, N=N_j[j], j=j, verbose=verbose) 
    full_data <- rbind(full_data, data_j)
  }
  
  ## randomize intervention
  X<- data.frame(cbind(
                  id=1:J, 
                  A=sample( c(rep(1, J/2), rep(0, J/2) ) )
                  ))
  
  ## assign the treatment 
  full_data$A = X$A[match(full_data$id, X$id)]
  ##  observed delta
  full_data$Delta<- ifelse(full_data$A == 1, full_data$Delta_1, full_data$Delta_0)
  ## set observed Y1 here 
  full_data$Y1 <- full_data$Y1_star*full_data$Delta
  
  ## observed Y2
  full_data$Y2 = ifelse(full_data$A == 1, full_data$Y2_1, full_data$Y2_0)
  
  full_data
}



## ==================================================================================================
## generate_cluster: Function to generate data for one cluster. 
## Input
# - dgp: simple or complex 
# - N : The cluster size i.e number of individuals within each cluster
# - j : jth cluster. j=1,2...J
# - verbose: whether to print summary stats and histograms for propensity scores 
## Output: 
# - Dataset containing individual-level variables and cluster-level variables for each cluster
#======================================================================================================
generate_cluster <- function(dgp, N, j, verbose){
  
  # N is the number of individuals in cluster j 
  UE1 <- rep(runif(1, 0, 1), N)
  UE2 <- rep(runif(1, 0, 0.5), N)
  
  ## cluster-level covariates(Ecs)
  E1 <- rep(abs(rnorm(1, UE1, 1)), N) 
  E2 <- rep(abs(rnorm(1, UE2, 1)), N)
  
  ## Individual level covariates
  W1 <- runif(N, 2, 2)
  W2 <- rbinom(N, 1, UE1)
  W3 <- rbinom(N, 1, UE2) 
  
  ## Underlying indicator of being in the focus population
  ## For OPAL, underlying indicator of having HIV riks
  UY1_star<-  runif(N, 0, 1)
  pY1_star <- plogis(-0.8+0.3*E1+0.15*E2+2*W1+1.5*W2-3*W3-0.55*UE1+0.15*UE2)
  if(verbose) { print(summary(pY1_star)); hist(pY1_star) }
  Y1_star <- as.numeric(UY1_star < pY1_star ) 
  
  ## Indicator fo screening
  UDelta <- runif(N, 0, 1)

  if(dgp=='simple'){
      Delta_1 <- get_delta_simple(A=1, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                  UE1=UE1, UE2=UE2, UDelta=UDelta, verbose=verbose)
      Delta_0 <- get_delta_simple(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                  UE1=UE1, UE2=UE2, UDelta=UDelta, verbose=verbose)
  } else{
      Delta_1 <- get_delta_complex(A=1, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                   UE1=UE1, UE2=UE2, UDelta=UDelta, verbose=verbose)
      Delta_0 <- get_delta_complex(A=0, E1=E1, E2=E2, W1=W1,W2=W2, W3=W3,
                                   UE1=UE1, UE2=UE2, UDelta=UDelta, verbose=verbose)
  }  
  
  #  generate outcome Y2 (uptake of PrEP/PEP)
  UY2 <- runif(N, 0, 1)
  
  Y2_0 <- get_Y2(A=0, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2,
                 Y1_star=Y1_star, Delta=Delta_0, UY2=UY2, verbose=verbose)

  Y2_1<- get_Y2(A=1, W1=W1, W2=W2, W3=W3, E1=E2, E2=E2,
                  Y1_star=Y1_star, Delta=Delta_1, UY2=UY2, verbose=verbose)

  
  dt<-data.frame( cbind(id=j, n=N, E1 = E1, E2 = E2,  
                        W1 = W1, W2 = W2, W3=W3,
                        Delta_1= Delta_1, Delta_0= Delta_0, 
                        Y1_star=Y1_star, Y2_1=Y2_1, Y2_0=Y2_0))
  dt
}



# "Simple" because once we stratify on cluster in Stage 1,
# this reduces to a simple linear function with main terms. 
get_delta_simple<- function(A, E1, E2, W1, W2, W3, UE1, UE2, UDelta, verbose=F) {
  
  pscore_delta <- A*plogis(1.5+ 0.5*W1 - .5*W2 -2*W3) + 
                  (1-A)*plogis(-1 + 0.5*W1 + .5*W2 + 2*W3)
  if(verbose) { print(summary(pscore_delta)); hist(pscore_delta) }
  Delta <- as.numeric(UDelta < pscore_delta)
  
  Delta
}

# "Complex" because even after stratifying on clsuter in Stage 1,
# this remains non-linear
get_delta_complex <- function(A, E1, E2, W1, W2, W3, UE1, UE2, UDelta, verbose=F) {

  pscore_delta <- W3*plogis(0.85+0.8*A-2*W1+0.8*W2*A-0.75*E1+0.75*E2)+
    (1-W3)*plogis(0.1+0.8*A+2*W1-0.8*W2*A+0.75*E1-0.75*E2 ) 
  
  # 
  # pscore_delta <- W3*plogis(1 + 0.5*W1 - 1*W2 + 0.25*A) + 
  #                 (1-W3)*plogis(-0.5 - 0.5*W1 + 1*W2 + 0.25*A)
  
  if(verbose) { print(summary(pscore_delta)); hist(pscore_delta) }
  Delta <- as.numeric(UDelta < pscore_delta) 
  
  Delta
  
}


# NOTE: Generating Y2 as a function of Y1_star and Delta is equivalent to generating Y2
# as a function of observed Y1= Delta*Y1_star

get_Y2 <- function(A, W1, W2, W3, E1, E2, Y1_star, Delta, UY2, verbose=F){
  
  # calculating the probability of the outcome if eligible
  # OPAL - probability of starting PrEP if at risk
  pY2 <- plogis(0.2 + 0.08*A + 0.5*W1 + 0.4*W2 -1*E1 + .87*E2 - 0.08*A*W3)   
  if(verbose) { print(summary( pY2)); hist(pY2) }
  
  t<- case_when (Y1_star==0 & Delta==1~0, # screened but not in focus pop
                 Y1_star==0 & Delta==0~0, # not in focus pop & not screened 
                 Y1_star==1 & Delta==0~0, # in focus pop but not screened
                 Y1_star==1 & Delta==1 ~ as.numeric(UY2 < pY2))  # 
  
  t
  
}



## ======================================================================================
## get_truth: function to calculate the true values for the causal estimand
## ======================================================================================
get_truth <- function(dgp, J, N_mean=200, N_sd=10,verbose=F){

  full_data <- get_full_data(dgp=dgp, J=J, N_mean=N_mean, N_sd=N_sd,
                             verbose=F)

  # get the cluster-level counterfactual Yc(ac)= P(Y2(ac)=1 | Y1*=1)
  Yc_1 <- Yc_0 <- rep(NA, nrow=J)
  for(j in 1:J){
    # grab the jth cluster
    temp <- full_data[full_data$id==j, ]
    # subset on Y1*=1
    temp <- temp[temp$Y1_star==1, ]
    # now take the counterfactual means
    Yc_1 <- mean(temp$Y2_1)
    Yc_0 <- mean(temp$Y2_0)
  }
  
  truth <- data.frame(cbind(
    EYc_1 = mean(Yc_1),
    EYc_0 = mean(Yc_0)
  ))
  truth
}
