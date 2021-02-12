
Getting Started
---------------

    library(CTGEIM4)

Generating time-course gene expression data
-----------------------------

    library(mvtnorm)

    n=100;m=10
    TT=seq(0,1,length=m)
    Y=mu=mu1=matrix(0,n,m)
    x1=rnorm(n);x2=rbinom(n,1,.1)
    b=rnorm(n)
    K=6
    index1=1:round(n/K)
    index2=(round(n/K)+1):(round(2*n/K))
    index3=(round(2*n/K)+1):round(3*n/K)
    index4=(round(3*n/K)+1):round(4*n/K)
    index5=(round(4*n/K)+1):round(5*n/K)
    index6=(round(5*n/K)+1):n

    cluster=rep(0,n)
    cluster[index1]=1
    cluster[index2]=2
    cluster[index3]=3
    cluster[index4]=4
    cluster[index5]=5
    cluster[index6]=6


    ar1_cor <- function(n, rho) {
      exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                        (1:n - 1))
      rho^exponent
    }
    R=ar1_cor(m,.5)

    for(i in index1){
      for(j in 1:m){
        mu[i,j]<-2*cos(2*pi*TT[j])}}
    for(i in index1){
      Y[i,]=rmvnorm(1,mu[i,],R)}


    for(i in index2){
      for(j in 1:m){
        mu[i,j]<-sqrt(2)*sin(4*pi*TT[j])}}
    for(i in index2){
      Y[i,]=rmvnorm(1,mu[i,],R)}


    for(i in index3){
      for(j in 1:m){
        mu[i,j]<-(TT[j]+.75)^(-3)-1}}
    for(i in index3){
      Y[i,]=rmvnorm(1,mu[i,],R)}


    for(i in index5){
      for(j in 1:m){
        mu[i,j]<--6*cos(2*pi*TT[j])}}
    for(i in index5){
      Y[i,]=rmvnorm(1,mu[i,],R)}


    for(i in index4){
          for(j in 1:m){
    mu[i,j]<--3*(2*TT[j]-1)}}
    for(i in index4){
      Y[i,]=rmvnorm(1,mu[i,],R)}

    for(i in index6){
      for(j in 1:m){
        mu[i,j]<-dbeta(TT[j],2,6)}}
    for(i in index6){
      Y[i,]=rmvnorm(1,mu[i,],R)}
  
  
  
### head(mExpression)
    |    | 1       | 2       | 3       | 4       | 5       | 6       | 7       | 8       | 9       | 10      |
    |----|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
    | 1  | 0.0178  | 0.0679  | -1.069  | -1.2376 | -0.0535 | -0.922  | 0.140   | 1.426   | 1.161   | 0.205   |
    | 2  | 0.1033  | 0.5072  | 0.252   | 1.7608  | 1.0259  | -1.538  | -0.948  | -1.362  | -0.309  | -1.401  |
    | 3  | -0.3261 | -0.1773 | 1.109   | 0.7118  | -0.0230 | 0.335   | -0.513  | -0.170  | -0.293  | 0.109   |
    | 4  | 2.1263  | 0.2035  | 0.400   | 0.0404  | 0.3727  | -0.630  | -2.042  | -0.257  | -0.861  | -1.571  |
    | 5  | -0.7039 | 0.1925  | -0.839  | 0.2509  | 0.1742  | 0.718   | -0.639  | 0.744   | -0.258  | 0.43    |
    
### runing CTGE_IM4 function

> For running CTGE_IM4, we consider 6 as maximum number of clusters 
> (maxcluster=6), one chains (nchains=1), the lengths of
> MCMC iteration and burn-in are set at 1000 and 500, respectively
> (niter=1000, burnin=500). Also, we consider
> tau0=1 and varsigma0=1  as initial values. 


     RESULTS=CTGE_IM4(mExpression=Y, maxcluster=10, niter=1000, nburnin=500, nthin=2, 
     nchains=1, tau0=1, varsigma0=1, genename=1:length(Y[,1]))
     
###  Cluster_Plot

> The time-course profiles for the genes of each cluster which detected by CTGE-IM4. 


    time=1:length(Y[1,])
    Cluster_Plot(mExpression=Y,time,RESULTS$MPEAR,Some=TRUE)
    !(Cluster.png)


###  PSM_Plot

> The complex Heatmap for the posterior similarity matrix of CTGE-IM4.

 
    PSM_plot(RESULTS$PSM)
    !(PSM.png)

