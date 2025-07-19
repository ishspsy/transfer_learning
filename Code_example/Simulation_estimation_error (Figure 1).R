
#setwd("transfer_learning-main_final")
path <- getwd()
source(paste0(path,"/Function/Functions_rankestimation_simul.R")) # load r function


##### Simulation model with n0=120, p=80, q=50, r=8, h=10


### Simulating target data
# rep: the number of replicates
# n=n0, r=rank of B
# corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
# cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
targestset <- generator_target_normal(rep=200,n=50,p=30,q=20,r=5, corx=0.5,cory=0)


### Simulating 100 source datasets
## rank source= rank of contrast
# corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
# cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
# generating 100 replicates for a source data with nk=300, transferring level=2, rank_source=5
sources <-generator_source_rankestimationsimul(rep=200, n=60,p=30,q=20,B=targestset$B,
                                               corx=0.5,cory=0,
                                               h=9,numsource=1, ranksource =5)

### Fitting algorithms to the first five replicates

## Ah-Trans-NR
# Ylist: contains the simulated Y
# Xlist: contains the simulated X
# auxYlist_list: contains the simulated auxYlist, wheer auxYlist is the list containing source response matrices
# auxXlist_list:  contains the simulated auxXlist, wheer auxYlist is the list containing source design matrices
# repstart, repend: simulation results are obtaeind from repstart to repend, 
#  e.g., when two values are set as 1 and 100, simulation results are obtained from the first 100 simulated data

# lamseq_w: candiates of tuning parameters for the firs step
# lamseq_delta: candiates of tuning parameters for the bias-correction step

# maxiter_inital: the maxium of the number of iterations in the ADMM algorithm for computing initial estimate
# tol_inital: tolerance error in the ADMM algorithm for computing initial estimate
# maxiter_biascorrection: the maxium of the number of iterations in the ADMM algorithm for computing contrast
# tol_inital: tolerance error in the ADMM algorithm for computing contrast
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)


Lassoerr_list = NULL
PooledNRerr_list = NULL
Twosteperr_list = NULL
for (ii in 1:200){
## NR
### Five-fold CV for NR
lamLassoseq <- c(0.001, 0.01, 0.05, 0.1,0.15,0.2, 0.3,0.35, 0.4, 0.5, 0.6, 0.7,  0.8, 0.9, 1, 1,1, 1.2, 1.5, 1.8, 2, 2.2 ,2.5)
# lamseq: candidates of lambda in algorithm in Section S2.1
cv.nuclear(Y=targestset$Ylist[[ii]], 
           X=targestset$Xlist[[ii]], 
           B=NULL,
           L=NULL,
           eta=1,lamseq = lamLassoseq,
           tol=1e-04,maxiter=100) -> cvLasso 
p <- ncol(targestset$Xlist[[ii]])
q <- ncol(targestset$Ylist[[ii]])
dLasso <- length(which(svd(matrix(cvLasso$B[-1,],p,q))$d > 10^(-6)))
Lassoerr <- norm(cvLasso$B[-1,] - targestset$B, "F") / sqrt(p*q)
Lassoerr_list= c(Lassoerr_list, Lassoerr)


### Five-fold CV for Pooled-NR
lamseq_w <- lamLassoseq
# auxYlist: a list contaring resposne matrix for source datasets
# auxXlist: a list contaring coviaraites matrix for source datasets
# lamseq: candidates of lambda in algorithm in Section S2.1
cv.pooledNR(Y=targestset$Ylist[[ii]], 
            X=targestset$Xlist[[ii]], 
            auxYlist = sources$auxYlist_list[[ii]],
            auxXlist = sources$auxXlist_list[[ii]],  
            B=NULL,
            L=NULL,eta=1,lamseq = lamLassoseq,
            tol=1e-04,maxiter=100) -> PooledNRres

dPooledNR <- length(which(svd(matrix(PooledNRres$B[-1,],p,q))$d > 10^(-6))) #  
PooledNRerr <- norm(PooledNRres$B[-1,] - targestset$B, "F") / sqrt(p*q)
PooledNRerr_list= c(PooledNRerr_list, PooledNRerr)


### Five-fold CV for [K]-Trans
lam_delta <- c(lamseq_w,3)
# auxYlist: a list contaring resposne matrix for source datasets
# auxXlist: a list contaring coviaraites matrix for source datasets
# lamseq_w: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
# lamseq_delta: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
cv.twostep(Y=targestset$Ylist[[ii]], 
           X=targestset$Xlist[[ii]], 
           B=NULL, L=NULL,eta=1,
           lamseq_w = lamseq_w,
           lamseq_delta=lam_delta,
           tol=1e-04,maxiter=100,
           auxYlist = sources$auxYlist_list[[ii]],
           auxXlist = sources$auxXlist_list[[ii]]) -> TwosteptransRes

dtwostep <- length(which(svd(matrix(TwosteptransRes$B[-1,],p,q))$d > 10^(-6)))  
Twosteperr <- norm(TwosteptransRes$B[-1,] - targestset$B, "F") / sqrt(p*q)
Twosteperr_list = c(Twosteperr, Twosteperr)
}


## Repeat the above settings by changing h and |A_h| values  

data <- data.frame(
  A_h = rep(c(5, 10), each = 4*3),
  h = rep(c(9, 12, 15, 18), times = 2, each = 3),
  Method = rep(c("NR", "A[h]-Pooled-NR", "A[h]-Trans-NR"), times = 8),
  EE = c(mean(Lassoerr_list),mean(PooledNRerr_list), mean(Twosteperr_list)),
  SD = c(sd(Lassoerr_list), sd(PooledNRerr_list), sd(Twosteperr_list))
  )

data$h <- factor(data$h, levels = c(9, 12, 15, 18))
data$Method <- factor(data$Method, levels = c("NR", "A[h]-Pooled-NR", "A[h]-Trans-NR"))
data$A_h <- factor(data$A_h, levels = c(5, 10), labels = c(expression(paste("|", A[h], "| = 5")), expression(paste("|", A[h], "| = 10"))))

save("data", file="f_1.RData")
