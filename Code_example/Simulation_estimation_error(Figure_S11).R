
#setwd("transfer_learning-main_final")
path <- getwd()
source(paste0(path,"/Function/Functions_rankestimation_simul.R")) # load r function


########################### Hetergoeneous design setting ###########################################
## generating target data
tgarget_n050p30r5q20 <- generator_target_simulsourcedetection(rep=100, n=50, p=30, q=20, r= 5, corx=0.5,cory=0)

# generating source datasets
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(60,8),hvec = c(5,5,5,20,20,20,30,30), 
                                                          p=30,q=20,ranksourcevec=rep(5,8),
                                                          B=tgarget_n050p30r5q20$B,cory=0,numsource=5)

# conduct simulation analysis
testrep1to3_het <- rep_forsimul(Ylist=tgarget_n050p30r5q20$Ylist,
                                Xlist=tgarget_n050p30r5q20$Xlist,
                                auxYlist_list=Sourceset_h5$auxYlist_list,
                                auxXlist_list=Sourceset_h5$auxXlist_list, repstart=1, repend=3, 
                                Btrue=tgarget_n050p30r5q20$B)


# Summary of source detection results
## truesupport: informative set
Summary_ft_detection_het(testrep1to3_het, truesupport=c(1,2,3))
# Summary of estimation results
Summary_ft_est(testrep1to3_het)


data <- data.frame(
  n_k = rep(c(20, 50), each=3*5),
  h = rep(rep(c(5,12,15), each=5), times=2),
  Method = rep(c("NR","Pooled-NR","[K]-Trans","MSD-Trans-NR","FSD-Trans-NR"), times=6),
  EE = mean(Summary_ft_est(testrep1to3_het)),
  SD = sd(Summary_ft_est(testrep1to3_het))
)

data$n_k <- factor(data$n_k, levels=c(20,50), labels=c(expression(n[k]==20),expression(n[k]==50)))
data$h <- factor(data$h)

save("data", file="f_S11.RData")
