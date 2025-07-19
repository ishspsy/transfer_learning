
## For source detection simulation results (2, Figure S1, Figure S2, Figure S4, Figure S5)
library(tictoc)


#setwd("transfer_learning-main_final")
path <- getwd()
source(paste0(path,"/Function/Functions_naiveapproaches.R")) # load r function
source(paste0(path,"/Function/Functions_MSDtrans.R")) # load r function
source(paste0(path,"/Function/Functions_FSDtrans.R")) # load r function
source(paste0(path,"/Function/Functions_nr-basedmethods_simul.R")) # load r function



########################### Hetergoeneous design setting ###########################################
# set h, K
tic()
h = 2; K=4

tgarget_n050p30r5q20 <- generator_target_simulsourcedetection(rep=100, n=120, p=80, q=20, r= 5, corx=0.5,cory=0)
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(300,4),hvec = c(h,20,20,50), 
                                                           p=30,q=20,ranksourcevec=rep(5,4),
                                                           B=tgarget_n050p30r5q20$B,cory=0,numsource=K)
testrep1to3_het <- rep_forsimul(Ylist=tgarget_n050p30r5q20$Ylist,
                            Xlist=tgarget_n050p30r5q20$Xlist,
                            auxYlist_list=Sourceset_h5$auxYlist_list,
                            auxXlist_list=Sourceset_h5$auxXlist_list, repstart=1, repend=3, 
                            Btrue=tgarget_n050p30r5q20$B)


# Summary of source detection results
## truesupport: informative set
Summary_ft_detection_het(testrep1to3_het, truesupport=c(1))
# Summary of estimation results
Summary_ft_est(testrep1to3_het)
toc()


########################### Hetergoeneous design setting ###########################################
# set h, K
tic()
h = 2; K=7

tgarget_n050p30r5q20 <- generator_target_simulsourcedetection(rep=100, n=120, p=80, q=20, r= 5, corx=0.5,cory=0)
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(300,7),hvec = c(h,h,h,h,20,20,50), 
                                                          p=30,q=20,ranksourcevec=rep(5,7),
                                                          B=tgarget_n050p30r5q20$B,cory=0,numsource=K)
testrep1to3_het <- rep_forsimul(Ylist=tgarget_n050p30r5q20$Ylist,
                                Xlist=tgarget_n050p30r5q20$Xlist,
                                auxYlist_list=Sourceset_h5$auxYlist_list,
                                auxXlist_list=Sourceset_h5$auxXlist_list, repstart=1, repend=3, 
                                Btrue=tgarget_n050p30r5q20$B)



# Summary of source detection results
## truesupport: informative set
Summary_ft_detection_het(testrep1to3_het, truesupport=c(1,2,3,4))
# Summary of estimation results
Summary_ft_est(testrep1to3_het)
toc()



df <- data.frame(
  Ah = rep(c(1,1,1,4,4,4), each=2),
  Ah_c = rep(3,12),
  h = rep(c(2,5,10,2,5,10), each=2),
  Method = rep(c("FSD-Trans-NR","MSD-Trans-NR"), times=6),
  PCI = Summary_ft_detection_het(testrep1to3_het, truesupport=c(1,2,3,4)),
  A_hat = mean(Summary_ft_est(testrep1to3_het)),
  SD = sd(Summary_ft_est(testrep1to3_het))
)

df$Scenario <- paste0("|Ah|=", df$Ah, ", h=", df$h)
save("data", file="f_S1.RData")


