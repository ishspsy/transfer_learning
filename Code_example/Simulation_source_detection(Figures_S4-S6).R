


#setwd("transfer_learning-main_final")
path <- getwd()
source(paste0(path,"/Function/Functions_naiveapproaches.R")) # load r function
source(paste0(path,"/Function/Functions_MSDtrans.R")) # load r function
source(paste0(path,"/Function/Functions_FSDtrans.R")) # load r function
source(paste0(path,"/Function/Functions_nr-basedmethods_simul.R")) # load r function





########################### Hetergoeneous design setting ###########################################
# Case 1
tic()
n0 = 50; p=30
tgarget_n050p30r5q20 <- generator_target_simulsourcedetection(rep=100, n=n0, p=p, q=20, r= 5, corx=0.01,cory=0)
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(50,5),hvec = rep(9,5), 
                                                           p=30,q=20,ranksourcevec=sample(3:12, size = 5, replace = TRUE),
                                                           B=tgarget_n050p30r5q20$B,cory=0,numsource=5)
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
toc()

# Case 2
tic()
n0=100; p=500
tgarget_n0100p500r5q20 <- generator_target_simulsourcedetection(rep=100, n=n0, p=p, q=20, r= 5, corx=0.01,cory=0)
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(50,5),hvec = rep(9,5), 
                                                          p=30,q=20,ranksourcevec=sample(3:12, size = 5, replace = TRUE),
                                                          B=tgarget_n050p30r5q20$B,cory=0,numsource=5)
testrep1to3_het <- rep_forsimul(Ylist=tgarget_n0100p500r5q20$Ylist,
                                Xlist=tgarget_n0100p500r5q20$Xlist,
                                auxYlist_list=Sourceset_h5$auxYlist_list,
                                auxXlist_list=Sourceset_h5$auxXlist_list, repstart=1, repend=3, 
                                Btrue=tgarget_n0100p500r5q20$B)


# Summary of source detection results
## truesupport: informative set
Summary_ft_detection_het(testrep1to3_het, truesupport=c(1,2,3))
# Summary of estimation results
Summary_ft_est(testrep1to3_het)
toc()

df <- data.frame(
  n0 = c(50, 50, 100, 100, 50, 50, 100, 100),
  p = c(30, 30, 500, 500, 30, 30, 500, 500),
  h = c(5, 9, 30, 50, 5, 9, 30, 50),
  Method = c("FSD-Trans-NR", "FSD-Trans-NR", "FSD-Trans-NR", "FSD-Trans-NR",
             "MSD-Trans-NR", "MSD-Trans-NR", "MSD-Trans-NR", "MSD-Trans-NR"),
  PCI = Summary_ft_detection_het(testrep1to3_het, truesupport=c(1,2,3)),
  A_hat = mean(Summary_ft_est(testrep1to3_het)),
  SD = sd(Summary_ft_est(testrep1to3_het))
)

df <- df %>% mutate(Scenario = paste0("n0=", n0, ", p=", p, ", h=", h))

save("df", file="f_S4.RData")
