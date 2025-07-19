## plots
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)



# Figure 1
library(ggplot2)
library(dplyr)
library(tibble)
load("Results/f_1.RData")

ggplot(data, aes(x = h, y = EE, fill = Method)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = EE - SD, ymax = EE + SD), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ A_h, labeller = label_parsed) +
  labs(x = "h", y = "Estimation Error (EE)", fill = "Method", title = "Estimation Errors with Standard Deviations") +
  scale_fill_discrete(labels = c("NR", expression(A[h]*"-Pooled-NR"), expression(A[h]*"-Trans-NR"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


# Figure 2
library(ggplot2)
library(patchwork)
load("Results/f_2.RData")

p1 <- ggplot(data, aes(x = h, y = PCI, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_wrap(~n_k, labeller = label_parsed) +
  labs(x = "h", y = "PCI (%)", fill = "Method",
       title = "PCI (%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14))

p2 <- ggplot(data, aes(x = h, y = A_hat_mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = A_hat_mean - A_hat_sd, ymax = A_hat_mean + A_hat_sd),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~n_k, labeller = label_parsed) +
  labs(x = "h", y = expression(paste("Estimated ", "|", hat(A), "|")),
       fill = "Method",
       title = expression(paste("Estimated ", "|", hat(A), "|"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14))

p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")







### Figure S1
library(tidyverse)
library(gridExtra)
load("Results/f_S1.RData")

p1 <- ggplot(df, aes(x=Scenario, y=PCI, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  ylim(0,100) +
  labs(title="Percentage of Correct Identification (PCI)",
       x="Scenario",
       y="PCI (%)") +
  theme_bw(base_size=15) +
  theme(axis.text.x = element_text(angle=15, hjust=1),
        legend.title=element_blank(),
        plot.title=element_text(size=17, face="bold"))

p2 <- ggplot(df, aes(x=Scenario, y=A_hat, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=A_hat-SD, ymax=A_hat+SD),
                width=0.2, position=position_dodge(0.8)) +
  labs(title=expression("Estimated Set Size ("~"|"*hat(A)*"|"~")"),
       x="Scenario",
       y=expression("|"*hat(A)*"|")) +
  theme_bw(base_size=15) +
  theme(axis.text.x = element_text(angle=15, hjust=1),
        legend.title=element_blank(),
        plot.title=element_text(size=17, face="bold"))

grid.arrange(p1, p2, ncol=2)



# Figure S2
library(ggplot2)
library(gridExtra)
load("Results/figure_S2.RData")

ggplot(data_long, aes(x = factor(h), y = Error, group = Method, 
                      color = Method, linetype = Method)) +
  geom_line(size = 1.1) + 
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Error - SE, ymax = Error + SE), width = 0.2, size = 1) +
  facet_wrap(~ Sources, scales = "free_x") +
  theme_bw(base_size = 16) +
  labs(
    x = "h",
    y = "Estimation Error",
    title = "Estimation Error with Standard Deviations by Method and Scenario"
  ) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )





##Figure S3
library(ggplot2)
load("Results/f_S3.RData")

ggplot(df, aes(x=Scenario, y=Time, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=Time - SD, ymax=Time + SD), 
                width=0.2, position=position_dodge(0.8)) +
  labs(title="Computational Times by Method and Scenario",
       x="Scenario (p, |Ah|)",
       y="Computational Time (minutes)") +
  theme_bw(base_size=16) +
  theme(axis.text.x = element_text(angle=15, hjust=1, size=14),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=14))



### Figure S4
library(ggplot2)
library(gridExtra)
load("Results/f_S4.RData")


p1 <- ggplot(df, aes(x=Scenario, y=PCI, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  ylim(0,100) +
  labs(y="Percentage of Correct Identification (PCI)", x="Scenario",
       title="Percentage of Correct Source Identification") +
  theme_bw(base_size=15) +
  theme(axis.text.x = element_text(angle=20, hjust=1),
        legend.title=element_blank(),
        plot.title=element_text(size=18, face="bold"))

p2 <- ggplot(df, aes(x=Scenario, y=A_hat, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  geom_errorbar(aes(ymin=A_hat-SD, ymax=A_hat+SD),
                position=position_dodge(width=0.8), width=0.2) +
  labs(x="Scenario",
       y=expression(paste("Estimated ", "|", hat(A), "|")),
       title="Estimated Size of Source Set") +
  theme_bw(base_size=15) +
  theme(axis.text.x = element_text(angle=20, hjust=1),
        legend.title=element_blank(),
        plot.title=element_text(size=18, face="bold"))

library(gridExtra)
grid.arrange(p1, p2, ncol=2)



# Figure S5
### 
library(ggplot2)
load("Results/f_S5.RData")

ggplot(df, aes(x = Method, y = EE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
  geom_errorbar(aes(ymin = EE - SD, ymax = EE + SD),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Scenario, scales = "free_y") +
  labs(title = "Estimation Errors by Scenario and Method",
       y = "Estimation Error",
       x = "Method") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_brewer(palette = "Set2")




##Figure S6
library(ggplot2)
load("Results/f_S6.RData")

ggplot(df, aes(x=Scenario, y=Time, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=Time - SD, ymax=Time + SD), width=0.2,
                position=position_dodge(0.8)) +
  labs(title="Mean Elapsed Times by Method and Scenario",
       x="Scenario (n0, p, h)",
       y="Elapsed Time (minutes)") +
  theme_bw(base_size=15) +
  theme(axis.text.x=element_text(angle=15, hjust=1, size=14),
        plot.title=element_text(size=18, face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=14))




### Figure S7
library(ggplot2)
library(gridExtra)
load("Results/f_S7.RData")

p1 <- ggplot(df, aes(x=factor(h), y=PCR, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  facet_wrap(~paste("True Rank =", rep(c(3,7), each=8)), scales="free_x") +
  ylim(0,100) +
  labs(x="h",
       y="PCR (%)",
       title="Percentage Correct Rank Recovery (PCR)") +
  theme_bw(base_size=15) +
  theme(legend.title=element_blank(),
        plot.title=element_text(size=18, face="bold"))

p2 <- ggplot(df, aes(x=factor(h), y=EE, color=Method, group=Method, linetype=Method)) +
  geom_line(size=1.1, position=position_dodge(width=0.5)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=EE-EE_sd, ymax=EE+EE_sd),
                position=position_dodge(width=0.5), width=0.2) +
  facet_wrap(~paste("True Rank =", rep(c(3,7), each=8)), scales="free") +
  labs(x="h", y="Estimation Error (Frobenius norm)",
       title="Estimation Error by Scenario and Method") +
  theme_bw(base_size=15) +
  theme(plot.title=element_text(size=18, face="bold"))

library(gridExtra)
grid.arrange(p1,p2,ncol=2)




# Figure S8
library(ggplot2)
library(patchwork)
load("Results/f_S8.RData")       

# Estimation Error plot
p1 <- ggplot(estimation_error_a, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Frobenius Norm"="dotted" #,
                                 # "Nuclear Norm + LASSO"  = "dotdash"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Frobenius Norm"="darkgreen"
                              #"Nuclear Norm + LASSO"  = "purple"
  )) +
  theme_bw() +
  labs(title="Estimation Error",
       x=expression(1/sigma^2),
       y="Estimation Error")  

# Estimation Error plot
p3 <- ggplot(estimation_error_b, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Frobenius Norm"="dotted"
                                 # "Nuclear Norm + LASSO"  = "dotdash"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Frobenius Norm"="darkgreen"
                              #"Nuclear Norm + LASSO"  = "purple"
  )) +
  theme_bw() +
  labs(title="Estimation Error",
       x=expression(1/sigma^2),
       y="Estimation Error")  #+


# Estimation Error plot
p5 <- ggplot(estimation_error_c, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Frobenius Norm"="dotted"
                                 #"Nuclear Norm + LASSO"  = "dotdash"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Frobenius Norm"="darkgreen"
                              #"Nuclear Norm + LASSO"  = "purple"
  )) +
  theme_bw() +
  labs(title="Estimation Error",
       x=expression(1/sigma^2),
       y="Estimation Error")  #+


library(patchwork)

final_plot <- (
  (p1) / 
    (p3) / 
    (p5) +
    plot_layout(guides = "collect")
) &
  theme_bw(base_size = 14) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )

final_plot




# Figure S9
library(ggplot2)
library(patchwork)
load("Results/f_S9.RData")  

# Estimation Error plot
p1 <- ggplot(estimation_error_a, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Nuclear Norm + Lasso"="dotted"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Nuclear Norm + Lasso"="darkgreen"
  )) +
  theme_bw() +
  labs(title=expression("Estimation Error when " * alpha == 0),
       x=expression(1/sigma^2),
       y="Estimation Error") 


# Estimation Error plot
p2 <- ggplot(estimation_error_b, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Nuclear Norm + Lasso"="dotted"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Nuclear Norm + Lasso"="darkgreen"
  )) +
  theme_bw() +
  labs(title=expression("Estimation Error when " * alpha == 0.25),
       x=expression(1/sigma^2),
       y="Estimation Error") 



# Estimation Error plot
p3 <- ggplot(estimation_error_c, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Nuclear Norm + Lasso"="dotted"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Nuclear Norm + Lasso"="darkgreen"
  )) +
  theme_bw() +
  labs(title=expression("Estimation Error when " * alpha == 0.5),
       x=expression(1/sigma^2),
       y="Estimation Error") 


# Estimation Error plot
p4 <- ggplot(estimation_error_d, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Nuclear Norm + Lasso"="dotted"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Nuclear Norm + Lasso"="darkgreen"
  )) +
  theme_bw() +
  labs(title=expression("Estimation Error when " * alpha == 0.75),
       x=expression(1/sigma^2),
       y="Estimation Error") 



# Estimation Error plot
p5 <- ggplot(estimation_error_e, aes(x=SNR, y=Error, 
                                     color=Method, linetype=Method, group=Method)) +
  geom_line(size=1) + 
  geom_point(size=3) +
  scale_linetype_manual(values=c("Nuclear Norm"="solid", 
                                 "LASSO"="dashed", 
                                 "Nuclear Norm + Lasso"="dotted"
  )) +
  scale_color_manual(values=c("Nuclear Norm"="blue",
                              "LASSO"="red",
                              "Nuclear Norm + Lasso"="darkgreen"
  )) +
  theme_bw() +
  labs(title=expression("Estimation Error when " * alpha == 1),
       x=expression(1/sigma^2),
       y="Estimation Error") 


library(patchwork)
row1 <- p1 + p2                                  # p1 | p2
row2 <- p3 + p4                                  # p3 | p4

row3 <- plot_spacer() + p5 + plot_spacer() +
  plot_layout(widths = c(0.25, 0.50, 0.25))

final_plot <- (row1 / row2 / row3) +             # 행 단위로 쌓기
  plot_layout(guides = "collect") &              # 범례 하나로 모으기
  theme_bw(base_size = 14) &
  theme(
    plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )

final_plot







# Figure S10
library(ggplot2)
library(patchwork)
load("Results/f_S10.RData")

p1 <- ggplot(estimation_error, aes(x=SNR, y=Error, color=Method, linetype=Method)) +
  geom_line(size=1) + geom_point(size=3) +
  labs(title="Estimation Error", y="Estimation Error", x=expression(1/sigma^2))

p2 <- ggplot(prediction_accuracy, aes(x=SNR, y=MSPE, color=Method, linetype=Method)) +
  geom_line(size=1) + geom_point(size=3) +
  labs(title="Prediction Accuracy (MSPE)", y="MSPE", x=expression(1/sigma^2))

p3 <- ggplot(tpr, aes(x=SNR, y=TPR, color=Method, linetype=Method)) +
  geom_line(size=1) + geom_point(size=3) +
  labs(title="Selection Accuracy (TPR)", y="True Positive Rate", x=expression(1/sigma^2))

p4 <- ggplot(fpr, aes(x=SNR, y=FPR, color=Method, linetype=Method)) +
  geom_line(size=1) + geom_point(size=3) +
  labs(title="Selection Accuracy (FPR)", y="False Positive Rate", x=expression(1/sigma^2))

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") &
  theme_bw(base_size = 16) &
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold")) &
  theme(legend.position = "bottom")





# Figure S11
library(ggplot2)
library(ggpattern)
load("Results/f_S11.RData")

ggplot(data, aes(x=h, y=EE, fill=Method, pattern=Method)) +
  geom_bar_pattern(stat="identity", position=position_dodge(0.9),
                   color="black", pattern_fill="black",
                   pattern_angle=45, pattern_density=0.1,
                   pattern_spacing=0.02,
                   pattern_key_scale_factor=0.6) +
  geom_errorbar(aes(ymin=EE-SD, ymax=EE+SD), width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~n_k, labeller = label_parsed) +
  labs(x="h", y="Estimation Error (EE)", fill="Method", pattern="Method",
       title="Estimation Errors with Standard Deviations") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        legend.position = "bottom")


