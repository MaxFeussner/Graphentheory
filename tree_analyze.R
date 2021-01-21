rm(list=ls())                                  # Clear Workspace to avoid data mistakes

list.of.packages <- c("rstudioapi",            # import find path script
                      "ggplot2",               # import plotting library
                      "ggpubr",                # Align plots in grid
                      "RColorBrewer",          # for beautiful colored 
                      "dplyr",                  # rename columns easely and more
                      "plotly",
                      "tidyverse"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(plotly)
library(tidyverse)

pathScript <- dirname(rstudioapi::getSourceEditorContext()$path) # den Pfad des Scriptes finden
setwd(pathScript)  # set Working directory

dataDir <- paste(pathScript, "01_Data", sep="/")
treeDataDf <- read.csv(paste(pathScript, 'Tree_data.csv', sep = '/'), sep = ',', dec = ".", header = TRUE, stringsAsFactors = F)

#i = 1
for (i in 1:length(treeDataDf$ID)) {
  treeDataDf$Group[i] <- strsplit(as.character(treeDataDf$ID[i]), "_")[[1]][1]
}
treeDataDf$Group <- as.factor(treeDataDf$Group)
names(treeDataDf)
levels(treeDataDf$Group)

# hist(as.numeric(treeDataDf$Fraction_of_Xenologs))

####  2.1 ####
# Spearman since data is not parametric # 'c("pearson", "kendall", "spearman")'
#### Preparation for SUMDF ####
# create an empty DF
sumDf <- data.frame(Gruppe = as.character(),
                    Duplication_Rate = as.numeric(),
                    Loss_Rate = as.numeric(),
                    HGT_Rate = as.numeric(),
                    Slope = as.numeric(),
                    Intercept = as.numeric(),
                    Spearman_Corr = as.numeric())
for (i in 1:length(levels(treeDataDf$Group))) {
  sumDf[i,1] <- 0  
}
######################
#### GENES VS HGT ####
######################
for (group in 1:length(levels(treeDataDf$Group))) {
  mod = lm(treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])] ~ 
             treeDataDf$Number_of_leaves_tgt[which(treeDataDf$Group == levels(treeDataDf$Group)[group])])
  
  coef <- cor(x = treeDataDf$Number_of_leaves_tgt[which(treeDataDf$Group == levels(treeDataDf$Group)[group])],
              y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])], method = 'spearman') 
  
  ### Plots ###
  png(paste("02_Plots/", levels(treeDataDf$Group)[group], "_Gene_vs_HGT_V2",".png", sep=""), width = 550, height = 350)
  
  plot(y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])], 
       x = treeDataDf$Number_of_leaves_tgt[which(treeDataDf$Group == levels(treeDataDf$Group)[group])],
       xlab = 'Number of Genes',
       ylab = 'Fraction of Xenologs',
       main = paste(levels(treeDataDf$Group)[group], '[',
                    'D:',
                    treeDataDf$dupl_rate[group*1000-1],
                    'L:',
                    treeDataDf$loss_rate[group*1000-1],
                    'H:',
                    treeDataDf$hgt_rate[group*1000-1], ']',
                    sep = ' '),
       pch = 19,
       col = '#30303080')
  abline(mod[[1]][1], mod[[1]][2], col = 'red', lwd = 2)
  
  dev.off()  
  
  ### Save to sumDF ###
  
  sumDf$Gruppe[group] <- levels(treeDataDf$Group)[group]
  sumDf$Duplication_Rate[group] <- treeDataDf$dupl_rate[group*1000-1]
  sumDf$Loss_Rate[group] <- treeDataDf$loss_rate[group*1000-1]
  sumDf$HGT_Rate[group] <- treeDataDf$hgt_rate[group*1000-1]
  sumDf$Slope[group] <-  round(mod[[1]][2], digits = 2)
  sumDf$Intercept[group] <- round(mod[[1]][1], digits = 2)
  sumDf$Spearman_Corr[group] <- round(coef, digits = 2)
}

write.csv(sumDf, 'Results_Gene_vs_HGT.csv' , dec = '.', sep = ';')



######################
#### SPECIES VS HGT ####
######################
group = 1

hist(treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])])
hist(treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])])


for (group in 1:length(levels(treeDataDf$Group))) {
  mod = lm(treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])] ~ 
             treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])])
  
  coef <- cor(x = treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])],
              y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])], method = 'spearman') 
  
  ### Plots ###
  png(paste("02_Plots/", levels(treeDataDf$Group)[group], "_Species_vs_HGT_V2",".png", sep=""), width = 550, height = 350)
  
  plot(y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])], 
       x = treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])],
       xlab = 'Number of Species',
       ylab = 'Fraction of Xenologs',
       main = paste(levels(treeDataDf$Group)[group], '[',
                    'D:',
                    treeDataDf$dupl_rate[group*1000-1],
                    'L:',
                    treeDataDf$loss_rate[group*1000-1],
                    'H:',
                    treeDataDf$hgt_rate[group*1000-1], ']',
                    sep = ' '),
       pch = 19,
       col = '#30303080')
  abline(mod[[1]][1], mod[[1]][2], col = 'red', lwd = 2)
  
  dev.off()  
  
  ### Save to sumDF ###
  
  sumDf$Gruppe[group] <- levels(treeDataDf$Group)[group]
  sumDf$Duplication_Rate[group] <- treeDataDf$dupl_rate[group*1000-1]
  sumDf$Loss_Rate[group] <- treeDataDf$loss_rate[group*1000-1]
  sumDf$HGT_Rate[group] <- treeDataDf$hgt_rate[group*1000-1]
  sumDf$Slope[group] <-  round(mod[[1]][2], digits = 2)
  sumDf$Intercept[group] <- round(mod[[1]][1], digits = 2)
  sumDf$Spearman_Corr[group] <- round(coef, digits = 2)
  
}

write.csv(sumDf, 'Results_Species_vs_HGT.csv' , dec = '.', sep = ';')

#############
#### 2.2 ####
#############
hh <- as.numeric(levels(as.factor(treeDataDf$dupl_rate)))
#### Plots ####
box_hgt_dupl <- ggplot(treeDataDf, aes(x = factor(dupl_rate), y = Fraction_of_Xenologs, group=factor(hgt_rate))) +
  geom_jitter(shape=16, position=position_jitter(0.15), aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction of Xenologs vs. Duplication Rate', x ='Duplication Rate', y = 'HGT Events', colour = 'HGT Rate') +
  theme_bw()
box_hgt_dupl

ggsave("02_Plots/Fraction_of_Xenologs_vs._Duplication_Rate.png", box_hgt_dupl, width = 8, height = 5.1)

getwd()

box_hgt_loss <- ggplot(treeDataDf, aes(x = as.factor(loss_rate), y = Fraction_of_Xenologs, group=as.factor(hgt_rate))) +
  geom_jitter(shape=16, position=position_jitter(0.15), aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction_of_Xenologs vs. Loss Rate', x ='Loss Rate', y = 'HGT Events', colour = 'HGT Rate') +
  theme_bw()
box_hgt_loss

ggsave("02_Plots/Fraction of Xenologs_vs._Loss_Rate.png", box_hgt_loss, width = 8, height = 5.1)

#### Signifikanzen ####

# Kruskal-Wallis-Test ???ber alle Gruppen
# wenn positiv, dann Mann-Whitney-U-Test jede Gruppe gegen jede

kruskal.test(data = treeDataDf)

#############
#### 2.3 ####
#############

plot(x = treeDataDf$loss_rate[which(treeDataDf$hgt_rate == '0.5')],
     y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$hgt_rate == '0.5')])

boxplot(Fraction_of_Xenologs~factor(hgt_rate), data = treeDataDf)

plot(x = treeDataDf$Number_of_leaves_tgt,
     y = treeDataDf$Fraction_of_Xenologs)
mod = lm(treeDataDf$Fraction_of_Xenologs ~ 
           treeDataDf$Number_of_leaves_tgt)
abline(mod[[1]][1], mod[[1]][2], col = 'red', lwd = 2)


#### true negatives ####
#TN = V(V-1)/2 - (TP +FN +FP)

# true positives: alle Kanten die im true fitch und im simulierten baum
# false negatives: alle Kanten die nur im wahren baum und nicht im simulierten
# false positives: alle Kanten welche nur im simulierten und nicht im wahren braum sind
# true negatives: alle möglichen Kantenkombinationen minus die in den Bäumen erkannten

treeDataDf$tn_cd_100 <- (treeDataDf$Number_of_Nodes_ldt_100 * (treeDataDf$Number_of_Nodes_ldt_100-1) / 2)  - (treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_negative_100 + treeDataDf$Edges_cd_false_positive_100)
treeDataDf$tn_cd_80 <- (treeDataDf$Number_of_Nodes_ldt_80 * (treeDataDf$Number_of_Nodes_ldt_80-1) / 2)  - (treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_negative_80 + treeDataDf$Edges_cd_false_positive_80)
treeDataDf$tn_cd_60 <- (treeDataDf$Number_of_Nodes_ldt_60 * (treeDataDf$Number_of_Nodes_ldt_60-1) / 2)  - (treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_negative_60 + treeDataDf$Edges_cd_false_positive_60)
treeDataDf$tn_cd_40 <- (treeDataDf$Number_of_Nodes_ldt_40 * (treeDataDf$Number_of_Nodes_ldt_40-1) / 2)  - (treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_negative_40 + treeDataDf$Edges_cd_false_positive_40)
treeDataDf$tn_cd_20 <- (treeDataDf$Number_of_Nodes_ldt_20 * (treeDataDf$Number_of_Nodes_ldt_20-1) / 2)  - (treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_negative_20 + treeDataDf$Edges_cd_false_positive_20)

treeDataDf$tn_rs_100 <- (treeDataDf$Number_of_Nodes_ldt_100 * (treeDataDf$Number_of_Nodes_ldt_100-1) / 2)  - (treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_negative_100 + treeDataDf$Edges_rs_false_positive_100)
treeDataDf$tn_rs_80 <- (treeDataDf$Number_of_Nodes_ldt_80 * (treeDataDf$Number_of_Nodes_ldt_80-1) / 2)  - (treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_negative_80 + treeDataDf$Edges_rs_false_positive_80)
treeDataDf$tn_rs_60 <- (treeDataDf$Number_of_Nodes_ldt_60 * (treeDataDf$Number_of_Nodes_ldt_60-1) / 2)  - (treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_negative_60 + treeDataDf$Edges_rs_false_positive_60)
treeDataDf$tn_rs_40 <- (treeDataDf$Number_of_Nodes_ldt_40 * (treeDataDf$Number_of_Nodes_ldt_40-1) / 2)  - (treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_negative_40 + treeDataDf$Edges_rs_false_positive_40)
treeDataDf$tn_rs_20 <- (treeDataDf$Number_of_Nodes_ldt_20 * (treeDataDf$Number_of_Nodes_ldt_20-1) / 2)  - (treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_negative_20 + treeDataDf$Edges_rs_false_positive_20)

#### Recall ####
# TP /(TP + FN)

treeDataDf$recall_cd_100 <- treeDataDf$Edges_cd_true_positive_100 / (treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_negative_100)
treeDataDf$recall_cd_80 <- treeDataDf$Edges_cd_true_positive_80/(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_negative_80)
treeDataDf$recall_cd_60 <- treeDataDf$Edges_cd_true_positive_60/(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_negative_60)
treeDataDf$recall_cd_40 <- treeDataDf$Edges_cd_true_positive_40/(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_negative_40)
treeDataDf$recall_cd_20 <- treeDataDf$Edges_cd_true_positive_20/(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_negative_20)

# recall_cd <- c(recall_cd_100, recall_cd_80, recall_cd_60, recall_cd_40, recall_cd_20) ????? HÄääscht dü läck gesöffen?

treeDataDf$recall_rs_100 <- treeDataDf$Edges_rs_true_positive_100/(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_negative_100)
treeDataDf$recall_rs_80 <- treeDataDf$Edges_rs_true_positive_80/(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_negative_80)
treeDataDf$recall_rs_60 <- treeDataDf$Edges_rs_true_positive_60/(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_negative_60)
treeDataDf$recall_rs_40 <- treeDataDf$Edges_rs_true_positive_40/(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_negative_40)
treeDataDf$recall_rs_20 <- treeDataDf$Edges_rs_true_positive_20/(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_negative_20)

# recall_rs <- c(recall_rs_100,recall_rs_80,recall_rs_60,recall_rs_40,recall_rs_20)

#### Precision ####
# TP/(TP + FP)

treeDataDf$precision_cd_100 <- treeDataDf$Edges_cd_true_positive_100/(treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_positive_100)
treeDataDf$precision_cd_80 <- treeDataDf$Edges_cd_true_positive_80/(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_positive_80)
treeDataDf$precision_cd_60 <- treeDataDf$Edges_cd_true_positive_60/(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_positive_60)
treeDataDf$precision_cd_40 <- treeDataDf$Edges_cd_true_positive_40/(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_positive_40)
treeDataDf$precision_cd_20 <- treeDataDf$Edges_cd_true_positive_20/(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_positive_20)

# precision_cd <- c(precision_cd_100,precision_cd_80,precision_cd_60,precision_cd_40,precision_cd_20)

treeDataDf$precision_rs_100 <- treeDataDf$Edges_rs_true_positive_100/(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_positive_100)
treeDataDf$precision_rs_80 <- treeDataDf$Edges_rs_true_positive_80/(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_positive_80)
treeDataDf$precision_rs_60 <- treeDataDf$Edges_rs_true_positive_60/(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_positive_60)
treeDataDf$precision_rs_40 <- treeDataDf$Edges_rs_true_positive_40/(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_positive_40)
treeDataDf$precision_rs_20 <- treeDataDf$Edges_rs_true_positive_20/(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_positive_20)

# precision_rs <- c(precision_rs_100,precision_rs_80,precision_rs_60,precision_rs_40,precision_rs_20)

#### Accuracy ####
# TP + TN / (TP + TN + FP + FN)
# number_ofnodes_percent

treeDataDf$accuracy_cd_100 <- (treeDataDf$Edges_cd_true_positive_100 + treeDataDf$tn_cd_100) / (treeDataDf$Edges_cd_true_positive_100 + treeDataDf$tn_cd_100 + treeDataDf$Edges_cd_false_positive_100 + treeDataDf$Edges_cd_false_negative_100)
treeDataDf$accuracy_cd_80 <- (treeDataDf$Edges_cd_true_positive_80 + treeDataDf$tn_cd_80) /(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$tn_cd_80 + treeDataDf$Edges_cd_false_positive_80 + treeDataDf$Edges_cd_false_negative_80)
treeDataDf$accuracy_cd_60 <- (treeDataDf$Edges_cd_true_positive_60 + treeDataDf$tn_cd_60) /(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$tn_cd_60 + treeDataDf$Edges_cd_false_positive_60 + treeDataDf$Edges_cd_false_negative_60)
treeDataDf$accuracy_cd_40 <- (treeDataDf$Edges_cd_true_positive_40 + treeDataDf$tn_cd_40) /(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$tn_cd_40 + treeDataDf$Edges_cd_false_positive_40 + treeDataDf$Edges_cd_false_negative_40)
treeDataDf$accuracy_cd_20 <- (treeDataDf$Edges_cd_true_positive_20 + treeDataDf$tn_cd_20) /(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$tn_cd_20 + treeDataDf$Edges_cd_false_positive_20 + treeDataDf$Edges_cd_false_negative_20)

# accuracy_cd <- c(accuracy_cd_100,accuracy_cd_80,accuracy_cd_60,accuracy_cd_40,accuracy_cd_20)

treeDataDf$accuracy_rs_100 <- (treeDataDf$Edges_rs_true_positive_100 + treeDataDf$tn_rs_100) /(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$tn_rs_100 + treeDataDf$Edges_rs_false_positive_100 + treeDataDf$Edges_rs_false_negative_100)
treeDataDf$accuracy_rs_80 <- (treeDataDf$Edges_rs_true_positive_80 + treeDataDf$tn_rs_80) /(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$tn_rs_80 + treeDataDf$Edges_rs_false_positive_80 + treeDataDf$Edges_rs_false_negative_80)
treeDataDf$accuracy_rs_60 <- (treeDataDf$Edges_rs_true_positive_60 + treeDataDf$tn_rs_60) /(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$tn_rs_60 + treeDataDf$Edges_rs_false_positive_60 + treeDataDf$Edges_rs_false_negative_60)
treeDataDf$accuracy_rs_40 <- (treeDataDf$Edges_rs_true_positive_40 + treeDataDf$tn_rs_40) /(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$tn_rs_40 + treeDataDf$Edges_rs_false_positive_40 + treeDataDf$Edges_rs_false_negative_40)
treeDataDf$accuracy_rs_20 <- (treeDataDf$Edges_rs_true_positive_20 + treeDataDf$tn_rs_20) /(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$tn_rs_20 + treeDataDf$Edges_rs_false_positive_20 + treeDataDf$Edges_rs_false_negative_20)

# accuracy_rs <- c(accuracy_rs_100,accuracy_rs_80,accuracy_rs_60,accuracy_rs_40,accuracy_rs_20)

write.csv(treeDataDf, 'Tree_Data_Full.csv' , dec = '.', sep = ';')

#################
#### Plot nützliche Größen ####

levels(treeDataDf$Group)
# for (i in do your thing) {}

# for Cluster Deletion
sumCDrecallDf <- data.frame(Group = character(c(7)))
sumCDprecisionlDf <- data.frame(Group = character(c(7)))
sumCDaccuracyDf <- data.frame(Group = character(c(7)))

# for Fitch 
sumRSrecallDf <- data.frame(Group = character(c(7)))
sumRSprecisionlDf <- data.frame(Group = character(c(7)))
sumRSaccuracyDf <- data.frame(Group = character(c(7)))
i = 1
for (i in 1:length(levels(treeDataDf$Group))) {
  sumRSrecallDf$Group[i] <- levels(treeDataDf$Group)[i]
  sumRSprecisionlDf$Group[i] <- levels(treeDataDf$Group)[i]
  sumRSaccuracyDf$Group[i] <- levels(treeDataDf$Group)[i]
 
  sumCDrecallDf$Group[i] <- levels(treeDataDf$Group)[i]
  sumCDprecisionlDf$Group[i] <- levels(treeDataDf$Group)[i]
  sumCDaccuracyDf$Group[i] <- levels(treeDataDf$Group)[i]
  # recall cd
  sumCDrecallDf$recall_cd_mean_100[i] <- mean(treeDataDf$recall_cd_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDrecallDf$recall_cd_mean_80[i] <- mean(treeDataDf$recall_cd_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDrecallDf$recall_cd_mean_60[i] <- mean(treeDataDf$recall_cd_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDrecallDf$recall_cd_mean_40[i] <- mean(treeDataDf$recall_cd_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDrecallDf$recall_cd_mean_20[i] <- mean(treeDataDf$recall_cd_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  #precision cd
  sumCDprecisionlDf$precision_cd_mean_100[i] <- mean(treeDataDf$precision_cd_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDprecisionlDf$precision_cd_mean_80[i] <- mean(treeDataDf$precision_cd_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDprecisionlDf$precision_cd_mean_60[i] <- mean(treeDataDf$precision_cd_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDprecisionlDf$precision_cd_mean_40[i] <- mean(treeDataDf$precision_cd_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDprecisionlDf$precision_cd_mean_20[i] <- mean(treeDataDf$precision_cd_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  # accuracy cd
  sumCDaccuracyDf$accuracy_cd_mean_100[i] <- mean(treeDataDf$accuracy_cd_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDaccuracyDf$accuracy_cd_mean_80[i] <- mean(treeDataDf$accuracy_cd_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDaccuracyDf$accuracy_cd_mean_60[i] <- mean(treeDataDf$accuracy_cd_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDaccuracyDf$accuracy_cd_mean_40[i] <- mean(treeDataDf$accuracy_cd_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumCDaccuracyDf$accuracy_cd_mean_20[i] <- mean(treeDataDf$accuracy_cd_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  
  #######################
  # recall rs
  sumRSrecallDf$recall_rs_mean_100[i] <- mean(treeDataDf$recall_rs_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSrecallDf$recall_rs_mean_80[i] <- mean(treeDataDf$recall_rs_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSrecallDf$recall_rs_mean_60[i] <- mean(treeDataDf$recall_rs_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSrecallDf$recall_rs_mean_40[i] <- mean(treeDataDf$recall_rs_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSrecallDf$recall_rs_mean_20[i] <- mean(treeDataDf$recall_rs_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  #precision rs
  sumRSprecisionlDf$precision_rs_mean_100[i] <- mean(treeDataDf$precision_rs_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSprecisionlDf$precision_rs_mean_80[i] <- mean(treeDataDf$precision_rs_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSprecisionlDf$precision_rs_mean_60[i] <- mean(treeDataDf$precision_rs_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSprecisionlDf$precision_rs_mean_40[i] <- mean(treeDataDf$precision_rs_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSprecisionlDf$precision_rs_mean_20[i] <- mean(treeDataDf$precision_rs_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  # accuracy rs
  #i=2
  
  sumRSaccuracyDf$accuracy_rs_mean_100[i] <- mean(treeDataDf$accuracy_rs_100[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSaccuracyDf$accuracy_rs_mean_80[i] <- mean(treeDataDf$accuracy_rs_80[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSaccuracyDf$accuracy_rs_mean_60[i] <- mean(treeDataDf$accuracy_rs_60[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSaccuracyDf$accuracy_rs_mean_40[i] <- mean(treeDataDf$accuracy_rs_40[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
  sumRSaccuracyDf$accuracy_rs_mean_20[i] <- mean(treeDataDf$accuracy_rs_20[which(treeDataDf$Group == levels(treeDataDf$Group)[i])], na.rm = T)
}

#PsumCDrecallDf <- pivot_longer(sumCDrecallDf, cols = -c("Group"), names_to = "Measure")

#### Correct stuff ####
cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumCDrecallDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}

PsumCDprecisionlDf <- pivot_longer(sumCDprecisionlDf, cols = -c("Group"), names_to = "Measure")
cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumCDprecisionlDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}
PsumCDaccuracyDf <- pivot_longer(sumCDaccuracyDf, cols = -c("Group"), names_to = "Measure")
cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumCDaccuracyDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}

PsumRSrecallDf <- pivot_longer(sumRSrecallDf, cols = -c("Group"), names_to = "Measure")
cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumRSrecallDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}
PsumRSprecisionlDf <- pivot_longer(sumRSprecisionlDf, cols = -c("Group"), names_to = "Measure")

cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumRSprecisionlDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}
PsumRSaccuracyDf <- pivot_longer(sumRSaccuracyDf, cols = -c("Group"), names_to = "Measure")

cor <- 0
for (i in 1:35) {
  cor2 <- 0
  for (j in 1:5) {
    PsumRSaccuracyDf$Measure[j+cor] <- 100 - cor2
    cor2 <- cor2 + 20
  }
  cor <- cor + 5
}

###############################
####  CD plot recall ####

recall_CD_Plot <- ggplot(PsumCDrecallDf, aes(x = as.numeric(Measure),
                                        y = value, 
                                        group = as.factor(Group))) +
  labs(title = 'Cluster Deletion: Mean Recall of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Recall', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
recall_CD_Plot

ggsave("02_Plots/Recall_CD_Groups.png", recall_CD_Plot, width = 8, height = 5.1)
# Je kleiner die Bäume sind desto geringer wird der Anteil der identischischen/richtigen Kanten (Xenologe)
# je weniger kannten desto mehr fallen die Fehlenden kanten ins Gewicht

###############################
#### CD plot accuracy  ####

acc_CD_Plot <- ggplot(PsumCDaccuracyDf, aes(x = as.numeric(Measure),
                                             y = value, 
                                             group = as.factor(Group))) +
  labs(title = 'Cluster Deletion: Mean Accuracy of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Accuracy', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
acc_CD_Plot

ggsave("02_Plots/Accuracy_CD_Groups.png", acc_CD_Plot, width = 8, height = 5.1)
# Jeh mehr Knoten, desto mehr 'Fehler' werden gemacht.
# Wie 

###############################
#### CD plot  Preciscion  ####

prec_CD_Plot <- ggplot(PsumCDprecisionlDf, aes(x = as.numeric(Measure),
                                            y = value, 
                                            group = as.factor(Group))) +
  labs(title = 'Cluster Deletion: Mean Precision of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Precision', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
prec_CD_Plot

ggsave("02_Plots/Prec_CD_Groups.png", acc_CD_Plot, width = 8, height = 5.1)

###############################
###############################

###############################
####  RS plot recall ####

recall_RS_Plot <- ggplot(PsumRSrecallDf, aes(x = as.numeric(Measure),
                                             y = value, 
                                             group = as.factor(Group))) +
  labs(title = 'RS Fitch: Mean Recall of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Recall', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
recall_RS_Plot

ggsave("02_Plots/Recall_RS_Groups.png", recall_RS_Plot, width = 8, height = 5.1)
# Je kleiner die Bäume sind desto geringer wird der Anteil der identischischen/richtigen Kanten (Xenologe)
# je weniger kannten desto mehr fallen die Fehlenden kanten ins Gewicht

###############################
#### RS plot accuracy  ####

acc_RS_Plot <- ggplot(PsumRSaccuracyDf, aes(x = as.numeric(Measure),
                                            y = value, 
                                            group = as.factor(Group))) +
  labs(title = 'RS Fitch: Mean Accuracy of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Accuracy', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
acc_RS_Plot

ggsave("02_Plots/Accuracy_RS_Groups.png", acc_RS_Plot, width = 8, height = 5.1)
# Jeh mehr Knoten, desto mehr 'Fehler' werden gemacht.
# Wie 

###############################
#### RS plot  Preciscion  ####

prec_RS_Plot <- ggplot(PsumRSprecisionlDf, aes(x = as.numeric(Measure),
                                               y = value, 
                                               group = as.factor(Group))) +
  labs(title = 'RS Fitch: Mean Precision of Groups', 
       x ='Percentage of Original Graph', 
       y = 'Precision', 
       colour = 'Group') +
  geom_line(aes(color = factor(Group))) +
  geom_point(aes(color = factor(Group))) +
  theme_bw()
prec_RS_Plot

ggsave("02_Plots/Prec_RS_Groups.png", acc_RS_Plot, width = 8, height = 5.1)


