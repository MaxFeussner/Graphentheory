rm(list=ls())                                  # Clear Workspace to avoid data mistakes

list.of.packages <- c("rstudioapi",            # import find path script
                      "ggplot2",               # import plotting library
                      "ggpubr",                # Align plots in grid
                      "RColorBrewer",          # for beautiful colored 
                      "dplyr",                  # rename columns easely and more
                      "plotly",
                      "tidyverse",
                      "kableExtra",
                      "gridExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(plotly)
library(tidyverse)
library(kableExtra)
library(gridExtra)

pathScript <- dirname(rstudioapi::getSourceEditorContext()$path) # den Pfad des Scriptes finden
setwd(pathScript)  # set Working directory

dataDir <- paste(pathScript, "01_Data", sep="/")
treeDataDf <- read.csv(paste(pathScript, 'Tree_data.csv', sep = '/'), sep = ',', dec = ".", header = TRUE, stringsAsFactors = F)

#i = 1
for (i in 1:length(treeDataDf$ID)) {
  treeDataDf$Group[i] <- strsplit(as.character(treeDataDf$ID[i]), "_")[[1]][1]
}
treeDataDf$Group <- as.factor(treeDataDf$Group)

####  2.1 ####

#Is there a dependence on the size of the gene tree, i.e., the number of species and genes? ?

#### Preparation for SUMDF ####
# create an empty DF
sumDfGenes <- data.frame(Group = as.character(),
                    Duplication_Rate = as.numeric(),
                    Loss_Rate = as.numeric(),
                    HGT_Rate = as.numeric(),
                    Slope = as.numeric(),
                    Intercept = as.numeric(),
                    Spearman_Corr = as.numeric())
for (i in 1:length(levels(treeDataDf$Group))) {
  sumDfGenes[i,1] <- 0  
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
  plot <- ggplot(subset(treeDataDf, Group == levels(treeDataDf$Group)[group]),aes(x = Number_of_leaves_tgt,
                                                       y = Fraction_of_Xenologs)) +
    geom_point(color = '#30303080') +
    labs(title = paste(levels(treeDataDf$Group)[group], 'Fraction of Xenelogs:', '[',
                       'D:',
                       treeDataDf$dupl_rate[group*1000-1],
                       'L:',
                       treeDataDf$loss_rate[group*1000-1],
                       'H:',
                       treeDataDf$hgt_rate[group*1000-1], ']',
                       sep = ' '), 
         x ='Number of Genes', 
         y = 'Fraction of Xenologs') +
    ylim(0,1) +
    geom_smooth(method='lm', colour = "red", size = 0.5) +
    theme_bw()
  plot
  
  ggsave(paste("02_Plots/", levels(treeDataDf$Group)[group], "_Gene_vs_HGT",".png", sep=""), plot, width = 8, height = 5.1)
  
  ### Save to sumDF ###
  
  sumDfGenes$Group[group] <- levels(treeDataDf$Group)[group]
  sumDfGenes$Duplication_Rate[group] <- treeDataDf$dupl_rate[group*1000-1]
  sumDfGenes$Loss_Rate[group] <- treeDataDf$loss_rate[group*1000-1]
  sumDfGenes$HGT_Rate[group] <- treeDataDf$hgt_rate[group*1000-1]
  sumDfGenes$Slope[group] <-  round(mod[[1]][2], digits = 4)
  sumDfGenes$Intercept[group] <- round(mod[[1]][1], digits = 2)
  sumDfGenes$Spearman_Corr[group] <- round(coef, digits = 2)
}

write.csv(sumDfGenes, 'Results_Gene_vs_HGT.csv' , dec = '.', sep = ';')

######################
#### SPECIES VS HGT ####
######################
group = 1

sumDfSpecies <- data.frame(Group = as.character(),
                    Duplication_Rate = as.numeric(),
                    Loss_Rate = as.numeric(),
                    HGT_Rate = as.numeric(),
                    Slope = as.numeric(),
                    Intercept = as.numeric(),
                    Spearman_Corr = as.numeric())
for (i in 1:length(levels(treeDataDf$Group))) {
  sumDfSpecies[i,1] <- 0  
}

for (group in 1:length(levels(treeDataDf$Group))) {
  mod = lm(treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])] ~ 
             treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])])
  
  coef <- cor(x = treeDataDf$Number_of_Species[which(treeDataDf$Group == levels(treeDataDf$Group)[group])],
              y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$Group == levels(treeDataDf$Group)[group])], method = 'spearman') 
  
  ### Plots ###
  plot <- ggplot(subset(treeDataDf, Group == levels(treeDataDf$Group)[group]),aes(x = Number_of_Species,
                                                                                  y = Fraction_of_Xenologs)) +
    geom_point(color = '#30303080') +
    labs(title = paste(levels(treeDataDf$Group)[group], 'Fraction of Xenelogs:', '[',
                       'D:',
                       treeDataDf$dupl_rate[group*1000-1],
                       'L:',
                       treeDataDf$loss_rate[group*1000-1],
                       'H:',
                       treeDataDf$hgt_rate[group*1000-1], ']',
                       sep = ' '), 
         x ='Number of Species', 
         y = 'Fraction of Xenologs') +
    ylim(0,1) +
    geom_smooth(method='lm', colour = "red", size = 0.5) +
    theme_bw()
  plot
  
  ggsave(paste("02_Plots/", levels(treeDataDf$Group)[group], "_Species_vs_HGT",".png", sep=""), plot, width = 8, height = 5.1)
  
  ### Save to sumDF ###
  
  sumDfSpecies$Group[group] <- levels(treeDataDf$Group)[group]
  sumDfSpecies$Duplication_Rate[group] <- treeDataDf$dupl_rate[group*1000-1]
  sumDfSpecies$Loss_Rate[group] <- treeDataDf$loss_rate[group*1000-1]
  sumDfSpecies$HGT_Rate[group] <- treeDataDf$hgt_rate[group*1000-1]
  sumDfSpecies$Slope[group] <-  round(mod[[1]][2], digits = 4)
  sumDfSpecies$Intercept[group] <- round(mod[[1]][1], digits = 2)
  sumDfSpecies$Spearman_Corr[group] <- round(coef, digits = 2)
  
}

write.csv(sumDfSpecies, 'Results_Species_vs_HGT.csv' , dec = '.', sep = ';')

#############
#### 2.2  Fixed HGT Rate ####
#############
# How does the fraction depend on the rate of duplications and losses for a fixed horizontal transfer rate?
#### Plots ####
box_hgt_dupl <- ggplot(treeDataDf, aes(x = factor(dupl_rate), 
                                       y = Fraction_of_Xenologs, 
                                       group=factor(hgt_rate))) +
  geom_jitter(shape=16, 
              position=position_jitter(0.15), 
              aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, 
               show.legend = FALSE, 
               aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction of Xenologs vs. Duplication Rate', 
       x = 'Duplication Rate', 
       y = 'Fraction of Xenelogs', 
       colour = 'HGT Rate') +
  theme_bw()
box_hgt_dupl

ggsave("02_Plots/Fraction_of_Xenologs_vs._Duplication_Rate.png", box_hgt_dupl, width = 8, height = 5.1)

box_hgt_loss <- ggplot(treeDataDf, aes(x = as.factor(loss_rate), 
                                       y = Fraction_of_Xenologs, 
                                       group=as.factor(hgt_rate))) +
  geom_jitter(shape=16, 
              position=position_jitter(0.15), 
              aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, 
               show.legend = FALSE, 
               aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction_of_Xenologs vs. Loss Rate', 
       x = 'Loss Rate', 
       y = 'Fraction of Xenelogs', 
       colour = 'HGT Rate') +
  theme_bw()
box_hgt_loss

ggsave("02_Plots/Fraction of Xenologs_vs._Loss_Rate.png", box_hgt_loss, width = 8, height = 5.1)

#### Signifikanzen ####

# Kruskal-Wallis-Test ???ber alle Gruppen
# wenn positiv, dann Mann-Whitney-U-Test jede Gruppe gegen jede
#kruskal.test(data = treeDataDf)

#############
#### 2.3 ####
#############

#### To Do ####
# How does the fraction depend on the horizontal transfer rate on the the rate of duplications and losses? original


# How does the fraction depend on the horizontal transfer rate with a fixed duplication and loss rate? --> fixed??

treeDataDf$d

test <- ggplot(treeDataDf, aes(x = as.factor(hgt_rate), 
                               y = Fraction_of_Xenologs, 
                               group=as.factor(dupl_rate))) +
  geom_jitter(shape=16, 
              position=position_jitter(0.15), 
              aes(color = factor(hgt_rate), group=factor(dupl_rate))) +
  geom_boxplot(outlier.shape = NA, 
               show.legend = FALSE, 
               aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction of Xenologs vs. HGT Rate', 
       x = 'HGT Rate', 
       y = 'Fraction of Xenelogs', 
       colour = 'Duplikation Rate') +
  theme_bw()
test

ggsave("02_Plots/09_Dupl_Fraction_of_Xenologs_vs_HGT_Rate.png", test, width = 8, height = 5.1)

test1 <- ggplot(treeDataDf, aes(x = as.factor(hgt_rate), 
                               y = Fraction_of_Xenologs, 
                               group=as.factor(loss_rate))) +
  geom_jitter(shape=16, 
              position=position_jitter(0.15), 
              aes(color = factor(hgt_rate), group=factor(loss_rate))) +
  geom_boxplot(outlier.shape = NA, 
               show.legend = FALSE, 
               aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction of Xenologs vs. HGT Rate', 
       x = 'HGT Rate', 
       y = 'Fraction of Xenelogs', 
       colour = 'Loss Rate') +
  theme_bw()
test1

ggsave("02_Plots/10_Loss_Fraction_of_Xenologs_vs._Loss_Rate.png", test1, width = 8, height = 5.1)

plot(x = treeDataDf$loss_rate[which(treeDataDf$hgt_rate == '0.5')],
     y = treeDataDf$Fraction_of_Xenologs[which(treeDataDf$hgt_rate == '0.5')])

boxplot(Fraction_of_Xenologs~factor(hgt_rate), data = treeDataDf)

plot(x = treeDataDf$Number_of_leaves_tgt,
     y = treeDataDf$Fraction_of_Xenologs)
mod = lm(treeDataDf$Fraction_of_Xenologs ~ 
           treeDataDf$Number_of_leaves_tgt)
abline(mod[[1]][1], mod[[1]][2], col = 'red', lwd = 2)

plot(x = treeDataDf$non_binary_prob,
     y = treeDataDf$Fraction_of_Xenologs)


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

## Save human readable tables ##

write.csv(sumCDrecallDf, 'Summary_CD_Recall.csv' , dec = '.', sep = ';')
write.csv(sumCDprecisionlDf, 'Summary_CD_Precision.csv' , dec = '.', sep = ';')
write.csv(sumCDaccuracyDf, 'Summary_CD_Accuracy.csv' , dec = '.', sep = ';')

write.csv(sumRSrecallDf, 'Summary_Fitch_Recall.csv' , dec = '.', sep = ';')
write.csv(sumRSprecisionlDf, 'Summary_Fitch_Precision.csv' , dec = '.', sep = ';')
write.csv(sumRSaccuracyDf, 'Summary_Fitch_Accuracy.csv' , dec = '.', sep = ';')

PsumCDrecallDf <- pivot_longer(sumCDrecallDf, cols = -c("Group"), names_to = "Measure")

#### Correct stuff ####
PsumCDrecallDf <- pivot_longer(sumCDrecallDf, cols = -c("Group"), names_to = "Measure")
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




#### Tripple ####
names(treeDataDf)
#### True Negatives ####

treeDataDf$T_LDT_True_Neg <- (treeDataDf$Number_of_Nodes_ldt_100 * (treeDataDf$Number_of_Nodes_ldt_100-1) / 2)  - (treeDataDf$T_ldt_true_positive + treeDataDf$T_ldt_false_negative + treeDataDf$T_ldt_false_positive)
treeDataDf$S_LDT_True_Neg <- (treeDataDf$Number_of_Nodes_ldt_100 * (treeDataDf$Number_of_Nodes_ldt_100-1) / 2)  - (treeDataDf$S_ldt_true_positive + treeDataDf$T_ldt_false_negative + treeDataDf$S_ldt_false_positive)

# recall 
treeDataDf$T_LDT_recall <- treeDataDf$T_ldt_true_positive / (treeDataDf$T_ldt_true_positive + treeDataDf$T_ldt_false_negative)
treeDataDf$S_LDT_recall <- treeDataDf$S_ldt_true_positive / (treeDataDf$S_ldt_true_positive + treeDataDf$S_ldt_false_negative)

# precision
treeDataDf$T_LDT_precision <- treeDataDf$T_ldt_true_positive/(treeDataDf$T_ldt_true_positive + treeDataDf$T_ldt_false_positive)
treeDataDf$S_LDT_precision <- treeDataDf$S_ldt_true_positive/(treeDataDf$S_ldt_true_positive + treeDataDf$S_ldt_false_positive)

# accuracy
treeDataDf$T_LDT_accuracy <- (treeDataDf$T_ldt_true_positive + treeDataDf$T_LDT_True_Neg) / (treeDataDf$T_ldt_true_positive + treeDataDf$T_LDT_True_Neg + treeDataDf$T_ldt_false_positive + treeDataDf$T_ldt_false_negative)
treeDataDf$S_LDT_accuracy <- (treeDataDf$S_ldt_true_positive + treeDataDf$S_LDT_True_Neg) / (treeDataDf$S_ldt_true_positive + treeDataDf$S_LDT_True_Neg + treeDataDf$S_ldt_false_positive + treeDataDf$S_ldt_false_negative)

#### Tripple T Plots ####
#### recall ####

gruppen <- levels(treeDataDf$Group)
sumTripple <- data.frame(Group = character(length(gruppen)))
for (i in 1:length(gruppen)) {
  sumTripple$Group[i] <- gruppen[i]
  sumTripple$T_LDT_Recall_Mean[i] <- round(mean(treeDataDf$T_LDT_recall[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTripple$T_LDT_Recall_Median[i] <- round(median(treeDataDf$T_LDT_recall[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTripple$S_LDT_Recall_Mean[i] <- round(mean(treeDataDf$S_LDT_recall[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTripple$S_LDT_Recall_Median[i] <- round(median(treeDataDf$S_LDT_recall[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
}

recall_T_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                       y = T_LDT_recall, 
                                       group = as.factor(Group))) +
  labs(title = 'Triple T: Mean Recall of Groups', 
       x ='Group', 
       y = 'Recall', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0,0.25) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
recall_T_LDT

ggsave("02_Plots/Tripple_T_Recall_Groups.png", recall_T_LDT, width = 8, height = 5.1)

#### precision ###

sumPrecision <- data.frame(Group = character(length(gruppen)))
for (i in 1:length(gruppen)) {
  sumPrecision$Group[i] <- gruppen[i]
  sumPrecision$T_LDT_Precision_Mean[i] <- round(mean(treeDataDf$T_LDT_precision[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumPrecision$T_LDT_Precision_Median[i] <- round(median(treeDataDf$T_LDT_precision[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumPrecision$S_LDT_Precision_Mean[i] <- round(mean(treeDataDf$S_LDT_precision[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumPrecision$S_LDT_Precision_Median[i] <- round(median(treeDataDf$S_LDT_precision[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
}

precision_T_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                       y = T_LDT_precision, 
                                       group = as.factor(Group))) +
  labs(title = 'Triple Gene Tree: Mean Precision of Groups', 
       x ='Group', 
       y = 'Precision', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  #ylim(0,0.25) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
precision_T_LDT

ggsave("02_Plots/Tripple_T_Precision_Groups.png", precision_T_LDT, width = 8, height = 5.1)

#### accuracy ###

sumAccuracy <- data.frame(Group = character(length(gruppen)))
for (i in 1:length(gruppen)) {
  sumAccuracy$Group[i] <- gruppen[i]
  sumAccuracy$T_LDT_Accuracy_Mean[i] <- round(mean(treeDataDf$T_LDT_accuracy[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumAccuracy$T_LDT_Accuracy_Median[i] <- round(median(treeDataDf$T_LDT_accuracy[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumAccuracy$S_LDT_Accuracy_Mean[i] <- round(mean(treeDataDf$S_LDT_accuracy[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumAccuracy$S_LDT_Accuracy_Median[i] <- round(median(treeDataDf$S_LDT_accuracy[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
}

accuracy_T_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                          y = T_LDT_accuracy, 
                                          group = as.factor(Group))) +
  labs(title = 'Triple T: Mean Precision of Groups', 
       x ='Group', 
       y = 'Accuracy', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0,0.002) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
accuracy_T_LDT

ggsave("02_Plots/Tripple_T_Accuracy_Groups.png", accuracy_T_LDT, width = 8, height = 5.1)

#### Tripple S Plots ####
#### recall ####
recall_S_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                       y = S_LDT_recall, 
                                       group = as.factor(Group))) +
  labs(title = 'Triple S: Mean Recall of Groups', 
       x ='Group', 
       y = 'Recall', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0.99,1.1) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
recall_S_LDT

ggsave("02_Plots/Tripple_S_Recall_Groups.png", recall_S_LDT, width = 8, height = 5.1)

#### precision ###
precision_S_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                          y = S_LDT_precision, 
                                          group = as.factor(Group))) +
  labs(title = 'Triple S: Mean Precision of Groups', 
       x ='Group', 
       y = 'Precision', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0,0.01) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
precision_S_LDT

ggsave("02_Plots/Tripple_S_Precision_Groups.png", precision_S_LDT, width = 8, height = 5.1)

#### accuracy ###
accuracy_S_LDT <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                         y = S_LDT_accuracy, 
                                         group = as.factor(Group))) +
  labs(title = 'Triple S: Mean Accuracy of Groups', 
       x ='Group', 
       y = 'Accuracy', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(-0.1,1) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=3) +
  theme_bw()
accuracy_S_LDT

ggsave("02_Plots/Tripple_S_Accuracy_Groups.png", accuracy_S_LDT, width = 8, height = 5.1)

#### Tripple Besonders ####

names(treeDataDf)

#TLDTripple / t-tripple


treeDataDf$T_Triple_Fraction <- treeDataDf$T_ldt_triple / treeDataDf$T_triple
treeDataDf$S_Triple_Fraction <- treeDataDf$S_ldt_triple / treeDataDf$S_triple

sumTriple_Fraction <- data.frame(Group = character(length(gruppen)))
for (i in 1:length(gruppen)) {
  sumTriple_Fraction$Group[i] <- gruppen[i]
  sumTriple_Fraction$T_Triple_Mean[i] <- round(mean(treeDataDf$T_Triple_Fraction[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTriple_Fraction$T_Triple_Median[i] <- round(median(treeDataDf$T_Triple_Fraction[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTriple_Fraction$S_Triple_Mean[i] <- round(mean(treeDataDf$S_Triple_Fraction[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
  sumTriple_Fraction$S_Triple_Median[i] <- round(median(treeDataDf$S_Triple_Fraction[which(treeDataDf$Group == gruppen[i])], na.rm = T), digits = 4)
}
sumTriple_Fraction

T_Triple_Fraction <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                         y = T_Triple_Fraction, 
                                         group = as.factor(Group))) +
  labs(title = 'Triple T: Fraction of Triples within the LDT Graph', 
       x ='Group', 
       y = 'Fraction of Triples in LDT-Graph', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0,0.70) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  theme_bw()
T_Triple_Fraction

ggsave("02_Plots/Tripple_T_Fractions_of_T_Tripples_Groups.png", T_Triple_Fraction, width = 8, height = 5.1)


S_Triple_Fraction <- ggplot(treeDataDf, aes(x = as.factor(Group),
                                            y = S_Triple_Fraction, 
                                            group = as.factor(Group))) +
  labs(title = 'Triple S: Fraction of Triples within the LDT Graph', 
       x ='Group', 
       y = 'Fraction of Triples in LDT-Graph', 
       colour = 'Group') +
  geom_boxplot(aes(color = factor(Group)), outlier.shape = NA, show.legend = FALSE) +
  ylim(0,0.70) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  theme_bw()
S_Triple_Fraction

ggsave("02_Plots/Tripple_S_Fractions_of_S_Tripples_Groups.png", S_Triple_Fraction, width = 8, height = 5.1)
