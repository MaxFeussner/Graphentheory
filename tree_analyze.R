
rm(list=ls())                                  # Clear Workspace to avoid data mistakes

list.of.packages <- c("rstudioapi",            # import find path script
                      "ggplot2",               # import plotting library
                      "ggpubr",                # Align plots in grid
                      "RColorBrewer",          # for beautiful colored 
                      "dplyr",                  # rename columns easely and more
                      "plotly"
                      
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(plotly)
pathScript <- dirname(rstudioapi::getSourceEditorContext()$path) # den Pfad des Scriptes finden
setwd(pathScript)  # set Working directory

dataDir <- paste(pathScript, "01_Data", sep="/")

treeDataDf <- read.table(paste(pathScript, 'Tree_data.csv', sep = '/'), sep = ',', header = TRUE)
# treeDataDf$Group <- as.factor(strsplit(treeDataDf$ID, '_')[[1]][1])

i = 1
for (i in 1:length(treeDataDf$ID)) {
  treeDataDf$Group[i] <- strsplit(as.character(treeDataDf$ID[i]), "_")[[1]][1]
}
treeDataDf$Group <- as.factor(treeDataDf$Group)

names(treeDataDf)

levels(treeDataDf$Group)

hist(treeDataDf$Fraction_of_Xenologs)

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

write.csv()


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

ggsave("Fraction_of_Xenologs_vs._Duplication_Rate.png", box_hgt_dupl)


box_hgt_loss <- ggplot(treeDataDf, aes(x = as.factor(loss_rate), y = Fraction_of_Xenologs, group=as.factor(hgt_rate))) +
  geom_jitter(shape=16, position=position_jitter(0.15), aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction_of_Xenologs vs. Loss Rate', x ='Loss Rate', y = 'HGT Events', colour = 'HGT Rate') +
  theme_bw()
box_hgt_loss

ggsave("Fraction of Xenologs_vs._Loss_Rate.png", box_hgt_loss)

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

tn_cd_100 <- V - (treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_negative_100 + treeDataDf$Edges_cd_false_positive_100)
tn_cd_80 <- V - (treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_negative_80 + treeDataDf$Edges_cd_false_positive_80)
tn_cd_60 <- V - (treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_negative_60 + treeDataDf$Edges_cd_false_positive_60)
tn_cd_40 <- V - (treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_negative_40 + treeDataDf$Edges_cd_false_positive_40)
tn_cd_20 <- V - (treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_negative_20 + treeDataDf$Edges_cd_false_positive_20)

tn_rs_100 <- V - (treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_negative_100 + treeDataDf$Edges_rs_false_positive_100)
tn_rs_80 <- V - (treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_negative_80 + treeDataDf$Edges_rs_false_positive_80)
tn_rs_60 <- V - (treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_negative_60 + treeDataDf$Edges_rs_false_positive_60)
tn_rs_40 <- V - (treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_negative_40 + treeDataDf$Edges_rs_false_positive_40)
tn_rs_20 <- V - (treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_negative_20 + treeDataDf$Edges_rs_false_positive_20)

#### Recall ####
# TP /(TP + FN)

recall_cd_100 <- mean(treeDataDf$Edges_cd_true_positive_100/(treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_negative_100), na.rm =TRUE)
recall_cd_80 <- mean(treeDataDf$Edges_cd_true_positive_80/(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_negative_80), na.rm =TRUE)
recall_cd_60 <- mean(treeDataDf$Edges_cd_true_positive_60/(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_negative_60), na.rm =TRUE)
recall_cd_40 <- mean(treeDataDf$Edges_cd_true_positive_40/(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_negative_40), na.rm =TRUE)
recall_cd_20 <- mean(treeDataDf$Edges_cd_true_positive_20/(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_negative_20), na.rm =TRUE)

recall_cd <- c(recall_cd_100,recall_cd_80,recall_cd_60,recall_cd_40,recall_cd_20)

recall_rs_100 <- mean(treeDataDf$Edges_rs_true_positive_100/(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_negative_100), na.rm =TRUE)
recall_rs_80 <- mean(treeDataDf$Edges_rs_true_positive_80/(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_negative_80), na.rm =TRUE)
recall_rs_60 <- mean(treeDataDf$Edges_rs_true_positive_60/(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_negative_60), na.rm =TRUE)
recall_rs_40 <- mean(treeDataDf$Edges_rs_true_positive_40/(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_negative_40), na.rm =TRUE)
recall_rs_20 <- mean(treeDataDf$Edges_rs_true_positive_20/(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_negative_20), na.rm =TRUE)

recall_rs <- c(recall_rs_100,recall_rs_80,recall_rs_60,recall_rs_40,recall_rs_20)

#### Precision ####
# TP/(TP + FP)

precision_cd_100 <- mean(treeDataDf$Edges_cd_true_positive_100/(treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_false_positive_100), na.rm =TRUE)
precision_cd_80 <- mean(treeDataDf$Edges_cd_true_positive_80/(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_false_positive_80), na.rm =TRUE)
precision_cd_60 <- mean(treeDataDf$Edges_cd_true_positive_60/(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_false_positive_60), na.rm =TRUE)
precision_cd_40 <- mean(treeDataDf$Edges_cd_true_positive_40/(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_false_positive_40), na.rm =TRUE)
precision_cd_20 <- mean(treeDataDf$Edges_cd_true_positive_20/(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_false_positive_20), na.rm =TRUE)

precision_cd <- c(precision_cd_100,precision_cd_80,precision_cd_60,precision_cd_40,precision_cd_20)

precision_rs_100 <- mean(treeDataDf$Edges_rs_true_positive_100/(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_false_positive_100), na.rm =TRUE)
precision_rs_80 <- mean(treeDataDf$Edges_rs_true_positive_80/(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_false_positive_80), na.rm =TRUE)
precision_rs_60 <- mean(treeDataDf$Edges_rs_true_positive_60/(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_false_positive_60), na.rm =TRUE)
precision_rs_40 <- mean(treeDataDf$Edges_rs_true_positive_40/(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_false_positive_40), na.rm =TRUE)
precision_rs_20 <- mean(treeDataDf$Edges_rs_true_positive_20/(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_false_positive_20), na.rm =TRUE)

precision_rs <- c(precision_rs_100,precision_rs_80,precision_rs_60,precision_rs_40,precision_rs_20)

#### Accuracy ####
# TP + TN / (TP + TN + FP + FN)
# number_ofnodes_percent

accuracy_cd_100 <- mean((treeDataDf$Edges_cd_true_positive_100 + treeDataDf$"truenegatives") /(treeDataDf$Edges_cd_true_positive_100 + treeDataDf$Edges_cd_true_negative_100 + treeDataDf$Edges_cd_false_positive_100 + treeDataDf$Edges_cd_false_negative_100), na.rm = TRUE)
accuracy_cd_80 <- mean((treeDataDf$Edges_cd_true_positive_80 + treeDataDf$"truenegatives") /(treeDataDf$Edges_cd_true_positive_80 + treeDataDf$Edges_cd_true_negative_80 + treeDataDf$Edges_cd_false_positive_80 + treeDataDf$Edges_cd_false_negative_80), na.rm = TRUE)
accuracy_cd_60 <- mean((treeDataDf$Edges_cd_true_positive_60 + treeDataDf$"truenegatives") /(treeDataDf$Edges_cd_true_positive_60 + treeDataDf$Edges_cd_true_negative_60 + treeDataDf$Edges_cd_false_positive_60 + treeDataDf$Edges_cd_false_negative_60), na.rm = TRUE)
accuracy_cd_40 <- mean((treeDataDf$Edges_cd_true_positive_40 + treeDataDf$"truenegatives") /(treeDataDf$Edges_cd_true_positive_40 + treeDataDf$Edges_cd_true_negative_40 + treeDataDf$Edges_cd_false_positive_40 + treeDataDf$Edges_cd_false_negative_40), na.rm = TRUE)
accuracy_cd_20 <- mean((treeDataDf$Edges_cd_true_positive_20 + treeDataDf$"truenegatives") /(treeDataDf$Edges_cd_true_positive_20 + treeDataDf$Edges_cd_true_negative_20 + treeDataDf$Edges_cd_false_positive_20 + treeDataDf$Edges_cd_false_negative_20), na.rm = TRUE)

accuracy_cd <- c(accuracy_cd_100,accuracy_cd_80,accuracy_cd_60,accuracy_cd_40,accuracy_cd_20)

accuracy_rs_100 <- mean((treeDataDf$Edges_rs_true_positive_100 + treeDataDf$"truenegatives") /(treeDataDf$Edges_rs_true_positive_100 + treeDataDf$Edges_rs_true_negative_100 + treeDataDf$Edges_rs_false_positive_100 + treeDataDf$Edges_rs_false_negative_100), na.rm = TRUE)
accuracy_rs_80 <- mean((treeDataDf$Edges_rs_true_positive_80 + treeDataDf$"truenegatives") /(treeDataDf$Edges_rs_true_positive_80 + treeDataDf$Edges_rs_true_negative_80 + treeDataDf$Edges_rs_false_positive_80 + treeDataDf$Edges_rs_false_negative_80), na.rm = TRUE)
accuracy_rs_60 <- mean((treeDataDf$Edges_rs_true_positive_60 + treeDataDf$"truenegatives") /(treeDataDf$Edges_rs_true_positive_60 + treeDataDf$Edges_rs_true_negative_60 + treeDataDf$Edges_rs_false_positive_60 + treeDataDf$Edges_rs_false_negative_60), na.rm = TRUE)
accuracy_rs_40 <- mean((treeDataDf$Edges_rs_true_positive_40 + treeDataDf$"truenegatives") /(treeDataDf$Edges_rs_true_positive_40 + treeDataDf$Edges_rs_true_negative_40 + treeDataDf$Edges_rs_false_positive_40 + treeDataDf$Edges_rs_false_negative_40), na.rm = TRUE)
accuracy_rs_20 <- mean((treeDataDf$Edges_rs_true_positive_20 + treeDataDf$"truenegatives") /(treeDataDf$Edges_rs_true_positive_20 + treeDataDf$Edges_rs_true_negative_20 + treeDataDf$Edges_rs_false_positive_20 + treeDataDf$Edges_rs_false_negative_20), na.rm = TRUE)

accuracy_rs <- c(accuracy_rs_100,accuracy_rs_80,accuracy_rs_60,accuracy_rs_40,accuracy_rs_20)

rm(list=ls())                                  # Clear Workspace to avoid data mistakes

list.of.packages <- c("rstudioapi",            # import find path script
                      "ggplot2",               # import plotting library
                      "ggpubr",                # Align plots in grid
                      "RColorBrewer",          # for beautiful colored 
                      "dplyr",                  # rename columns easely and more
                      "plotly"
                      
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(plotly)
pathScript <- dirname(rstudioapi::getSourceEditorContext()$path) # den Pfad des Scriptes finden
setwd(pathScript)  # set Working directory

dataDir <- paste(pathScript, "01_Data", sep="/")

treeDataDf <- read.table(paste(pathScript, 'Tree_Data_NEW.csv', sep = '/'), sep = ',', header = TRUE)
# treeDataDf$Group <- as.factor(strsplit(treeDataDf$ID, '_')[[1]][1])

i = 1
for (i in 1:length(treeDataDf$ID)) {
  treeDataDf$Group[i] <- strsplit(treeDataDf$ID[i], '_')[[1]][1]
}
treeDataDf$Group <- as.factor(treeDataDf$Group)

names(treeDataDf)

levels(treeDataDf$Group)

hist(treeDataDf$Fraction_of_Xenologs)

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

write.csv()


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

ggsave("Fraction_of_Xenologs_vs._Duplication_Rate.png", box_hgt_dupl)


box_hgt_loss <- ggplot(treeDataDf, aes(x = as.factor(loss_rate), y = Fraction_of_Xenologs, group=as.factor(hgt_rate))) +
  geom_jitter(shape=16, position=position_jitter(0.15), aes(color = factor(hgt_rate), group=factor(hgt_rate))) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, aes(alpha = 0.8, group=factor(dupl_rate))) + 
  labs(title = 'Fraction_of_Xenologs vs. Loss Rate', x ='Loss Rate', y = 'HGT Events', colour = 'HGT Rate') +
  theme_bw()
box_hgt_loss

ggsave("Fraction of Xenologs_vs._Loss_Rate.png", box_hgt_loss)

#### Signifikanzen ####

# Kruskal-Wallis-Test über alle Gruppen
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


