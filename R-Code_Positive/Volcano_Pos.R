#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

# Read the data into R
dset<-read.csv("~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/12102022_Normalised_Quant_pos.csv", header = TRUE)
data<- as.matrix(dset[,2:72])
featureIDs <- dset[,1]
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames
ko = data[,1:19]
wt_lim = data[,20:35]
wt = data[,36:55]
a7 = data[,56:71]

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
wt_lim.mean = apply(wt_lim, 1, mean)
a7.mean = apply(a7, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, ko.mean, wt_lim.mean, a7.mean)
means <- cbind.data.frame(wt.mean, ko.mean, wt_lim.mean, a7.mean)

###############################################################
######### WT VS KO ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / ko.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = ko[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Volcano Plot of Features"){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             "Significant, downregulated", 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    "Significant, upregulated", "Normally expressed")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

install.packages('svglite')
ggsave(file="~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/pos_volcano.svg", plot=plot, width=10, height=8)
write.csv(DEres.df, "~/Desktop/volcano_t-test_pos.csv", row.names=FALSE)


############################################################################
#####box plot with wilcoxon and krustal-wallis stats for feature m/z 797####
library("ggpubr")
library("datatable")
library("dplyr")
library("tidyverse")
library("readr")
my_data2 <- read.csv("~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/Normalized_Imp_Quant_799.csv")
##note: this was created from copying feature ID 10001_797.404_8.868 into the new table above##
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = .0012)
fig
ggsave(file="wilcox_boxplot.svg", plot=fig, width=10, height=8)

############################################################################
#####box plot with anova and stats for feature m/z 797####

library("ggpubr")
library("data.table")
my_data2 <- read.csv("~/Dropbox/2023-Lanthanide_Nature/R-Code_Negative/Volcano_Inputs/Normalized_Imp_Quant_797.csv")
##note: this was created from copying feature ID 10001_797.404_8.868 into the new table above##
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(method = "anova", label.y = .0012)
fig

ggsave(file="anova_boxplot.svg", plot=fig, width=10, height=8)
