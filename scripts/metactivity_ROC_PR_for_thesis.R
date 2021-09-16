library(readxl)
library(vsn)
library(limma)
library(readr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(scico)
library(scales)

source("scripts/metactivity_function.R")
source("scripts/support_functions.R") 

PLASMAX <- as.data.frame(
  read_excel("data/Metabolomics Plasmax N=3.xlsx", 
             sheet = "Sheet2"))

targets <- PLASMAX[,c(2,1)]
names(targets) <- c("sample","condition")
targets$sample <- paste(targets$condition,c(1:length(PLASMAX[,1])), sep = "_")

batches <- as.data.frame(t(PLASMAX[,-1]))
names(batches) <- targets$sample
row.names(batches) <- gsub(" ", "_", row.names(batches))
row.names(batches) <- tolower(row.names(batches))

batches <- as.data.frame(apply(batches,2,function(x){as.numeric(as.character(x))}))
batches[batches == 0] <- NA
batches <- log2(batches)

row.names(batches) <- names(PLASMAX[,-1])

comparisons <- list("OvHK2" = c(2,-1), 
                    "M1AvHK2" = c(3,-1), 
                    "M2AvHK2" = c(4,-1),
                    "M1AvO" = c(3,-2))

limmaRes <- runLimma(batches, targets, comparisons = comparisons)

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes, 
                         comp_names = names(comparisons),
                         number = length(batches[,1]),
                         adjust.method = "fdr"))

to_dotplot <- t_table
to_dotplot <- to_dotplot[-45,] #remove itaconitate misslabeled
to_dotplot <- to_dotplot[order(to_dotplot$OvHK2, decreasing = F),]
to_dotplot_long <- melt(to_dotplot[,c(1,2,3)])

to_dotplot_long$ID <- factor(to_dotplot_long$ID, levels = unique(to_dotplot_long$ID))

ggplot(to_dotplot_long, aes(x = value, y = ID, size = abs(value))) +
  geom_point(aes(colour = variable)) +
  # scale_color_gradient2(low="blue", high="blue", midpoint = 0, mid = "blue") +
  theme_minimal() + ggtitle("786-0 and 786-M1A versus HK2 metabolomic") + 
  geom_vline(xintercept = -2) +
  geom_vline(xintercept = 2) 

to_dotplot$ID <- factor(to_dotplot$ID, levels = to_dotplot$ID)

ggplot(to_dotplot, aes(x = M1AvO, y = ID, size = abs(M1AvO))) +
  geom_point(aes(colour = M1AvO)) +
  # scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  scale_color_gradient2(low="blue", high="blue", midpoint = 0, mid = "blue") +
  theme_minimal() + ggtitle("786-M1A versus 786-O metabolomic") + 
  geom_vline(xintercept = -2) +
  geom_vline(xintercept = 2) 

mapping_table <- as.data.frame(read_csv("support/plasmax_name_to_kegg.txt"))

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))

load("support/reaction_set_list_merged_SUC.Rdata")
load("support/reaction_network_ASS1_SUC_no_cofac.Rdata")

##########################
metactivity_list <- list()
for(i in 1:10)
{
  print(i)
  penalty <- i
  
  n <- length(unique(reaction_set_list_merged[[penalty]][,1]))
  
  regulonNames = unique(reaction_set_list_merged[[penalty]][,1])[1:n]
  regulons_df <- reaction_set_list_merged[[penalty]]
  
  regulons_df <- regulons_df[regulons_df[,1] %in% unique(regulons_df[,1])[1:n],]
  
  metactivity_list[[i]] <- metactivity(metabolomic_t_table = t_table, 
                                       regulons_df = regulons_df, 
                                       compartment_pattern = "_[a-z]$", 
                                       k = 1000)
}

ttop_786vHK2 <- as.data.frame(read_csv("data/proteo/ttop_786vHK2.csv"))
ttop_M1AvHK2 <- as.data.frame(read_csv("data/proteo/ttop_M1AvHK2.csv"))
ttop_M1Av786 <- as.data.frame(read_csv("data/proteo/ttop_M1Av786.csv"))

cor_coef_list <- list()
c786_yardstick_list <- list()
M1A_yardstick_list <- list()
M1Av786_yardstick_list <- list()
for(i in 1:10)
{
  metactivity <- metactivity_list[[i]]$NES
  # metactivity <- metactivity[!grepl("^SLC",metactivity$KEGG),]
  metactivity$KEGG <- gsub(">.*","",metactivity$KEGG)
  
  metactivity <- metactivity %>% group_by(KEGG) %>% summarise_each(funs(mean(., na.rm = TRUE)))
  metactivity <- as.data.frame(metactivity)
  
  metactivty_786 <- metactivity[,c(1,2)]
  names(metactivty_786) <- c("ID","NES")
  metactivty_M1A <- metactivity[,c(1,3)]
  names(metactivty_M1A) <- c("ID","NES")
  metactivty_M1Av786 <- metactivity[,c(1,4)]
  names(metactivty_M1Av786) <- c("ID","NES")
  
  c786_comp <- merge(ttop_786vHK2[,c(1,4)],metactivty_786,by = "ID")
  c786_comp_yardstick <- c786_comp
  c786_comp_yardstick$penalty <- paste("p",i,sep = "")
  # c786_comp <- c786_comp[abs(c786_comp$NES) > 0.3,]
  
  M1A_comp <- merge(ttop_M1AvHK2[,c(1,4)],metactivty_786,by = "ID")
  M1A_comp_yardstick <- M1A_comp
  M1A_comp_yardstick$penalty <- paste("p",i,sep = "")
  # M1A_comp <- M1A_comp[abs(M1A_comp$NES) > 0.3,]
  
  M1Av786_comp <- merge(ttop_M1Av786[,c(1,4)],metactivty_786,by = "ID")
  M1Av786_comp_yardstick <- M1Av786_comp
  M1Av786_comp_yardstick$penalty <- paste("p",i,sep = "")
  # M1Av786_comp <- M1Av786_comp[abs(M1Av786_comp$NES) > 0.3,]
  
  c786_yardstick_list[[i]] <- c786_comp_yardstick
  M1A_yardstick_list[[i]] <- M1A_comp_yardstick
  M1Av786_yardstick_list[[i]] <- M1Av786_comp_yardstick
  
  c786_comp_cor <- cor.test(as.numeric(c786_comp[,2]),as.numeric(c786_comp[,3]), method = "spearman")$estimate
  M1A_comp_cor <- cor.test(as.numeric(M1A_comp[,2]),as.numeric(M1A_comp[,3]), method = "spearman")$estimate
  M1Av786_comp_cor <- cor.test(as.numeric(M1Av786_comp[,2]),as.numeric(M1Av786_comp[,3]), method = "spearman")$estimate
  
  cor_coef_list[[i]] <- c(c786_comp_cor, M1A_comp_cor, M1Av786_comp_cor)
}

cor_coef_df <- as.data.frame(do.call(rbind,cor_coef_list))

##########YARDSTICK#########
library(yardstick)

###c786 split up
c786_yardstick_df <- as.data.frame(do.call(rbind,c786_yardstick_list))
c786_yardstick_df <- c786_yardstick_df[,-1]
names(c786_yardstick_df)[1] <- "t_val"
# c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$NES > 0,]
c786_yardstick_df$t_val <- ifelse(c786_yardstick_df$t_val < 1.7, 0, c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- sign(c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- abs(c786_yardstick_df$t_val)
# c786_yardstick_df$NES <- abs(c786_yardstick_df$NES)
c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$penalty != "p10",]

c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_auc(t_val, NES)
c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_auc(t_val, NES)
pr_baseline <- length(c786_yardstick_df[c786_yardstick_df$t_val == 1,1])/length(c786_yardstick_df[,1])
pr_baseline

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_curve(t_val, NES) 

ggplot(to_ggplot, aes(x = 1 - specificity, y = sensitivity, color = penalty)) + geom_path(size = 1) + geom_abline(intercept = 0, slope = 1) +
  theme_minimal()

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_curve(t_val, NES) 
to_ggplot <- to_ggplot[complete.cases(to_ggplot),]
pr_baseline <- length(c786_yardstick_df[c786_yardstick_df$t_val == 1,1])/length(c786_yardstick_df[,1])

ggplot(to_ggplot, aes(x = recall, y = precision, color = penalty)) + geom_path(size = 1) + coord_equal() + geom_abline(intercept = pr_baseline, slope = 0) +
  theme_minimal()

###c786 split down
c786_yardstick_df <- as.data.frame(do.call(rbind,c786_yardstick_list))
c786_yardstick_df <- c786_yardstick_df[,-1]
names(c786_yardstick_df)[1] <- "t_val"
# c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$NES < 0,]
c786_yardstick_df$t_val <- ifelse(c786_yardstick_df$t_val > -1.7, 0, c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- sign(c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- abs(c786_yardstick_df$t_val)
c786_yardstick_df$NES <- c786_yardstick_df$NES * -1
c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$penalty != "p10",]

c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_auc(t_val, NES)
c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_auc(t_val, NES)
pr_baseline <- length(c786_yardstick_df[c786_yardstick_df$t_val == 1,1])/length(c786_yardstick_df[,1])
pr_baseline


to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_curve(t_val, NES) 

ggplot(to_ggplot, aes(x = 1 - specificity, y = sensitivity, color = penalty)) + geom_path(size = 1) + geom_abline(intercept = 0, slope = 1) +
  theme_minimal()

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_curve(t_val, NES) 
to_ggplot <- to_ggplot[complete.cases(to_ggplot),]

pr_baseline <- mean(as.numeric(unlist(to_ggplot[to_ggplot$recall == 1,4])))
pr_baseline

to_ggplot <- to_ggplot[to_ggplot$recall != 0 & to_ggplot$precision != 0,]

ggplot(to_ggplot, aes(x = recall, y = precision, color = penalty)) + geom_path(size = 1) + coord_equal() + geom_abline(intercept = pr_baseline, slope = 0) +
  theme_minimal()


############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

ttop_786vHK2 <- as.data.frame(read_csv("data/RNAseq_ttop_metabs_OvHK2.csv"))

cor_coef_list <- list()
c786_yardstick_list <- list()
M1A_yardstick_list <- list()
M1Av786_yardstick_list <- list()
for(i in 1:10)
{
  metactivity <- metactivity_list[[i]]$NES
  # metactivity <- metactivity[!grepl("^SLC",metactivity$KEGG),]
  metactivity$KEGG <- gsub(">.*","",metactivity$KEGG)
  
  metactivity <- metactivity %>% group_by(KEGG) %>% summarise_each(funs(mean(., na.rm = TRUE)))
  metactivity <- as.data.frame(metactivity)
  
  metactivty_786 <- metactivity[,c(1,2)]
  names(metactivty_786) <- c("ID","NES")
  metactivty_M1A <- metactivity[,c(1,3)]
  names(metactivty_M1A) <- c("ID","NES")
  metactivty_M1Av786 <- metactivity[,c(1,4)]
  names(metactivty_M1Av786) <- c("ID","NES")
  
  c786_comp <- merge(ttop_786vHK2[,c(1,4)],metactivty_786,by = "ID")
  c786_comp_yardstick <- c786_comp
  c786_comp_yardstick$penalty <- paste("p",i,sep = "")
  c786_comp <- c786_comp[abs(c786_comp$NES) > 2,]
  
  M1A_comp <- merge(ttop_M1AvHK2[,c(1,4)],metactivty_786,by = "ID")
  M1A_comp_yardstick <- M1A_comp
  M1A_comp_yardstick$penalty <- paste("p",i,sep = "")
  M1A_comp <- M1A_comp[abs(M1A_comp$NES) > 2,]
  
  M1Av786_comp <- merge(ttop_M1Av786[,c(1,4)],metactivty_786,by = "ID")
  M1Av786_comp_yardstick <- M1Av786_comp
  M1Av786_comp_yardstick$penalty <- paste("p",i,sep = "")
  M1Av786_comp <- M1Av786_comp[abs(M1Av786_comp$NES) > 2,]
  
  c786_yardstick_list[[i]] <- c786_comp_yardstick
  M1A_yardstick_list[[i]] <- M1A_comp_yardstick
  M1Av786_yardstick_list[[i]] <- M1Av786_comp_yardstick
  
  c786_comp_cor <- cor.test(as.numeric(c786_comp[,2]),as.numeric(c786_comp[,3]), method = "spearman")$estimate
  M1A_comp_cor <- cor.test(as.numeric(M1A_comp[,2]),as.numeric(M1A_comp[,3]), method = "spearman")$estimate
  M1Av786_comp_cor <- cor.test(as.numeric(M1Av786_comp[,2]),as.numeric(M1Av786_comp[,3]), method = "spearman")$estimate
  
  cor_coef_list[[i]] <- c(c786_comp_cor, M1A_comp_cor, M1Av786_comp_cor)
}

cor_coef_df <- as.data.frame(do.call(rbind,cor_coef_list))

##########YARDSTICK#########

###c786 split up
c786_yardstick_df <- as.data.frame(do.call(rbind,c786_yardstick_list))
c786_yardstick_df <- c786_yardstick_df[,-1]
names(c786_yardstick_df)[1] <- "t_val"
# c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$NES > 0,]
c786_yardstick_df$t_val <- ifelse(c786_yardstick_df$t_val < 1.7, 0, c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- sign(c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- abs(c786_yardstick_df$t_val)
# c786_yardstick_df$NES <- abs(c786_yardstick_df$NES)
c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$penalty != "p10",]

c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_auc(t_val, NES)
c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_auc(t_val, NES)

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_curve(t_val, NES) 

ggplot(to_ggplot, aes(x = 1 - specificity, y = sensitivity, color = penalty)) + geom_path(size = 3) + geom_abline(intercept = 0, slope = 1) +
  theme_minimal()

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_curve(t_val, NES) 
to_ggplot <- to_ggplot[complete.cases(to_ggplot),]

pr_baseline <- mean(as.numeric(unlist(to_ggplot[to_ggplot$recall == 1,4])))
pr_baseline

ggplot(to_ggplot, aes(x = recall, y = precision, color = penalty)) + geom_path(size = 3) + coord_equal() + geom_abline(intercept = pr_baseline, slope = 0) +
  theme_minimal()

###c786 split down
c786_yardstick_df <- as.data.frame(do.call(rbind,c786_yardstick_list))
c786_yardstick_df <- c786_yardstick_df[,-1]
names(c786_yardstick_df)[1] <- "t_val"
# c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$NES < 0,]
c786_yardstick_df$t_val <- ifelse(c786_yardstick_df$t_val > -1.7, 0, c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- sign(c786_yardstick_df$t_val)
c786_yardstick_df$t_val <- abs(c786_yardstick_df$t_val)
c786_yardstick_df$NES <- c786_yardstick_df$NES * -1
c786_yardstick_df <- c786_yardstick_df[c786_yardstick_df$penalty != "p10",]

c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_auc(t_val, NES)
c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_auc(t_val, NES)



to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% roc_curve(t_val, NES) 

ggplot(to_ggplot, aes(x = 1 - specificity, y = sensitivity, color = penalty)) + geom_path(size = 3) + geom_abline(intercept = 0, slope = 1) +
  theme_minimal()

to_ggplot <- c786_yardstick_df %>% as_tibble() %>% group_by(penalty) %>% mutate(t_val = factor(as.numeric(t_val), levels = c(1,0))) %>% pr_curve(t_val, NES) 
to_ggplot <- to_ggplot[complete.cases(to_ggplot),]

pr_baseline <- mean(as.numeric(unlist(to_ggplot[to_ggplot$recall == 1,4])))
pr_baseline

to_ggplot <- to_ggplot[to_ggplot$recall != 0 & to_ggplot$precision != 0,]

ggplot(to_ggplot, aes(x = recall, y = precision, color = penalty)) + geom_path(size = 3) + coord_equal() + geom_abline(intercept = pr_baseline, slope = 0) +
  theme_minimal()

###
