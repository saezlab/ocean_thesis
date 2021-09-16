
metactivity <- function(metabolomic_t_table, regulons_df, compartment_pattern = "", k = 10000)
{
  enzymes_nes_list <- list()
  enzymes_es_list <- list()
  for(i in 2:length(metabolomic_t_table[1,]))
  {

    if(compartment_pattern != "")
    {
      t_table <- metabolomic_t_table
      t_table$unique_metab <- gsub(compartment_pattern, "", metabolomic_t_table[,1])
      t_table_reduced <- t_table[,c(length(t_table[1,]),2:(length(t_table[1,])-1))]
      t_table_reduced <- unique(t_table_reduced)
      
      t_null <- as.data.frame(replicate(k, sample(t_table_reduced[,i], length(t_table_reduced[,i]))))
      t_null <- as.data.frame(cbind(t_table_reduced[,1],t_null))
      names(t_null)[1] <- names(t_table_reduced)[1]
      
      t_table_with_null <- merge(t_table[,c(1,i,length(t_table[1,]))], 
                                 t_null, 
                                 by = names(t_table_reduced)[1])
      t_table_with_null <- t_table_with_null[,-1]
      row.names(t_table_with_null) <- t_table_with_null[,1]
      t_table_with_null <- t_table_with_null[,-1]
    }
    
    regulons_df <- regulons_df[regulons_df[,2] %in% row.names(t_table_with_null),]

    regulons_mat <- dcast(regulons_df, targets~set, value.var = "weight")
    row.names(regulons_mat) <- regulons_mat[,1]
    regulons_mat <- regulons_mat[,-1]
    
    t_table_with_null <- t_table_with_null[row.names(t_table_with_null) %in% row.names(regulons_mat),]

    regulons_mat <- regulons_mat[row.names(t_table_with_null),]
    regulons_mat[is.na(regulons_mat)] <- 0 
    
    
    metabolites <- row.names(t_table_with_null)
    enzymes <- names(regulons_mat)
    t_table_with_null <- t(t_table_with_null)
    
    enzyme_ES <- as.matrix(t_table_with_null) %*% as.matrix(regulons_mat) 
    
    enzymes_es_list[[i]] <- enzyme_ES[1,]
    
    enzyme_NES <- enzyme_ES[1,]
    null_means <- colMeans(enzyme_ES[-1,])
    null_sds <- apply(enzyme_ES[-1,],2,sd)
    
    enzyme_NES <- (enzyme_NES - null_means) / null_sds
    
    enzymes_nes_list[[i]] <- enzyme_NES
  }
  
  enzymes_nes_df <- as.data.frame(do.call(cbind,enzymes_nes_list))
  enzymes_nes_df$enzyme <- row.names(enzymes_nes_df)
  enzymes_nes_df <- enzymes_nes_df[,c(length(enzymes_nes_df[1,]),1:(length(enzymes_nes_df[1,])-1))]
  names(enzymes_nes_df) <- names(metabolomic_t_table)
  
  enzymes_es_df <- as.data.frame(do.call(cbind,enzymes_es_list))
  enzymes_es_df$enzyme <- row.names( enzymes_es_df)
  enzymes_es_df <-  enzymes_es_df[,c(length( enzymes_es_df[1,]),1:(length( enzymes_es_df[1,])-1))]
  names( enzymes_es_df) <- names(metabolomic_t_table)
  
  return(list("ES" = enzymes_es_df, "NES" = enzymes_nes_df))
}

target_set_from_forest <- function(forest, measured_species, penalty = 0.5)
{
  target_set <- as.data.frame(matrix(NA, 1, 3))
  names(target_set) <- c("set","targets","weight")
  
  i <- 1
  for(tree in forest)
  {
    print(paste0(i, "/", length(names(forest)), " ", "pen = ",penalty))
    direction_sign <- -1
    #print(names(tree))
    for (direction in tree)
    {
      j <- 1
      for (layer in direction)
      {
        if(length(layer) > 0)
        {
          layer_set <- as.data.frame(matrix(NA,length(layer), 3))
          names(layer_set) <- c("set","targets","weight")
          
          layer_set$set <- names(forest)[i]
          layer_set$targets <- layer
          layer_set$weight <- j*direction_sign
          
          j <- j*penalty
          
          layer_set <- layer_set[layer_set$targets %in% measured_species,]
          
          if(length(layer_set[,1] > 0))
          {
            target_set <- rbind(target_set, layer_set)
          }
        }
      }
      direction_sign <- direction_sign*-1
    }
    i <- i+1
  }
  
  return(target_set[-1,])
}

target_set_from_forest_2 <- function(forest, measured_species, penalty = 0.5)
{
  # target_set <- as.data.frame(matrix(NA, 1, 3))
  # names(target_set) <- c("set","targets","weight")
  # 
  target_set <- list()
  i <- 1
  l <- 1
  for(tree in forest)
  {
    print(paste0(i, "/", length(names(forest)), " ", "pen = ",penalty))
    direction_sign <- -1
    #print(names(tree))
    for (direction in tree)
    {
      j <- 1
      for (layer in direction)
      {
        if(length(layer) > 0)
        {
          layer_set <- as.data.frame(matrix(NA,length(layer), 3))
          names(layer_set) <- c("set","targets","weight")
          
          layer_set$set <- names(forest)[i]
          layer_set$targets <- layer
          layer_set$weight <- j*direction_sign
          
          j <- j*penalty
          
          layer_set <- layer_set[layer_set$targets %in% measured_species,]
          
          if(length(layer_set[,1] > 0))
          {
            target_set[[l]] <- layer_set
            l <- l+1
          }
        }
      }
      direction_sign <- direction_sign*-1
    }
    i <- i+1
  }
  target_set <- as.data.frame(do.call(rbind,target_set))
  return(target_set)
}

ttop_list_to_t_table <- function(ttop_list)
{
  if(length(ttop_list) > 1)
  {
    t_table <- merge(ttop_list[[1]][,c(1,4)], ttop_list[[2]][,c(1,4)], by = "ID")
    if(length(ttop_list) > 2)
    {
      for(i in 3:length(ttop_list))
      {
        t_table <- merge(t_table, ttop_list[[i]][,c(1,4)], by = "ID")
      }
    }
  }
  else
  {
    t_table <- ttop_list[[1]][,c(1,4)]
  }
  names(t_table) <- c("ID",names(ttop_list))
  return(t_table)
}

ttop_list_to_log2FC_table <- function(ttop_list)
{
  if(length(ttop_list) > 1)
  {
    log2FC_table <- merge(ttop_list[[1]][,c(1,2)], ttop_list[[2]][,c(1,2)], by = "ID")
    if(length(ttop_list) > 2)
    {
      for(i in 3:length(ttop_list))
      {
        log2FC_table <- merge(log2FC_table, ttop_list[[i]][,c(1,2)], by = "ID")
      }
    }
  }
  else
  {
    log2FC_table <- ttop_list[[1]][,c(1,2)]
  }
  names(log2FC_table) <- c("ID",names(ttop_list))
  return(log2FC_table)
}

limma_res_to_ttop_list <- function(limma_res, comp_names, number, adjust.method = "fdr")
{
  ttop_list <- list()
  n_comp <- length(limma_res[[2]][1,])
  for(i in 1:n_comp)
  {
    ttop_list[[i]] <- ttopFormatter(topTable(limmaRes[[1]], coef = i, number = number, adjust.method = adjust.method))
    ttop_list[[i]] <- ttop_list[[i]][complete.cases(ttop_list[[i]]),]
  }
  names(ttop_list) <- comp_names
  return(ttop_list)
}

t_table_metactivity_input_formater <- function(metabolomic_t_table, mapping_table, affixes = c("c","l","x","m","e","n","r"))
{
  names(metabolomic_t_table)[1] <- "metabolite"
  metabolomic_t_table[,1] <- gsub(",","_",metabolomic_t_table[,1])
  metabolomic_t_table[,1] <- gsub(" ","",metabolomic_t_table[,1])
  
  names(mapping_table)[1] <- "metabolite" 
  mapping_table$metabolite <- gsub(",","_",mapping_table$metabolite)
  mapping_table$metabolite <- gsub(" ","",mapping_table$metabolite)
  
  mapping_table$KEGG <- paste("cpd:",mapping_table$KEGG, sep = "")
  
  ####WITH RECON
  affixes <- c("c","l","x","m","e","n","r")
  
  temp_mapping_table <- mapping_table
  
  for (aff in affixes)
  {
    new_mapping_table <- temp_mapping_table
    new_mapping_table$KEGG <- paste(new_mapping_table$KEGG, aff, sep = "_")
    mapping_table <- as.data.frame(rbind(mapping_table, new_mapping_table))
  }
  
  mapping_table[,1] <- tolower(mapping_table[,1])
  metabolomic_t_table$metabolite <- tolower(metabolomic_t_table$metabolite)
  
  metabolomic_t_table <- merge(metabolomic_t_table, mapping_table, by = "metabolite")   
  table_length <- length(metabolomic_t_table[1,])
  metabolomic_t_table <- metabolomic_t_table[,c(table_length,2:(table_length-1))]
  
  return(metabolomic_t_table)
}

prepare_metabolite_set <- function(penalty_range, reaction_tree, metab_list)
{
  reaction_set_list <- list()
  for(i in penalty_range)
  {
    print(i)
    pen <- i/10
    reaction_set_list[[i]] <- target_set_from_forest_2(reaction_tree, metab_list, penalty = pen)
  }
  return(reaction_set_list)
}

plotMetaboliteContribution <- function(enzyme, stat_df, metabolite_sets, contrast_index, stat_name = "stat", scaling_factor = 1, nLabels = 10)
{
  
  metabolite_sets <- metabolite_sets[metabolite_sets$set == enzyme,]
  
  stat_df <- stat_df[stat_df[,1] %in% metabolite_sets$targets,]
  stat_df <- stat_df[,c(1,contrast_index+1)]
  
  
  stat_df <- unique(stat_df)
  
  names(stat_df) <- c("ID","C1")
  names(metabolite_sets) <- c("set","ID","weight")
  
  df <- merge(metabolite_sets,stat_df, by = "ID")
  
  df$contribution <- df$weight*df$C1*scaling_factor
  # df$contribution <- rescale(df$contribution, to = c(-3,3))
  df <- df[order(abs(df$contribution),decreasing = T),]

  save_labels <- df$ID  
  
  if(length(df$ID) > nLabels)
  {
    df$ID <- c(df[1:nLabels,'ID'],rep(NA,length(df$ID) - nLabels)) 
  }

  # df$ID <- ifelse(abs(df$contribution) >=  1,df$ID,"")
  # df$ID <- ifelse(abs(df$C1) >=  2,df$ID,"")
  
  ggscat <- ggplot(df, aes(x = weight, y = C1, label = ID)) + 
    geom_point(aes(fill = contribution),colour = "black",pch=21, size = abs(df$contribution)) +
    xlim(c(-1,1)) +
    geom_text_repel(point.padding = unit(max(abs(df$contribution)), "points")) +
    # scale_color_scico(palette = "vik", limits = c(-1, 1) * max(abs(df$contribution))) +
    scale_fill_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
    theme_minimal() + 
    geom_abline(intercept = 0, slope = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    ylim(c(-max(abs(df$C1)),max(abs(df$C1)))) + 
    ggtitle(paste0(enzyme," metabolic consumption/production profile")) +
    labs(x = "consuption <==> production", y = stat_name)
  
  df <- df[order(df$contribution,decreasing = F),]
  df$running_sum_contribution <- cumsum(df$contribution)
  
  df$ID <- save_labels
  df$ID <- paste0(c(1:length(df$ID)), " ", df$ID)
  df$ID <- factor(df$ID, level = df$ID)
  df$group <- "A"
  df$label <- gsub("[0-9]+ ","",df$ID)
  
  ylim_bot <- ifelse(min(df$running_sum_contribution) < 0, min(df$running_sum_contribution)*1.25, 0 ) 
  ylim_bot <- ifelse(ylim_bot < min(df$contribution), ylim_bot, min(df$contribution) * 1.25)
  
  ylim_top <- ifelse(max(df$running_sum_contribution) > 0, max(df$running_sum_contribution)*1.25, 0 ) 
  ylim_top <- ifelse(ylim_top > max(df$contribution), ylim_top, max(df$contribution) * 1.25)
  
  ggcumSum <- ggplot(df, aes(x = ID, y = running_sum_contribution, group = group, label = label)) + 
    geom_point(aes(x = ID, y = contribution)) +
    geom_hline(yintercept = 0) +
    geom_line() + 
    theme_minimal() + 
    ylim(c(ylim_bot,ylim_top)) +
    scale_x_discrete(guide = guide_axis(n.dodge=8)) +
    theme(plot.margin = unit(c(2,2,2,2), "cm"))
  
  return(list("scatter" = ggscat, "cumsumPlot" = ggcumSum))
}

pathway_HM <- function(mean_NES_df, pathway_name, pathways, sort_by = 1, manual_pathway = F)
{
  current_pathway <- pathways[pathways$term == pathway_name,"gene"]
  
  mean_NES_df_to_subset <- mean_NES_df
  
  if(manual_pathway == F)
  {
    mean_NES_df_to_subset$ID_short <- gsub(">.*","",mean_NES_df$KEGG)
    mean_NES_df_to_subset$ID_short <- gsub("_.*","",mean_NES_df_to_subset$ID_short)
  } else
  {
    mean_NES_df_to_subset$ID_short <- mean_NES_df_to_subset$KEGG
  }
  
  current_pathway_activities <- mean_NES_df[mean_NES_df_to_subset$ID_short %in% current_pathway,]
  
  row.names(current_pathway_activities) <- current_pathway_activities$KEGG
  current_pathway_activities <- current_pathway_activities[order(current_pathway_activities[,sort_by+1], decreasing = F),]
  
  
  current_pathway_act_melt <- melt(current_pathway_activities)
  current_pathway_act_melt$KEGG <- factor(current_pathway_act_melt$KEGG,levels = unique(current_pathway_act_melt$KEGG))
  names(current_pathway_act_melt) <- c("Enzyme","Contrast","NES")
  
  gp <- ggplot(current_pathway_act_melt, aes(x=Contrast, y = Enzyme, fill = NES)) + geom_tile() +
    scale_fill_gradient2(low="blue", high="red", midpoint = 0, mid = "white") + 
    theme_minimal() +
    labs(y = "Enzyme", x = "Contrast") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(gp)
}
