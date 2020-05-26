library(dplyr)
library(ggplot2)
library(gridExtra)

assign_colnames_bed <- function(df){
  colnames(df) = c("chrom", "chromStart", "chromEnd", "name", "score",
                   "strand", "thickStart", "thickEnd", "itemRgb")
  return(df)
}

fix_vir_coordinates <- function(df, vir_res){
  df[,c("chromStart", "chromEnd", "thickStart", "thickEnd")] = 
    df[,c("chromStart", "chromEnd", "thickStart", "thickEnd")] * vir_res
  return(df)
}

read_bedgraph <- function(file_path, vir_res){
  df <- read.table(file_path)
  df <- assign_colnames_bed(df)
  df <- fix_vir_coordinates(df, vir_res)
  df <- format(df, scientific = FALSE)
  df$chromStart = as.numeric(df$chromStart)
  df$chromEnd = as.numeric(df$chromEnd)
  return(df)
}

seq_annotation <- function(df, chr_name, start, end, resolution){
  df = filter(df, chrom == chr_name & start < chromEnd & chromStart < end)
  df[1,"chromStart"] = start
  df[nrow(df), "chromEnd"] = end
  seq_annot = rep(df$name, ceiling((df$chromEnd - df$chromStart)/resolution))
  return(seq_annot)
}

seq_plot <- function(seq_annot){
  ggplot() +
    geom_line(aes(x = seq(length(seq_annot)), y = 1, color = as.integer(seq_annot)),
              size = 6) + 
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_text(size=12, face = "bold"),
          legend.position="none", panel.background=element_blank(),
          panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
}

label_freq <- function(df, chr_name, start, end, resolution){
  df$chromStart = as.numeric(df$chromStart)
  df$chromEnd = as.numeric(df$chromEnd)
  start = as.numeric(start)
  end = as.numeric(end)
  num_labels = length(unique(df$name[!is.na(df$name)]))
  labels = seq(num_labels)
  names(labels) = as.vector(unique(df$name[!is.na(df$name)]))
  label_freq = rep(0, num_labels)
  df = filter(df, chrom == chr_name & chromEnd > start & chromStart < end)
  label_indices = labels[as.vector(df$name)]
  label_freq[label_indices[1]] = label_freq[label_indices[1]] + 
    (min(df$chromEnd[1], end) - start)/resolution
  if (nrow(df) > 1){
    label_freq[label_indices[nrow(df)]] = label_freq[label_indices[nrow(df)]] + 
      (end - df$chromStart[nrow(df)])/resolution
    if(nrow(df) > 2){
      for (i in c(2:(nrow(df)-1))){
        label_freq[label_indices[i]] = label_freq[label_indices[i]] + 
          (df$chromEnd[i] - df$chromStart[i])/resolution
      }
    }
  }
  return(label_freq)
}

enrichment <- function(df, chrom_name, start, end){
  df$chromStart = as.numeric(df$chromStart)
  df$chromEnd = as.numeric(df$chromEnd)
  start = as.numeric(start)
  end = as.numeric(end)
  df = filter(df, chrom == chr_name & chromEnd > start & chromStart < end)
  df[1,"chromStart"] = start
  if(nrow(df) == 1){
    df[1,"chromEnd"] = end
  }
}

#####make_ratio_table
ratio_table <- function(table){
  row_sum = apply(table, 1, FUN = sum)
  col_sum = apply(table, 2, FUN = sum)
  row_prop = row_sum / sum(row_sum)
  col_prop = col_sum / sum(col_sum)
  expected_prop = row_prop %*% t(col_prop)
  expected_table = expected_prop * sum(table)
  return(table/expected_table)
}

fold_enrichment <- function(df1, df2, resolution){
  labels_name = as.vector(unique(df2$name[!is.na(df2$name)]))
  num_labels = length(labels_name)
  df1[,10:(9+num_labels)] <- t(apply(df1[,1:3],1,
                                    function(x) label_freq(df2,x[1],
                                                           x[2],x[3],resolution)))
  colnames(df1)[10:(9+num_labels)] = labels_name
  enrichment_table = data.frame()
  df1_labels_name = as.vector(unique(df1$name[!is.na(df1$name)]))
  for (i in df1_labels_name){
    f_df1 = df1[df1$name==i & !is.na(df1$name),]
    domain_freq = apply(f_df1[,10:(9+num_labels)], 2, function(x) sum(x))
    enrichment_table = rbind(enrichment_table, domain_freq)
  }
  rownames(enrichment_table) = df1_labels_name
  colnames(enrichment_table) = labels_name
  enrichment_table = round(ratio_table(enrichment_table), 3)
  return(enrichment_table)
}

HiC_mat <- function(COO, resolution){
  COO_data <- read.table(COO)
  COO_data$V1 = COO_data$V1 / resolution
  COO_data$V2 = COO_data$V2 / resolution
  print(head(COO_data))
  HiC_mat <- as.matrix(sparseMatrix(i=COO_data[,1], j=COO_data[,2],
                                      x=COO_data[,3], symmetric = TRUE))
  HiC_mat[is.na(HiC_mat)] = 0
  return(HiC_mat)
}

OtoE_HiC_mat <- function(COO, resolution){
  COO_data <- read.table(COO)
  COO_data$V1 = COO_data$V1 / resolution
  COO_data$V2 = COO_data$V2 / resolution
  expected_contact <- COO_data %>% group_by(V2 - V1) %>% 
    summarise(expected_contact = mean(V3, na.rm = T))
  exp_contact <- expected_contact$expected_contact
  names(exp_contact) = expected_contact$`V2 - V1`
  COO_data$OtoE = COO_data$V3 / exp_contact[as.character(COO_data$V2 - COO_data$V1)]
  print(head(COO_data))
  OtoE_mat <- as.matrix(sparseMatrix(i=COO_data[,1], j=COO_data[,2],
                                           x=COO_data[,4], symmetric = TRUE))
  OtoE_mat[is.na(OtoE_mat)] = 0
  return(OtoE_mat)
}
