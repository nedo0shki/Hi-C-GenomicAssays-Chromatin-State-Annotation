early_rep <- read.table("../../data/RT/GM12878-RT.bedGraph")
late_rep <- read.table("../../data/RT/GM12878_RT_late.bedGraph")
EL_rep <- early_rep
colnames(EL_rep) = c("chrom", "chromStart", "chromEnd", "E")
EL_rep$L = late_rep$V4
rm(early_rep, late_rep)
EL_rep_seq = read.table("../../data/RT/RT_GM12878_Lymphocyte_Int90901931_hg19.bedgraph")
noise = 0.01
EL_rep$log_ratio = log(((EL_rep$E+noise)/(EL_rep$L+noise)),2)
EL_rep$log_ratio2 = log(((EL_rep$E)/(EL_rep$L)),2)
colnames(EL_rep_seq) = c("chrom", "chromStart", "chromEnd", "EL_log_ratio")
EL_rep_seq$label = apply(EL_rep_seq[,1:3],1,
                          function(x) get_domain(segway_func_assays,x[1],x[2],x[3]))
EL_rep_seq$label = as.numeric(EL_rep_seq$label)
EL_rep_seq$label_agg = apply(EL_rep_seq[,1:3],1,
                         function(x) get_domain(segway_aggregation,x[1],x[2],x[3]))
EL_rep_seq$label_agg = as.numeric(EL_rep_seq$label_agg)
library(ggplot2)
ggplot(EL_rep_seq,aes(x=EL_log_ratio))+geom_histogram()+facet_grid(~label)+theme_bw()
ggplot(EL_rep_seq,aes(x=EL_log_ratio))+geom_histogram()+facet_grid(~label_agg)+theme_bw()
get_domain <- function(df, chrom_name, start, end){
  coord = mean(as.numeric(start), as.numeric(end))
  df = filter(df, chrom == chrom_name & chromStart <= coord & coord < chromEnd)
  return(df$name[1])
}