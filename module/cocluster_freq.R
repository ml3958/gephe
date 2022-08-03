#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This is a highly specialized script, espetially the plotting function

suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(NbClust))
# install.packages("tidyinftheo")
suppressMessages(library(infotheo))
suppressMessages(library(optparse))


option_list = list(
  # FILE INPUT
  make_option(c( "--dir_pog"), type="character", default=NULL),
  make_option(c( "--dir_module"), type="character", default=NULL),
  make_option(c( "--prefix_pog_pp"), type="character", default=NULL, help="prefix for pog phylogenetic profiles"),
  
  # COCLUSTER PARAMS (DEFAULT)
  make_option(c("--k_min"), type="double", default=0.1, help="cocluster k minimum (percentile)" ),
  make_option(c("--k_max"), type="double", default=0.6,help="cocluster k maximum (percentile)"),
  make_option(c("--k_step"), type="double", default=0.1,help="cocluster k maximum (percentile)"),
  make_option(c("--k_times"), type="integer", default=20, help="times for perform clustering for each k")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
dir_pog = opt$dir_pog
dir_module = opt$dir_module
prefix_pog_pp = opt$`prefix_pog_pp`
k_min = opt$k_min
k_max = opt$k_max
k_step = opt$k_step
k_times = opt$k_times


get_cocluster_freq = function(pog_pp, k_min, k_max, k_step, k_times){
  freq = matrix(0,ncol=nrow(pog_pp),nrow=nrow(pog_pp))
  rownames(freq) =colnames(freq) = row.names(pog_pp)
  
  k_range = round(quantile(0:nrow(pog_pp), seq(k_min,k_max,k_step)))
  print(paste0('Kmeans with K=: ', paste(k_range,collapse = ","),'; Repeat each K for: ', k_times, ' times'))
  k_range = rep(k_range,k_times)
  
  for (k in k_range){
    cluster = kmeans(pog_pp,k)
    for (i in unique(cluster$cluster)){
      member = names(cluster$cluster[cluster$cluster==i])
      freq[member,member] = freq[member,member] + 1
    }
  }
  freq = freq/length(k_range)
  return(freq)
}

dist_fun = function(m) {
  return( as.dist(1-as.matrix(m)))
}

# [input]
file_pog_pp = paste0(dir_pog, prefix_pog_pp, '.pp')

# [output]
fileout_pog_pp_coclusterfreq = paste0(dir_module, prefix_pog_pp,"_cocluster.txt")
figureout_pog_pp_coclusterfreq = paste0(dir_module, prefix_pog_pp,"_cocluster.pdf")

# read input
pog_pp = read.delim(file_pog_pp,row.names = 1,sep=','); row.names(pog_pp) = paste('pog',row.names(pog_pp),sep="")
print(paste0('calculating cocluster frequency for: ', file_pog_pp))

# calculate cocluster freq
freq = get_cocluster_freq(pog_pp,k_min, k_max, k_step, k_times)
write.table(freq, fileout_pog_pp_coclusterfreq,sep='\t',quote=F)
print(paste0('Saved to: ',fileout_pog_pp_coclusterfreq))

# plot
pdf(file = figureout_pog_pp_coclusterfreq,width = 20,height = 20)
Heatmap(as.matrix(freq),
        name='Coclustering\nfrequency',
        show_column_names = FALSE,
        clustering_distance_rows = dist_fun, clustering_distance_columns=dist_fun,
        width = ncol(freq)*unit(2.5, "mm"), height = nrow(freq)*unit(3, "mm"))
dev.off()
print(paste0('figure saved to: ',figureout_pog_pp_coclusterfreq))