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
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="metadata file"),
  make_option(c("-c", "--metadata_colname"), type="character", default=NULL, help="column name for metadata in the metadata file"),

  make_option(c( "--dir_pog"), type="character", default=NULL),
  make_option(c( "--dir_module"), type="character", default=NULL),

  make_option(c( "--prefix_protein"), type="character", default=NULL, help="prefix"),
  make_option(c( "--prefix_pog"), type="character", default=NULL, help="prefix"),
  make_option(c( "--prefix_pog_pp"), type="character", default=NULL, help="prefix for pog phylogenetic profiles"),
  make_option(c( "--prefix_module"), type="character", default=NULL, help="prefix"),

  # PARAMS
  make_option(c("-n", "--module_n"), type="character", default=NULL, help="# of modules")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_metadata = opt$metadata
metadata_colname = opt$metadata_colname

dir_pog = opt$dir_pog
dir_pog_pp = opt$dir_pog_pp
dir_module = opt$dir_module

prefix_protein = opt$`prefix_protein`
prefix_pog = opt$`prefix_pog`
prefix_pog_pp = opt$`prefix_pog_pp`
prefix_module = opt$`prefix_module`

MODULE_N = opt$module_n
# POG_PROTEIN_N_MIN=opt$min_protein_N
POG_PROTEIN_N_MIN=2

dist_fun = function(m) {
  return( as.dist(1-as.matrix(m)))
}

# [input]
file_result = paste0(dir_pog, prefix_protein, '.csv')
file_protein_annot = paste0(dir_pog, prefix_protein, '.annot')
file_pog_annot = paste0(dir_pog, prefix_pog, '.annot')
file_pog_pp = paste0(dir_pog, prefix_pog_pp, '.pp')
file_coclusterfreq = paste0(dir_module, prefix_pog_pp, '_cocluster.txt')

# [output]
# ----
figureout_pog_pp = paste0(dir_module, prefix_pog,"_pp.pdf")
figureout_pog_cocluster = paste0(dir_module, prefix_module, '.pdf')
figureout_module = paste0(dir_module, prefix_module, '_cuttree.pdf')
fileout_module = paste0(dir_module, prefix_module, '.txt')
# ----


# Import
# ----
metadata = read.delim(file_metadata)
result <- read.csv(file_result,row.names=1)
freq = read.table(file_coclusterfreq, sep='\t',check.names = F); row.names(freq) = colnames(freq) = gsub('family','pog',row.names(freq))
pp = read.csv(file_pog_pp,row.names = 1);
if (sum(colnames(pp) %in% metadata$taxon_oid) ==0) {
  metadata$taxon_oid = paste0('genome',metadata$taxon_oid)# jgi specific 
}
pp = pp[,metadata$taxon_oid]; rownames(pp) = paste0('pog',rownames(pp))
annot_protein = read.delim(file_protein_annot); annot_protein$pog = paste0('pog',annot_protein$pog)
annot = read.delim(file_pog_annot,row.names = 1);rownames(annot) = paste0('pog',rownames(annot));annot = annot[row.names(freq),]

print(head(annot_protein))
print(head(annot))
print(head(result))
print(head(freq))
print(head(pp[1:5]))


print(dim(annot))
print(dim(pp))
print(dim(result))
print(dim(freq))
print(dim(annot_protein))

annot = annot[annot$N>POG_PROTEIN_N_MIN,]
pp = pp[row.names(annot),]
freq = freq[row.names(annot),row.names(annot)]

print(dim(annot))
print(dim(pp))
print(dim(result))
print(dim(freq))
print(dim(annot_protein))

# ----

# Annotation
# ----
row_annot = unlist(lapply(as.character(annot$keywords),function(x) strsplit(x,";")[[1]][1]))
row_annot2 = gsub(" [(]n=[0-9]*[)]","",row_annot)
row_annot[row_annot2==''] = ''
annot_mi_z = tapply(result[as.character(annot_protein$protein),'mi_z'],
                    annot_protein$pog,
                    mean)[row.names(freq)]
# print(annot_mi_z)
row_annot = data.frame('pog' = row.names(freq),
                       'N' = annot$N,
                       'keyword'=row_annot ,
                       # 'keywrod' = annot$keywords,
                       'mi' = apply(pp,1,function(x) mutinformation(x,metadata[,metadata_colname])),
                       'mi_z' = annot_mi_z[row.names(freq)],
                       'spearman_r' = apply(pp,1,function(x) cor(x,metadata[,metadata_colname],method = 'spearman')))
# print(row_annot)
# print(result[as.character(annot_protein$protein),'mi_z'])
# print(row_annot$mi_z)
# print(row_annot$spearman_r)
# ----

# -- plot pog phylogenetic profiles --
# ----
row_labels = paste0(row_annot$pog," (n=",row_annot$N,")",": ",row_annot$keyword)
pdf(file = figureout_pog_pp,width = 18, height = 11)
Heatmap(pp,
        row_labels = row_labels, #row_labels = row_annot$keyword,,
        row_names_gp = gpar(fontsize=8),row_names_side = 'right',
        width = ncol(pp)*unit(0.1, "mm"), height = nrow(pp)*unit(3, "mm"),
        show_heatmap_legend = F,
        col = colorRampPalette(c("lightgray", "blue"))(2) ,
        column_split= metadata[,metadata_colname],
        show_column_names = F, show_row_dend = F,
        left_annotation = rowAnnotation('MI'=anno_barplot(row_annot$mi,
                                                          width = unit(2, "cm"),
                                                          axis_param =c(side='bottom')),
                                        'MI zscore'=anno_barplot(row_annot$mi_z,
                                                                 width = unit(2, "cm"),
                                                                 axis_param =c(side='bottom')),
                                        'Spearman r'=anno_barplot(row_annot$spearman_r,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                        annotation_name_side='top',
                                        annotation_name_rot=0,
                                        gap = unit(5, "points")
                                        ))
dev.off()
print(paste0('figure saved to: ',figureout_pog_pp))
# ----

# -- plot POG cocluster heatmap
# ----
plot_pog_cocluster = function(figureout_pog_cocluster)

pdf(file = figureout_pog_cocluster, width = 20,height = 20)
mat = as.matrix(freq)
ha=Heatmap(mat,name='Coclustering\nfrequency',

           show_column_names = T,column_names_gp = gpar(fontsize=8),

           clustering_distance_rows = dist_fun, clustering_distance_columns=dist_fun,
           clustering_method_rows = "average",clustering_method_columns= "average",

           width = ncol(mat)*unit(3.5, "mm"), height = nrow(mat)*unit(3.5, "mm"),

           row_labels = row_labels, #row_labels = row_annot$keyword,
           show_row_dend = F, row_names_gp = gpar(fontsize=8), row_names_side = 'right',
           row_names_max_width = max_text_width(row_annot$keyword, gp = gpar(fontsize = 8)),

           left_annotation = rowAnnotation('MI'=anno_barplot(row_annot$mi,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           'MI zscore'=anno_barplot(row_annot$mi_z,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           'Spearman r'=anno_barplot(row_annot$spearman_r,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           annotation_name_side='top',annotation_name_rot=0,gap = unit(5, "points"))
)

ha
dev.off()
print(paste0('figure saved to: ',figureout_pog_cocluster))
# ----

# -- plot module heatmap
# ----
pdf(file = figureout_module, width = 20,height = 20)
mat = as.matrix(freq)
ha=Heatmap(mat,name='Coclustering\nfrequency',

           show_column_names = T,column_names_gp = gpar(fontsize=8),

           clustering_distance_rows = dist_fun, clustering_distance_columns=dist_fun,
           clustering_method_rows = "average",clustering_method_columns= "average",
           #         cluster_rows = dend,cluster_columns=dend,

           width = ncol(mat)*unit(3.5, "mm"), height = nrow(mat)*unit(3.5, "mm"),
           row_split = as.numeric(MODULE_N),  row_title = "module %s",
           row_title_rot = 0, row_title_side = 'left',


           row_labels = row_labels, #row_labels = row_annot$keyword,
           show_row_dend = F, row_names_gp = gpar(fontsize=8), row_names_side = 'right',
           row_names_max_width = max_text_width(row_annot$keyword, gp = gpar(fontsize = 8)),

           left_annotation = rowAnnotation('MI'=anno_barplot(row_annot$mi,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           'MI zscore'=anno_barplot(row_annot$mi_z,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           'Spearman r'=anno_barplot(row_annot$spearman_r,width = unit(2, "cm"),axis_param =c(side='bottom')),
                                           annotation_name_side='top',annotation_name_rot=0,gap = unit(5, "points"))
           )

ha
dev.off()
print(paste0('figure saved to: ',figureout_module))
# ----

# module-pog relationship
# ----
module_pog = row_order(ha)
module_pog = cbind(rep(paste0('module',c(1:length(module_pog))),unlist(lapply(module_pog,length))),
                   unlist(lapply(module_pog,function(x) row_annot$pog[x])),
                   unlist(lapply(module_pog,function(x) row_annot$keyword[x])))
write.table(module_pog,fileout_module,row.names = F, quote=F,sep='\t')
# ----

dim(pp)
length(row_labels)
