library(Seurat)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Pitt/OV")
n = 11
source('aux_spatial.R')
source('aux_genelist.R')

file.dirs = paste0('seurat_raw/GSE189843_Stur/II214', seq(72,83))
file.dirs[8:12] = file.dirs[c(9:12,8)]
sample.ids = c(paste0('ER_', 1:6), paste0('PR_', 1:5))
ov = readRDS(sprintf("seurat_spatial/gse189843_stur_sct%s.rds", n))
group.ids = sprintf('Stur_%.2d',1:n)
# group.ids = c(rep('ER', 6), rep('PR', 5))

FONTSIZE = 25
RELSIZE = 0.8

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

group2color = setNames(
  c(inlmisc::GetColors(21,'BuRd') %>% .[3:7], inlmisc::GetColors(30, scheme = 'bpy') %>% .[24:29]),
  sprintf('Stur_%.2d',1:n))
genes = c('CXCL12', 'COL1A1', 'ACTA2')

df.cor = data.frame(matrix(ncol = n, nrow = 2))
colnames(df.cor) = group.ids
rownames(df.cor) = c('COL1A1', 'ACTA2')

df.cor.cxcl12 = df.cor
for (sample in sample.ids){
  obj = ov[[sample]]
  if ('CXCL12' %in% rownames(obj@assays$Spatial)){
  df.cor.cxcl12['COL1A1', sample2group[[sample]]] = cor(
    x=obj@assays$Spatial@data['CXCL12', ], y=obj@assays$Spatial@data['COL1A1', ], method = 'spearman')
  df.cor.cxcl12['ACTA2', sample2group[[sample]]] = cor(
    x=obj@assays$Spatial@data['CXCL12', ], y=obj@assays$Spatial@data['ACTA2', ], method = 'spearman')
  } 
}

df2colum = function(df, target){
  df.column = df %>% mutate(ind = rownames(.)) %>% gather(variable, value, -ind)
  colnames(df.column) = c('gene', 'sample', 'corr')
  df.column$group[df.column$sample %in% group.ids[1:5]] = 'High stromal'
  df.column$group[df.column$sample %in% group.ids[6:11]] = 'Low stromal'
  df.column$sample = factor(x = df.column$sample, levels = group.ids %>% rev())
  df.column$target = target
  return(df.column)
}

df.cor.column = df2colum(df.cor.cxcl12, 'CXCL12')
df.cor.column[is.na(df.cor.column)] = -1

# p = ggplot(df.cor.column, aes(y=sample,x=gene,size = 3)) +
#   geom_point(aes(size = corr, color=corr)) + 
#   scale_color_gradientn('Correlation',
#                         colors = brewer.pal(9,'YlOrRd'), na.value = 'white',
#                         breaks = seq(0,0.75,0.25), limits=c(0,0.75)) +
#   theme_classic() +
#   facet_grid(group~target, scales = "free", space='free', switch = "y") +
#   theme(panel.grid = element_blank(),
#         strip.placement = "outside",
#         strip.background = element_rect(colour="grey", fill="white"),
#         strip.text.x = element_text(face='bold', size = rel(RELSIZE)),
#         strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 90),
#         legend.title = element_text(size = rel(RELSIZE^2)),
#         legend.text = element_text(size = rel(RELSIZE^2)),
#         legend.key.height = unit(0.6,"cm"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         text = element_text(size=FONTSIZE)) 
# 
# pdf('plots_scrna_v3/st_pairwise_cor.pdf', width=7.5, height=6)
# print(p)
# dev.off()

ha = HeatmapAnnotation(Group = c(rep('High stromal', 5), rep('Low stromal', 6)),
                   col = list(Group = c('High stromal' = "#437DBF", 'Low stromal' = '#F7CB45')))

p = Heatmap(df.cor.cxcl12 %>% as.matrix(), 
        col = colorRamp2(seq(0,0.8,0.1), brewer.pal(9,'Reds')),
        cluster_rows = F, cluster_columns = F,
        name = 'Correlation\nwith\nCXCL12',
        # row_names_gp = gpar(fontsize = FONTSIZE*RELSIZE^3),
        # column_names_gp = gpar(fontsize = FONTSIZE*RELSIZE^3),
        # heatmap_legend_param = gpar(fontsize = 10),
        top_annotation = ha)
pdf('plots_scrna_v3/st_pairwise_cor.pdf', width=8, height=2)
print(p)
dev.off()

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

group2color = setNames(
  c(inlmisc::GetColors(21,'BuRd') %>% .[3:7], inlmisc::GetColors(30, scheme = 'bpy') %>% .[24:29]),
  sprintf('Stur_%.2d',1:n))
genes = c('CXCL12', 'COL1A1', 'ACTA2')
df.scatter = data.frame(matrix(ncol = 0, nrow = length(genes)+2))
rownames(df.scatter) = c(genes, 'color')
for (sample in sample.ids){
  obj = ov[[sample]]
  df.exp = data.frame(matrix(ncol = dim(obj)[2], nrow = length(genes)+2))
  colnames(df.exp) = colnames(obj)
  rownames(df.exp) = c(genes, 'sample', 'color')
  df.exp[genes, ] = obj@assays$Spatial@data[genes, ]
  df.exp['sample',] = sample2group[[sample]]
  df.exp['color',] = group2color[[sample2group[[sample]]]]
  df.scatter = cbind(df.scatter, df.exp)
}

df.plot = df.scatter %>% t() %>% as.data.frame()
df.plot$CXCL12 = df.plot$CXCL12  %>% as.numeric()
df.plot$COL1A1 = df.plot$COL1A1 %>% as.numeric()
df.plot$ACTA2 = df.plot$ACTA2 %>% as.numeric()

s = cor(df.plot$CXCL12 %>% as.numeric(), df.plot$COL1A1  %>% as.numeric())
ggplot(df.plot, aes(x = CXCL12, y = COL1A1, color = sample)) +
  geom_point(size = 1) + theme_bw() +
  scale_color_manual(name = NULL,
                     values = group2color) +
  xlab(sprintf('%s expression', 'CXCL12')) + ylab(sprintf('%s expression', 'COL1A1')) + 
  ggtitle(sprintf('corr = %.4f', s)) +
  theme(plot.margin = margin(t=10, r=20, b=10, l=10),
        # axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE), face = 'bold'),
        legend.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.text = element_text(size = rel(RELSIZE)),
        axis.title = element_text(size = rel(RELSIZE), face = 'bold'),
        text = element_text(size=FONTSIZE)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) + guides(fill=guide_legend(ncol=2))
