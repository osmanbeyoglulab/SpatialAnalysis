library(Seurat)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

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

genelist = c()
for (i in seq(n)){ genelist = union(genelist, rownames(ov[[i]]))}
genelist = genelist %>% sort()
df.morans = MakeMoransTables(ov, genelist)

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

# cytokines
threshold = 0.15
cytokines = GeneList(obj=NULL) %>% unlist() %>% unname()
cor.max.1 = matrixStats::rowMaxs(df.morans$coef[cytokines, ] %>% replace(is.na(.), -Inf) %>% as.matrix())
names(cor.max.1) = rownames(df.morans$coef[cytokines, ])
mat1 = MergeMoransTables(mat.coef = df.morans$coef[cytokines, ] %>% .[cor.max.1>threshold, ], 
                        mat.pval = df.morans$pval[cytokines, ] %>% .[cor.max.1>threshold, ])
mat1$type = 'Cytokines/SPs'

# marker genes
marker.genes = c('PTPRC', 'CD79A', 'CD3D', 'CD4', 'CD8A', "NCAM1", 'KLRB1', 'CD14', 'CD68', 'FCGR1A', 
                 "CD40", 'FLT3', 'ITGAX', 'PECAM1', 'COL1A1', 'ACTA2', "ALDH1A3", 'EPCAM', "MUC1")
cor.max.2 = matrixStats::rowMaxs(df.morans$coef[marker.genes, ] %>% replace(is.na(.), -Inf) %>% as.matrix())
names(cor.max.2) = rownames(df.morans$coef[marker.genes, ])
mat2 = MergeMoransTables(mat.coef = df.morans$coef[marker.genes, ] %>% .[cor.max.2>threshold, ], 
                         mat.pval = df.morans$pval[marker.genes, ] %>% .[cor.max.2>threshold, ])
mat2$type = 'Markers'

mat = rbind(mat1, mat2)
mat$group[mat$sample %in% group.ids[1:5]] = 'High stromal'
mat$group[mat$sample %in% group.ids[6:11]] = 'Low stromal'
mat$sample = factor(x = mat$sample, levels = group.ids %>% rev())
mat$log10pval = -log10(mat$pval+1e-8)

p = ggplot(mat, aes(y=sample,x=gene)) +
  geom_point(aes(size=log10pval, color=coef)) + 
  scale_color_gradientn('Moran\nCorrelation',
                        colors = c(brewer.pal(9,'Greens') %>% .[1:2] %>% rev(),
                                   brewer.pal(9,'Reds')),
                        breaks = seq(-0.1,0.4,0.1), limits=c(-0.1,0.4)) +
  theme_classic() +
  facet_grid(group~type, scales = "free", space='free', switch = "y") +
  labs(size='-log10(p)') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.x = element_text(face='bold', size = rel(RELSIZE)),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 90),
        legend.title = element_text(size = rel(RELSIZE^2)),
        legend.text = element_text(size = rel(RELSIZE^2)),
        legend.key.height = unit(0.6,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("ns","*","**", '****')) 

pdf('plots_scrna_v3/st_moran_cor.pdf', width=13.5, height=6)
print(p)
dev.off()

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

df.morans.cor = rbind(df.morans$coef[marker.genes, ] %>% .[cor.max.2>threshold, ],
                      df.morans$coef[cytokines, ] %>% .[cor.max.1>threshold, ])

significance = function(p.value){
  if(p.value > 0.1){ s = 'ns' }
  else if(p.value > 0.05){ s = '.' }
  else if(p.value > 0.01){ s = '*' }
  else if(p.value > 0.001){ s = '**' }
  else {s = '***'}
  return(s)
}              

genes = rownames(df.morans$coef) %>% sort()
samples.high = sprintf('Stur_%.2d',1:5)
samples.low = sprintf('Stur_%.2d',6:11)
pvals = c()
for (gene in genes){
  x = df.morans$coef[gene, samples.high] %>% as.numeric()
  y = df.morans$coef[gene, samples.low] %>% as.numeric()
  if (sum(is.na(x))<5 & sum(is.na(y))<6){
    wt = wilcox.test(x = x, y = y, alternative = "two.sided")
    s = significance(wt$p.value)
    if (s<0.1){
      print(sprintf('%s %.4f %s', gene, wt$p.value, s))
    }
    pvals = c(pvals, s)
  }
}

# # marker genes
# marker.genes = c('PTPRC', 'CD79A', 'CD3D', 'CD4', 'CD8A', "NCAM1", 'KLRB1', 'CD14', 'CD68', 'FCGR1A', 
#                  "CD40", 'FLT3', 'ITGAX', 'PECAM1', 'COL1A1', 'ACTA2', "ALDH1A3", 'EPCAM', "MUC1")
# cor.max.2 = matrixStats::rowMaxs(df.morans$coef[marker.genes, ] %>% replace(is.na(.), -Inf) %>% as.matrix())
# names(cor.max.2) = rownames(df.morans$coef[marker.genes, ])
# mat2 = MergeMoransTables(mat.coef = df.morans$coef[marker.genes, ] %>% .[cor.max.2>threshold, ], 
#                         mat.pval = df.morans$pval[marker.genes, ] %>% .[cor.max.2>threshold, ])
# MyDotPlot(mat)
# pdf('plots_st/moran_cor_markers.pdf', width=5, height=9)
# print(MyDotPlot(mat))
# dev.off()

# # ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
#
# # Spatially variable genes
# threshold = 0.5
# cor.max = matrixStats::rowMaxs(df.morans$coef %>% replace(is.na(.), -Inf) %>% as.matrix())
# names(cor.max) = rownames(df.morans$coef)
# mat = MergeMoransTables(mat.coef = df.morans$coef[cor.max>threshold, ], 
#                         mat.pval = df.morans$pval[cor.max>threshold, ])
# MyDotPlot(mat)
# pdf('plots_st/moran_cor_svg.pdf', width=5, height=9)
# print(MyDotPlot(mat))
# dev.off()
# 
# marker.stromal = c('COL1A1', 'ACTA2')
# tab = MakeBoxTables(marker.stromal)
# tab$variable = factor(x = tab$variable, levels = c('COL1A1', 'ACTA2'))
# MyBoxPlot(tab)
# pdf('plots_st/box_ave_exp.pdf', width=8, height=4)
# print(MyBoxPlot(tab))
# dev.off()
# 
# # ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
# 
# groups = c("PR_4" = "High", "ER_6" = "High", "ER_3" = "High", "PR_5" = "High", "ER_5" = "High",
#            "ER_4" = "Low", "PR_2" = "Low", "PR_3" = "Low", "ER_1" = "Low", "PR_1" = "Low", "ER_2" = "Low")
# tab = MakeHeatmapTables('CD24', marker.genes)
# MyHeatmap(tab)
# pdf('plots_st/cor_cd24.pdf', width=5.4, height=8)
# print(MyHeatmap(tab))
# dev.off()

