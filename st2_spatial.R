library(Seurat)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)

setwd("~/Desktop/Pitt/OV")
source('aux_genelist.R')
genelist = GeneList() %>% unlist() %>% unname()
n = 11

source('aux_spatial.R')
sample.ids = c(paste0('ER_', 1:6), paste0('PR_', 1:5))
ov = readRDS(sprintf("seurat_spatial/gse189843_stur_sct%s.rds", n))
# group.ids = c(rep('ER', 6), rep('PR', 5))
# group.ids = sprintf('Stur_%.2d',1:n)
# for (i in seq(n)){ ov[[sample.ids[i]]]$Group = group.ids[i] }

ratios = list()
for (sample in sample.ids){
  data = ov[[sample]]@images[[sample]]@coordinates
  ratios[[sample]]$aspect = (max(data$imagerow)-min(data$imagerow))/(max(data$imagecol)-min(data$imagecol))
  ratios[[sample]]$font = 1e4/(max(data$imagecol)-min(data$imagecol))
}

FONTSIZE = 25
RELSIZE = 0.8

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

marker.stromal = c('COL1A1', 'ACTA2')
tab = MakeBoxTables(marker.stromal)
p = MyBoxPlot(tab) + theme(legend.position = 'bottom')
pdf('plots_scrna_v3/st_box_ave_exp.pdf', width=8, height=6)
print(p)
dev.off()

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

markers = c('COL1A1', 'ACTA2', 'EPCAM', 'CXCL12', 'CD8A', 'KLRB1')
tab.cols = as.data.frame(matrix(nrow=n, ncol=6))
rownames(tab.cols) = sample.ids
colnames(tab.cols) = markers

for (sample in sample.ids){
  obj = ov[[sample]]
  for (gene in markers %>% intersect(rownames(obj))){
    obj[[gene]] = ov[[sample]]@assays$SCT[gene, ] %>% as.numeric() %>% round(1)
    tab.cols[sample, gene] = max(obj[[gene]])
  }
  ov[[sample]] = obj
}
tab.cols[is.na(tab.cols)] = 0

colors = list()
for(gene in markers){
  n = max(tab.cols[[gene]])*10+1
  colors[[gene]] = setNames(inlmisc::GetColors(n = n, scheme = 'jet'), 
                            seq(0, max(tab.cols[[gene]]), 0.1))
}

tab.plot = as.data.frame.matrix(tab.cols) %>%
  .[, sort(colnames(.))] %>%
  mutate(ind = rownames(.)) %>% gather(variable, value, -ind)

dir.create('plots_scrna_v3/plots_st_legend')
MakeLegend = function(gene, threshold, skip){
  p = ggplot(tab.plot, aes(x=ind, y=variable)) + 
    geom_point(aes(color=value)) +
    scale_color_gradientn(gene,
                          colors=colors[[gene]] %>% as.character(),
                          breaks=seq(0,threshold,skip), limits=c(0,threshold)) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = rel(RELSIZE)),
          legend.text = element_text(size = rel(RELSIZE^2)),
          legend.key.width = unit(0.8,"cm"),
          legend.position = 'bottom',
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          text = element_text(size=FONTSIZE))
  leg = get_legend(p) %>% as_ggplot()
  pdf(sprintf('plots_scrna_v3/plots_st_legend/%s.pdf', gene), width=4, height=0.8)
  print(leg)
  dev.off()
}
gene = markers[6]
threshold = names(colors[[gene]]) %>% as.numeric() %>% max()
print(threshold)
skip = 0.5
MakeLegend(gene, threshold, skip)

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

SingleSpatialDimPlot = function(obj, feature, ratio){
  if (!(feature %in% union(rownames(obj@assays$SCT), rownames(obj@assays$Spatial)))){
    print(sprintf('%s is not found.', feature))
    p = SpatialPlot(obj, pt.size.factor = 1e-4) + ggtitle(feature) + NoLegend()
  } else {
    Idents(obj) = feature
    levels(obj) = sort(levels(obj))
    p = SpatialDimPlot(obj, pt.size.factor = ratio$font, cols = colors[[feature]][levels(obj)])
  }
  p = p + theme(aspect.ratio = ratio$aspect,
                plot.title = element_text(hjust = 0.5),
                legend.position = 'none',
                text = element_text(size=20)) + ggtitle(feature)
  return(p)
}

ncol = 6
nrow = 1
wrapper = function(markers){
  plots = list()
  for (sample in sample.ids){
    p = ggarrange(
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[1]),
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[2]),
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[3]),
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[4]),
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[5]),
      SingleSpatialDimPlot(obj = ov[[sample]], ratio = ratios[[sample]], feature = markers[6]),
      ncol = ncol, nrow = nrow
    )
    plots[[sample]] = annotate_figure(p, top = text_grob(sample2group[sample], size = FONTSIZE, face = 'bold'))
  }
  return(plots)
}

plots = wrapper(markers)
dir.create('plots_scrna_v3/plots_st')
for (sample in names(plots)){
  pdf(sprintf('plots_scrna_v3/plots_st/%s.pdf', sample2group[sample]), width=4*ncol, height=(4*ratios[[sample]]$aspect+1)*nrow)
  print(plots[[sample]])
  dev.off()
}

for (sample in c('ER_5', 'ER_2')){
  pdf(sprintf('plots_scrna_v3/plots_st_legend/%s.pdf', sample2group[sample]), width=4*ncol, height=(4*ratios[[sample]]$aspect+2)*nrow)
  print(plots[[sample]])
  dev.off()
}
