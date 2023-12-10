Load10X_Spatial_v2 = function(data.dir, image.dir = NULL, assay = 'Spatial', slice = 'slice1', filter.matrix = F){
  object = Read10X(data.dir = data.dir) %>%
    CreateSeuratObject(min.features = 400, min.cells = 3, assay = assay)
  image = Read10X_Image(image.dir = data.dir, image.name = "tissue_lowres_image.png", filter.matrix = filter.matrix)
  image = image[Cells(x = object)]
  DefaultAssay(object = image) = assay
  object[[slice]] = image
  return(object)
}

sample2group = setNames(sprintf('Stur_%.2d',1:n), 
                        c('PR_4', 'ER_6', 'ER_3', 'PR_5', 'ER_5', 'ER_4',
                          'PR_2', 'PR_3', 'ER_1', 'PR_1', 'ER_2'))

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

MakeMoransTables = function(ov, genelist){
  df.coef = data.frame(matrix(ncol = n, nrow = length(genelist)))
  rownames(df.coef) = genelist
  colnames(df.coef) = group.ids
  df.pval = df.coef
  for (sample in names(ov)){
    df = ov[[sample]]@assays[["SCT"]]@meta.features
    genes = rownames(df)[!is.na(df$MoransI_observed)]
    df.coef[genes, sample2group[sample]] = df[genes, 'MoransI_observed']
    df.pval[genes, sample2group[sample]] = df[genes, 'MoransI_p.value']
  }
  return(list(coef=df.coef, pval=df.pval))
}

MergeMoransTables = function(mat.coef, mat.pval){
  mat.coef = mat.coef %>% mutate(ind = rownames(.)) %>% gather(variable, value, -ind)
  colnames(mat.coef) = c('gene', 'sample', 'coef')
  mat.pval = mat.pval %>% mutate(ind = rownames(.)) %>% gather(variable, value, -ind)
  colnames(mat.pval) = c('gene', 'sample', 'pval')
  mat = mat.coef %>% merge(mat.pval, by = c('gene', 'sample'))
  return(mat)
}

# ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

MakeBoxTables = function(markers){
  n = length(markers)
  cell_list = c()
  sample_list = c()
  cell_sample_list = c()
  for (sample in sample.ids){
    cells = colnames(ov[[sample]])
    samples = ov[[sample]]$Sample
    cell_list = c(cell_list, cells)
    sample_list = c(sample_list, samples)
    cell_sample_list = c(cell_sample_list, sprintf("%s (%s)", cells, samples))
  }
  ncells = length(cell_list)
  tab.temp = as.data.frame(matrix(nrow=ncells, ncol=n))
  rownames(tab.temp) = cell_sample_list
  colnames(tab.temp) = markers
  tab.ave = as.data.frame(matrix(nrow=length(sample.ids), ncol=n))
  rownames(tab.ave) = sample.ids
  colnames(tab.ave) = markers
  idx = 0
  for (sample in sample.ids){
    ncells = length(colnames(ov[[sample]]))
    for (gene in markers){
      if (gene %in% rownames(ov[[sample]][["SCT"]]@data)){
        tab.temp[(idx+1):(idx+ncells), gene] = ov[[sample]][["SCT"]]@data[gene, ]
        tab.ave[sample, gene] = mean(ov[[sample]][["SCT"]]@data[gene, ])
      }
    }
    idx = idx + ncells
  }
  tab = tab.temp%>%
    .[, sort(colnames(.))] %>%
    mutate(ind = rownames(.)) %>% gather(variable, value, -ind)
  tab$sample = str_sub(tab$ind, start = -5, end = -2)
  tab$sample = factor(
    tab$sample,
    levels = tab.ave[order(tab.ave$COL1A1),] %>% rownames() %>% rev())
  tab$group = sample2group[tab$sample] %>% unname()
  tab$variable = factor(x = tab$variable, levels = markers)
  return(tab)
}

MyBoxPlot = function(tab){
  colors = c(inlmisc::GetColors(21,'BuRd') %>% .[3:7],
             inlmisc::GetColors(30, scheme = 'bpy') %>% .[24:29])
  p = ggplot(tab, aes(x = variable, y = value, fill = group)) + 
    geom_boxplot() + 
    scale_fill_manual(values = colors) +
    facet_grid(~variable, scales = "free", space='free') +
    theme_classic() + ylab('Expression') + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size=16, face='bold'),
      legend.title = element_blank(),
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      # axis.text.x = element_text(angle = 90, face='bold'),
      text = element_text(size=FONTSIZE))
  return(p)
}

MakeHeatmapTables = function(target_gene, markers){
  tab = as.data.frame(matrix(nrow=length(sample.ids), ncol=length(markers)))
  rownames(tab) = sample.ids
  colnames(tab) = markers
  for (sample in sample.ids){
    exp = ov[[sample]][['Spatial']]@data %>% as.data.frame() %>% t()
    genes = markers %>% intersect(colnames(exp))
    s = cor(x=exp[, target_gene], y=exp[, genes], method = 'spearman')
    tab[sample, genes] = s
  }
  return(tab)
}

MyHeatmap = function(tab, threshold=0.4){
  tab.filtered = replace(tab, is.na(tab), 0) %>% .[, (.>threshold) %>% apply(2, any)]
  p = Heatmap(t(tab.filtered), na_col = 'white',
              col = colorRamp2(c(-0.2, 0, 1), c("#377EB8", "white", "#E41A1C")),
              cluster_rows = T, cluster_columns = T,
              show_column_names = T, column_names_side = c("top"), column_names_rot = 45, column_names_centered = T,
              name = 'Correlation \nwith CD24',
              top_annotation = HeatmapAnnotation(
                Group = groups[rownames(tab.temp)],
                col = list(Group = c("High" = "#44BB99", "Low" = "#FFAABB"))))
  return(p)
}
