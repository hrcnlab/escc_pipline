# escc_gsva
setwd("~/R/ESCC/escc_seq")
load("escc_seq.Rdata")
# part1 gsva

# load the packages 
pacman::p_load(clusterProfiler,ggplot2,GSVA,GSEABase,BiocParallel)

ac<-colData(escc_seq)
escc_counts<-assay(escc_seq,"counts")
escc_fpkm<-assay(escc_seq,"fpkm")
escc_fpkm<-escc_fpkm[,ac_sub$sample_id]
boxplot(escc_fpkm)

 
load("gmt_pathways.gmt"） # loading the gmt list from collections
gmt_list_raw_sp

gmt_types<- read.csv('gmt_types.cvs',header=TRUE)


head(define_gmt_types);head(nrDEGs)
nrDEGs_hsa_class<-inner_join(nrDEGs,define_gmt_types,
                             by=c("gslist"="pathways"))

ssgsea_mats <- 
  GSVA::gsva(expr=as.matrix(escc_zfpkm), 
             method='ssgsea',
             kcdf="Gaussian", # Poisson
             gset.idx.list= mylsits, 
             max.sz=500,
             min.sz = 3, # min.sz=5, max.sz=500
             verbose=TRUE, 
             parallel.sz=4L,
             BPPARAM=SerialParam()
  )



ssgsea_mats<-ssgsea_scores_filter 

head(ssgsea_mats)

# DE pathways using limma

library(limma)
group_list<-droplevels(ac_sub$his)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=ssgsea_mats
rownames(design)=colnames(exprSet)
design
colnames(design)
contrast.matrix<-makeContrasts(contrasts=c("low_grade-Normal",
                                           "high_grade-Normal",
                                           "cancer-Normal"
                                           ),
                               levels = design)
contrast.matrix 

fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  ## default no trend 

contranst_list<-c("low_grade-Normal",
  "high_grade-Normal",
  "cancer-Normal"
)
library(rgl)
library(tidyverse)

nrDEGs<-purrr::map_dfr(contranst_list,function(x){
  topTable(fit2, coef=x, n=Inf) %>%
    rownames_to_column(var="gslist") %>%
    mutate(groups=x)
})

head(define_gmt_types);head(nrDEGs)
nrDEGs_hsa_class<-inner_join(nrDEGs,gmt_types,
                             by=c("gslist"="pathways"))
							 
  nrDEGs_hsa_class %>% 
  mutate(groups_reg=if_else(logFC>0,"up","down")) %>%
  filter(adj.P.Val<0.01 & abs(logFC)>0, groups_reg=="up" ) %>%
  group_by(groups,source) %>%
  arrange(adj.P.Val, desc(logFC), .by_group = TRUE) %>%
  dplyr::slice(1:10) %>% # 
  pull(gslist) %>%
  unique() -> plot_glists

# part2 visualization
# Heatmap visualization

library(ComplexHeatmap)
library(circlize)
library(paletteer) 
library(ggsci)
library(scales)

set.seed(372)
groups_sample<-droplevels(ac_sub$his)

dend1 = cluster_between_groups(ssgsea_mats, groups_sample)
groups.cor<-mypal[1:nlevels(groups_sample)] # define the color
groups.cor<-c("#FBC1B8","#87C6BB","#FDD073","#ADA3C6")
names(groups.cor)<-levels(factor(groups_sample))


gmt_types_sub<-define_gmt_types[match(plot_glists,
                                             gmt_types$pathways),]

gmt_types_sub<-
  gmt_types_sub |> mutate(new_types=case_when(
  source == "paper_cell" ~ "Immune_related",
  source == "Immnue17" ~ "Immune_related",
  source == "Immune_combine" ~ "Immune_related",
  .default = as.character(source)
))

row_split <- factor(gmt_types_sub$new_types) %>% droplevels()
unique(row_split)
row_split<-factor(row_split,levels = c(
  "GO","kegg_all","Abnormal_metab","reactome_selected","Immune_related"
))
left_annotation = left_anno
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.85)(9)
mypal<-mypal[-5]
library("scales")
show_col(mypal)
type.cor<-rev( mypal[1:nlevels(row_split)]) # Adjust color
names(type.cor)<-levels(row_split)

left_anno = rowAnnotation(
  cluster = anno_block(gp = gpar(fill = type.cor,
                                 col = NA),
                       labels = names(type.cor),
                       labels_gp = gpar(col = "white",
                                        fontsize = 8)))
										
ht<-
  ComplexHeatmap::Heatmap(sce.all[plot_glists,],
        # col=col_fun, # define your ht color
        # rect_gp = gpar(col = "white", lwd = 0.5), # add the border
        cluster_columns = dend1, 
        border_gp=FALSE,
        row_names_max_width = unit(10, "cm"), 
        column_split = 4,
        row_split =row_split,
        left_annotation = left_anno,
        show_row_dend = T, show_column_dend = F,
        row_dend_width = unit(0.2, "cm"), 
        cluster_row_slices = F,cluster_column_slices = F,
        show_heatmap_legend = T,
        border = F,
        row_title = NULL,column_title=NULL,
        # row_title = "cluster_between_groups",
        # col_title = NULL,
        show_column_names = F,
        row_names_gp = grid::gpar(fontsize = 10,
                                  col=type.cor,
                                  fontface ="bold.italic"),
        top_annotation = HeatmapAnnotation(

        foo=anno_block(gp = gpar(fill = groups.cor, 
                                  col = NA),
        labels_gp = gpar(col = "white", fontsize = 18) )
        ))

ht_lg <-Legend(col_fun = colorRamp2(seq(-4,4),
                                   rev(colfunc(9)) ), title = "fooht")
 

lgd <- Legend(labels = names(groups.cor),
             title = "foo",
             grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
             legend_gp = gpar(fill = groups.cor,
                              col = "Black", fontsize = 14),
             title_position = "leftcenter") #leftcenter-rot
draw(lgd)
dev.off()
draw(ht, heatmap_legend_side = "left", 
     annotation_legend_list = list(lgd),
     merge_legend = TRUE)

# patr 3 
# SCORPIUS analysis

library(SCORPIUS)
ac<-colData(escc_seq)
escc_fpkm<-assay(escc_seq,"fpkm")
escc_fpkm<-escc_fpkm[,ac_sub$sample_id]
boxplot(escc_fpkm)

expression <- escc_fpkm
group_name <- ac$his
table(group_name)
dim(expression)
expression[1:4,1:4]

 table(group_name)
 dim(expression)
expression[1:4,1:4]


space <- reduce_dimensionality(expression, "spearman")
draw_trajectory_plot(space, group_name, contour = TRUE)
space <- reduce_dimensionality(expression, "spearman")
draw_trajectory_plot(space, group_name, contour = TRUE)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, group_name, traj$path, contour = TRUE)

# warning: setting num_permutations to 10 requires a long time (~30min) to run!
# set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
gimp <- gene_importances(
  expression, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expression[,gene_sel])

# Draw a time series heatmap
time <- traj$time
draw_trajectory_heatmap(expr_sel, time)

## Also show the progression groupings
draw_trajectory_heatmap(expr_sel, time, 
                        progression_group=group_name)

## define your colors
table(group_name)
draw_trajectory_heatmap(
  expr_sel, time, progression_group=group_name,
  progression_group_palette = setNames(RColorBrewer::brewer.pal(3, "Set2"),
                                       levels(group_name))
)

## extact models
modules <- extract_modules(scale_quantile(expr_sel))
draw_trajectory_heatmap(expr_sel, time,
                        progression_group=group_name, modules=modules)

# checke you genes
genes_check_list<-"your gene lists"

library(ggplot2)
draw_trajectory_plot(space, expression[,"mygene"], 
                     progression_group_palette = 
                       setNames(RColorBrewer::brewer.pal(3, "Set1"),
                                                  levels(group_name)),
                     point_size =1.2,
                                      traj$path) + 

  labs(colour = "your Gene\nexpression")
