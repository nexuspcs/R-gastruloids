packages <- c(
  #Data manipulation and analytics
  "SummarizedExperiment","dplyr","scater","scran",
  #Visualization
  "RColorBrewer","dendsort","ggplot2","ComplexHeatmap","gridExtra","circlize")
lapply(packages, require, character.only = TRUE)

# I. DATA PROCESSING ----------

## Input Persues file   ----------

#Processed Persues files in a MQ run
gasBulk.P <- read.delim("BulkData/gasBulk Persues.txt")
gasBulk.P <- gasBulk.P[-c(1:2), ]

#Sox17 - Mt - Bra
RePro <- c("Q61473","P02802","P20293")
Rege <- c("Sox17", "Mt1", "Tbxt", "T", "Mt")


#Cleaning
#Remove >2nd leading proteins
gasBulk.P$Majority.protein.IDs <- sub(';.*','',gasBulk.P$Majority.protein.IDs)
#Remove >2nd leading genes
gasBulk.P$Gene.names <- sub(';.*','',gasBulk.P$Gene.names)

#Find proteins that have no correspoding gene
blankgene <- gasBulk.P[which(gasBulk.P[,"Gene.names"] == ""),c("Majority.protein.IDs","Gene.names")]

#Check with uniprot
blankgene

#Manual add genes
gasBulk.P[which(gasBulk.P[,"Gene.names"] == ""),c("Gene.names")] <- c("Ldah","0610012G03Rik","Bmt2","2210016L21Rik","Pnkd","Plekha7","Zfp654",
                                                              "2310057M21Rik","Gm11639","Gcna","4931406C07Rik","Bclaf3","Arxes2","Mfap1b",
                                                              "Ndufb1","Uncharacterized protein C1orf112 homolog","Dipk2a","Hnrnpll","UPF0711 protein C18orf21 homolog","Carnmt1","UPF0489 protein C5orf22 homolog",
                                                              "Timm29","Maip1","UPF0600 protein C5orf51 homolog","Hpf1","UPF0415 protein C7orf25 homolog","UPF0729 protein C18orf32 homolog","Ashwin",
                                                              "Uncharacterized protein C4orf3 homolog","Uncharacterized protein C7orf50 homolog","Trir")

## 1. Reconstructing the dataset to Create singleCellExperiment  object ----------

#singleCellExperiment object (SCE) help to perform PCA with scran and Marker detection in later steps

#Transfer into numeric matrix
gasBulk.sce <- apply(gasBulk.P[1:12],2, as.numeric)
#Change rownames of matrix
rownames(gasBulk.sce) <- gasBulk.P$Gene.names
#colData for SCE object
Population <- c("Mt1-BFP+","Mt1-BFP+","Mt1-BFP+",
                "Bra-GFP+","Bra-GFP+","Bra-GFP+",
                "Triple-Neg","Triple-Neg","Triple-Neg",
                "Sox17-RFP+","Sox17-RFP+","Sox17-RFP+")

#Create SCE object
#Scale matrix is convinient for Heatmap plot
gasBulk.sce <- SingleCellExperiment(
  assays = list(assay = gasBulk.sce, #log value matrix without scaling
                scaleassay = t(scale(t(gasBulk.sce),center = TRUE, scale = TRUE))),#log value matrix withscaling
  colData=DataFrame(Population=Population))

## 2. Quality control at Protein level  -----

#Already done with Perseus

## 3. Heatmap for significant expressed proteins  -----------

#Assessing the significant proteins that passed ANOVA test
#Settings: Multiway anova in Perseus with s0=2 and FDR=0.05
#There were 991 significant proteins that passed the ANOVA test which note as "+" in ANOVA.Significant columns

gasBulk.aov <- gasBulk.P[which(gasBulk.P[,"ANOVA.Significant"] == "+"),]

#Transfer into numeric matrix
gasBulk.aov.hm <- as.matrix(sapply(gasBulk.aov[1:12], as.numeric))
#Scale matrix
gasBulk.aov.hm = t(scale(t(gasBulk.aov.hm)))

#Colour
#For heatmap
colour.hm <-  c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
#For Channels
colour.Chan <- c("dodgerblue","#E03E3E","#53BB79","#EF8228","#937264","#3043A2","#C25D7B")

#Manually create legend of Cells Population
lgChan = Legend(labels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg"), 
                legend_gp = gpar(fill = colour.Chan), 
                title = "Population",direction = "horizontal")

### Hierarchical clustering -----

# Clustering with different distance methods:
#dist(1-cor(t(gasBulk.aov.hm))) = euclidean
#"pearson"
#spearman
#kendall
dist.method <- "euclidean"
gasBulk.aov.hm.hc <- Heatmap(gasBulk.aov.hm, name = "Scaled Log Abundance",
                             clustering_distance_rows = dist.method,
                             clustering_method_rows = "ward.D",
                             cluster_columns = dendsort(hclust(dist(t(gasBulk.aov.hm)))),
                             #row_split = 4,
                             #column_split = 4,
                             col = colour.hm,
                             heatmap_legend_param = list(direction = "horizontal"),
                             top_annotation = HeatmapAnnotation(Population = as.character(Population),
                                                                border = TRUE,gap = unit(2, "points"),
                                                                col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                          "Sox17-RFP+" = "#E03E3E",
                                                                                          "Bra-GFP+" = "#53BB79",
                                                                                          "Triple-Neg" = "#EF8228")),
                                                                show_legend = FALSE),
                             heatmap_width = unit(8, "cm"), heatmap_height = unit(18, "cm"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 30),
                             #column_title = paste0("All ", nrow(gasBulk.aov) ," significant proteins - ",dist.method)
)
gasBulk.aov.hm.hc

png(paste0("Figure/I.3. Heatmap on All ", nrow(gasBulk.aov) ," significant proteins - ",dist.method," .png"),width = 5.25,height = 8.25,units="in",res = 320)
draw(gasBulk.aov.hm.hc, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
draw(lgChan, x = unit(0.93, "npc"), y = unit(0.665, "npc"), just = c("right", "top"))
dev.off()

### k-means clustering -----

n <- 4
kclus <- kmeans(gasBulk.aov.hm, n)
#Heatmap of 991 significant hits
split <- paste0("Cluster\n", kclus$cluster)
gasBulk.aov.hm.km <- Heatmap(gasBulk.aov.hm, name = "Scale Log \nAbundance",
                             row_dend_reorder = TRUE,
                             cluster_row_slices = FALSE,
                             split=split, 
                             clustering_method_rows = "ward.D",
                             cluster_columns = dendsort(hclust(dist(t(gasBulk.aov.hm)))),
                             #column_split = 4,
                             col = colour.hm,
                             top_annotation = HeatmapAnnotation(Population = as.character(Population),
                                                                border = TRUE,gap = unit(2, "points"),
                                                                col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                     "Sox17-RFP+" = "#E03E3E",
                                                                                     "Bra-GFP+" = "#53BB79",
                                                                                     "Triple-Neg" = "#EF8228")),
                                                                show_legend = FALSE),
                             heatmap_width = unit(8, "cm"), heatmap_height = unit(18, "cm"),
                             column_names_gp = grid::gpar(fontsize = 10),
                             row_names_gp = grid::gpar(fontsize = 30),
                             column_title = paste0("All ", nrow(gasBulk.aov) ," significant proteins - k-means"))

gasBulk.aov.hm.km

png(paste0("Figure/I.3. Heatmap on All ", nrow(gasBulk.aov) ," significant proteins - k-means .png"),width = 5.25,height = 8.25,units="in",res = 320)
gasBulk.aov.hm.km
draw(lgChan, x = unit(0.93, "npc"), y = unit(0.665, "npc"), just = c("right", "top"))
dev.off()

# II. DIMENSIONALITY REDUCTION AND VISUALIZATION ----------

#colour PCA
colour.PCA <- c("#53BB79","dodgerblue","#E03E3E","#EF8228")

## 1. Principal component analysis  -------

# Exploring 1000 most variable protein: ntop = 1000
set.seed(350)
PCAtop <- 1000
gasBulk.sce <- runPCA(gasBulk.sce,
                      ncomponents = 5,
                      ntop = PCAtop,
                      scale = FALSE,
                      exprs_values = "assay",
                      name = "PCA")

PCAplot <- plotReducedDim(gasBulk.sce,
                          dimred = "PCA",
                          colour_by = "Population",
                          point_alpha = 0.8) +
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gasBulk.sce), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gasBulk.sce), "percentVar")[2],1), "%)"))
PCAplot

ggsave(paste0("Figure/II.1.- PCA for top ",PCAtop," PCA determinants - PC.png"), 
       plot=PCAplot, width = 3, height = 2, dpi = 320)

## 2. Percentage of variance explained  -------------

pv <- attr(reducedDim(gasBulk.sce), "percentVar")
elbow <- PCAtools::findElbowPoint(pv)
pvplot <- data.frame(c(seq(1,length(pv))),pv)
pvPlot <- pvplot %>% ggplot(aes(x=pvplot[,1],y=pvplot[,2])) +
  geom_line(color="blue") +
  labs(title="Plot % variance explained",
       x ="PC",
       y = "% variance explained")+
  geom_vline(xintercept = elbow,
             linetype="dashed",
             color = "red", size=1)+
  theme_classic()
pvPlot

ggsave("Figure/II.2. pvPlot.png", 
       plot=pvPlot, width = 3,height = 3, dpi = 320)

#Top 3PCs explain for %

pvplot$pv[1:3] %>% sum
pvplot$pv[1:elbow] %>% sum

#Take only 3 PC
reducedDim(gasBulk.sce, "PCA.elbow") <- reducedDim(gasBulk.sce)[,1:elbow]
reducedDimNames(gasBulk.sce)

#Check first 3 PC correlation
PCAtoplot <- plotPCA(gasBulk.sce, 
                     ncomponents = elbow, 
                     colour_by = "Population")+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")
  #+  theme(legend.position="none")
PCAtoplot
ggsave("Figure/II.2.PCAtoplot.png", 
       plot=PCAtoplot, width = 4,height = 3, dpi = 320)

## 3. Uniform Manifold Approximation and Projection  ------------

set.seed(500)
gasBulk.sce <- runUMAP(gasBulk.sce,
                 ncomponents = elbow,
                 exprs_values = "assay",
                 n_neighbors = 4,
                 dimred = "PCA.elbow",
                 name = "UMAP")
UMAPplot <- plotReducedDim(gasBulk.sce,
                           dimred = "UMAP",
                           colour_by = "Population",
                           point_alpha = 1)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")
UMAPplot

ggsave("Figure/II.3.UMAPplot.png", 
       plot=UMAPplot, width = 3, height = 2, dpi = 320)

## 4. Detect the PCA determinants  ---------

#Detect PCA determinants/gene that contribute to the most variance PC - Based on Rotation 
pcrot <- attr(reducedDim(gasBulk.sce),"rotation")
#Number of PC
npc <- 3
genenum <- 107
pclist <- list()
for (i in 1:npc) {
  pc <- round(pvplot$pv[i]/sum(pvplot$pv[1:npc])*genenum)
  genePC <- rownames(pcrot[order(abs(pcrot[,i]), decreasing = T)[1:pc],])
  pclist[[i]] <- genePC
}
#list of protein contribute to PC1, PC2, PC3
listPC <- c(pclist[[1]],pclist[[2]],pclist[[3]])
length(listPC)
#Find the top PCA determinants
PCAdeterminants <- which(match(rownames(gasBulk.sce), listPC)>0)
length(PCAdeterminants)

#Plot Expression ON PCAdeterminant - non-Scale 
PCAEplot <- plotExpression(gasBulk.sce,
                           rownames(gasBulk.sce)[head(PCAdeterminants,20)],
                           exprs_values = "assay",
                           log2_values = FALSE,
                           colour_by = "Population",)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population") +
  ggtitle("Expression of top 20 PCA determinent")
PCAEplot

ggsave("Figure/II.4. Top 20 PCA determinent - non-Scale.png", 
       plot=PCAEplot, width = 7,height = 6, dpi = 320)

#Plot Expression ON PCAdeterminant - Scale 
PCAEscaleplot <- plotExpression(gasBulk.sce, 
                           rownames(gasBulk.sce)[head(PCAdeterminants,20)],
                           exprs_values = "scaleassay",
                           log2_values = FALSE,
                           colour_by = "Population",)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle("Expression of top 20 PCA determinent - Scale logcounts")
PCAEscaleplot

ggsave("Figure/II.4. 20 PCA determinent - Scale.png", 
       plot=PCAEscaleplot, width = 7,height = 6, dpi = 320)


## 5. Heatmap for top PCA determinants/genes  -----------

#Choose Scale data for Heatmap
gasBulk.PCA.hm <- assay(gasBulk.sce,"scaleassay")[PCAdeterminants,]

#Or do it manually
#gasBulk.PCA.hm <- assay(gasBulk.sce,"logcounts")[PCAdeterminants,]
#scale_gasBulk = t(scale(t(gasBulk), center = TRUE, 
#                       scale = TRUE))

### Hierarchical clustering -----

# Clustering with different distance methods
#dist(1-cor(t(gasBulk.aov.hm))) 
#"pearson"
#spearman
#kendall
#euclidean

dist.method <- "euclidean"
gasBulk.PCA.hm.hc <- Heatmap(gasBulk.PCA.hm, name = "Scale Log Abundance",
                         clustering_distance_rows = dist.method,
                         clustering_method_rows = "ward.D2",
                         cluster_columns = dendsort(hclust(dist(t(gasBulk.PCA.hm)))),
                         #column_split = 4,
                         col = colour.hm,
                         heatmap_legend_param = list(direction = "horizontal"),
                         top_annotation = HeatmapAnnotation(Population = Population,
                                                            border = TRUE,gap = unit(0, "points"),
                                                            col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                 "Sox17-RFP+" = "#E03E3E",
                                                                                 "Bra-GFP+" = "#53BB79",
                                                                                 "Triple-Neg" = "#EF8228")),
                                                            show_legend = FALSE),
                         heatmap_width = unit(8, "cm"), 
                         heatmap_height = unit(18, "cm"),
                         column_names_gp = grid::gpar(fontsize = 5),
                         row_names_gp = grid::gpar(fontsize = 5),
                         column_title = paste("First", length(PCAdeterminants),"proteins contribute\n to PC1 and PC2"))

gasBulk.PCA.hm.hc

#
png(paste0("Figure/II.5. Heatmap on First ", length(PCAdeterminants)," proteins contribute to PC1 and PC2 - ",dist.method,".png"),width = 5.25,height = 8.25,units="in",res = 320)
draw(gasBulk.PCA.hm.hc, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
draw(lgChan, x = unit(1, "npc"), y = unit(0.665, "npc"), just = c("right", "top"))
dev.off()

### k-means clustering -----

n <- 4
kclus <- kmeans(gasBulk.PCA.hm, n)
split2 <- paste0("Cluster\n", kclus$cluster)
gasBulk.PCA.hm.km <- Heatmap(gasBulk.PCA.hm, name = "Scale Log \nAbundance",
                     row_dend_reorder = TRUE,
                     cluster_row_slices = FALSE,
                     split=split2, 
                     clustering_method_rows = "ward.D2",
                     cluster_columns = dendsort(hclust(dist(t(gasBulk.PCA.hm)))),
                     col = colour.hm,
                     top_annotation = HeatmapAnnotation(Population = as.character(Population),
                                                        border = TRUE,gap = unit(2, "points"),
                                                        col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                             "Sox17-RFP+" = "#E03E3E",
                                                                             "Bra-GFP+" = "#53BB79",
                                                                             "Triple-Neg" = "#EF8228")),
                                                        show_legend = FALSE),
                     column_names_gp = grid::gpar(fontsize = 6),
                     row_names_gp = grid::gpar(fontsize = 5),
                     column_title = paste("First", length(PCAdeterminants),"proteins contribute\n to PC1 and PC2 - k-means"))

gasBulk.PCA.hm.km

png(paste0("Figure/II.5.Heatmap on First", length(PCAdeterminants),"proteins contribute to PC1 and PC2 - k-means.png"),width = 5.25,height = 8.25,units="in",res = 320)
gasBulk.PCA.hm.km
draw(lgChan, x = unit(1, "npc"), y = unit(0.665, "npc"), just = c("right", "top"))
dev.off()

# III. DOWNSTREAM ANALYSES  ------------

#Choose Scale data for more dynamic
#.Find Gene Marker with scran
#Test <- c("t", "wilcox", "binom")
#pval.type <- c("any", "some", "all")
#t-test

markers <- findMarkers(gasBulk.sce,
                       assay.type = "assay",
                       groups=gasBulk.sce$Population,
                       pval.type="any",
                       lfc=1)

## 1. Heatmap for Marker proteins/genes ---------
#Layer
layer <- c("Mt1-BFP+","Bra-GFP+","Triple-Neg","Sox17-RFP+")

#choose n top genes
ntop <- 10

#Colour
col.hm <- colorRamp2(seq(-4, 4, length = 14), 
                     colour.hm)
#Legend
lgLFC = Legend(col_fun = col.hm, title = "Log2 FC",  direction = "horizontal")

#Generate loop for Heatmap
for (i in layer){
  layer.marker <- markers[[i]]
  top.marker <- layer.marker[layer.marker$Top <= ntop,]
  LFCsig <- getMarkerEffects(top.marker[1:ntop,])
  plot <- Heatmap(LFCsig,
                  col = col.hm,
                  column_title = paste0("Top ", ntop, " differential expressed proteins of ",i," \nover other populations"),
                  heatmap_legend_param = list(
                    #legend_height = unit(4, "cm")))
                    #title_position = "lefttop-rot"),
                    direction = "horizontal",
                    title = "Log2FC"),
                  show_heatmap_legend = T)
  assign(paste0("Marker",i),plot)
}

png(paste0("Figure/III.1.1. Heatmap on top ", ntop," significant Log2 Fold Change of Mt1-BFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerMt1-BFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. Heatmap on top ", ntop," significant Log2 Fold Change of Sox17-RFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerSox17-RFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. Heatmap on top ", ntop," significant Log2 Fold Change of Bra-GFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerBra-GFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. Heatmap on top ", ntop," significant Log2 Fold Change of Triple-Neg - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerTriple-Neg`, heatmap_legend_side = "bottom")
dev.off()


## 2. Detect all the marker genes--------

#List all marker gene
genemarker <- list()
for(i in layer){
  layer.marker <- markers[[i]]
  top.marker <- layer.marker[layer.marker$Top <= ntop,]
  genemarker[[i]] <- rownames(top.marker[,])
}

#Marker gene list 
EMK <- which(match(rownames(gasBulk.sce), genemarker %>% unlist)>0)
#Remove duplicated genes
#EMK <- EMK[-which(duplicated(rownames(gasBulk.sce)[EMK]))]

## 3. Plot all marker genes--------

### Genemarker plot -------
#non-Scale
GMplot <- plotExpression(gasBulk.sce, 
                         rownames(gasBulk.sce)[EMK],
                         exprs_values = "assay",#non-Scale
                         log2_values = FALSE, #already log2
                         one_facet=TRUE,
                         colour_by = "Population")+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ", ntop," marker per cell population"))
GMplot
ggsave("Figure/III.1.2. Expression of Marker genes - non-Scale.png",plot =GMplot,
       width = 10.25,height = 4.25,units="in",dpi = 320)



#scale
GMscaleplot <- plotExpression(gasBulk.sce, 
                         rownames(gasBulk.sce)[EMK],
                         exprs_values = "scaleassay", #Scale
                         log2_values = FALSE, #already log2
                         one_facet = TRUE,
                         colour_by = "Population",)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ", ntop," marker per cell population - Scale"))

GMscaleplot
ggsave("Figure/III.1.3. Expression of Marker genes - Scale.png", plot =GMscaleplot,
       width = 10.25,height = 4.25,units="in",dpi = 320)

### Dot plot -------

#non-Scale
dotEMK <- plotDots(gasBulk.sce, 
                   rownames(gasBulk.sce)[EMK],
                   exprs_values = "assay",
                   group="Population") + 
  ggtitle(paste0("Top ", ntop," markers per cell population"))+
  theme_bw()
dotEMK
ggsave("Figure/III.1.4. Dot Expression of Marker genes  - non-Scale.png",
       plot=dotEMK, width = 5,height = 8,units="in",dpi = 320)

#Scale
dotscaleEMK <- plotDots(gasBulk.sce, 
                   rownames(gasBulk.sce)[EMK],
                   exprs_values = "scaleassay",
                   group="Population") + 
  ggtitle(paste0("Top ", ntop," markers per cell population - Scale"))+
  theme_bw()
dotscaleEMK

ggsave("Figure/III.1.5. Dot Expression of Marker genes  - Scale.png",
       plot=dotscaleEMK, width = 5,height = 8,units="in",dpi = 320)

# IV. Integrate with scRNA-seq dataset ------

## 1. Reproduce UMAP from van der Brink et al. study --------
DEmarkers <- read.csv("BulkData/Fig1_41586_2020_2024_MOESM18_ESM.csv")
p <- ggplot(data = DEmarkers) +
  geom_point(mapping = aes(x = u1, y = u2, color = as.factor(celltype)),
             size=0.5)
p + theme_void()

ggsave("Figure/IV.1. UMAP plot form Brink - ALL .png", plot = p + theme_void(),
       width=7,height=4,units="in",dpi=320)

#DE gene list
ge1 <- read.csv("BulkData/geneExtendedDataFig.2.csv")
ge2 <- read.csv("BulkData/geneExtendedDataFig.2b.csv")
#Combine datatsets
ge <- cbind(ge1, ge2[,4:26])

# Setting colour 
colfunc <- colorRampPalette(c("white","blue"))

#Manual

#Sox17
Level <- ge$Sox17
Sox17 <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Sox17 +labs(title = "Sox17") + theme_void()

ggsave("Figure/IV.1. UMAP plot form Brink Sox17.png", plot = Sox17 +labs(title = "Sox17") + theme_void(),
       width=6,height=4,units="in",dpi=320)

#Sox2
Level <- ge$Sox2
Sox2 <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Sox2 +labs(title = "Sox2") + theme_void()

ggsave("Figure/IV.1. UMAP plot form Brink Sox2.png", plot = Sox2 +labs(title = "Sox2") + theme_void(),
       width=6,height=4,units="in",dpi=320)


#Mt1
Level <- ge$Mt1
Mt1 <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Mt1 + theme_void()+ labs(title = "Mt1")

ggsave("Figure/IV.1. UMAP plot form Brink Mt1.png", plot = Mt1 +labs(title = "Mt1") + theme_void(),
       width=6,height=4,units="in",dpi=320)



#Utf1
Level <- ge$Utf1
Utf1 <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Utf1 + theme_void()+ labs(title = "Utf1")

ggsave("Figure/IV.1. UMAP plot form Brink Utf1.png", plot = Utf1 +labs(title = "Utf1") + theme_void(),
       width=6,height=4,units="in",dpi=320)

#Tbx6
Level <- ge$Tbx6
Tbx6 <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Tbx6 + theme_void()+ labs(title = "Tbx6")

ggsave("Figure/IV.1. UMAP plot form Brink Tbx6.png", plot = Tbx6 +labs(title = "Tbx6") + theme_void(),
       width=6,height=4,units="in",dpi=320)

#T
Level <- ge$T
Bra <- ggplot(data = ge) +
  geom_point(mapping = aes(x = u1, y = u2, color = Level), size=0.5) +
  scale_colour_gradientn(colours = colfunc(10)[2:10])
Bra + theme_void()+ labs(title = "T/Bra")

ggsave("Figure/IV.1. UMAP plot form Brink T or Bra .png", plot = Bra +labs(title = "T/Bra") + theme_void(),
       width=6,height=4,units="in",dpi=320)


## 2.Integrate with scRNA-seq dataset from van der Brink study --------

### 2.1. Integrate Data-----

# Load lists from four cluster in 120h mouse gastruloids 
DEG_endoderm <- read.delim("BulkData/DEG_endoderm.txt", header = F)
DEG_ectoderm <- read.delim("BulkData/DEG_ectoderm.txt", header = F)
DEG_PSM <- read.delim("BulkData/DEG_PSM.txt", header = F)
DEG_NMP <- read.delim("BulkData/DEG_NMPs.txt", header = F)

#Create data frame contains all clusters
DEG_all <- data.frame(gene = c(DEG_ectoderm$V1,
                               DEG_endoderm$V1,
                               DEG_PSM$V1,
                               DEG_NMP$V1),
                      cluster = c(rep("Ecto", nrow(DEG_ectoderm)),
                                  rep("Endo", nrow(DEG_endoderm)),
                                  rep("PSM", nrow(DEG_PSM)),
                                  rep("NMP", nrow(DEG_NMP))))

#Choose only unique gene in van den Brink
DEG_all<- DEG_all[!duplicated(DEG_all[ , 1]),]


# Annotated each file with cluster
DEG_endoderm <- DEG_all %>% filter(cluster == "Ecto")
DEG_ectoderm <- DEG_all %>% filter(cluster == "Endo") 
DEG_PSM <- DEG_all %>% filter(cluster == "PSM") 
DEG_NMP <- DEG_all %>% filter(cluster == "NMP")

#Find overlapped genes of two datasets for each clusters

## Endoderm
tmp12 <- which(match(gasBulk.aov[,"Gene.names"], DEG_endoderm$gene)>0)
length(tmp12)
tmp12a <- which(match(DEG_endoderm$gene, gasBulk.aov[,"Gene.names"])>0)
length(tmp12a)
dim(DEG_endoderm[tmp12a, ])

## Ectoderm
tmp13 <- which(match(gasBulk.aov[,"Gene.names"], DEG_ectoderm$gene)>0)
length(tmp13)
tmp13a <- which(match(DEG_ectoderm$gene, gasBulk.aov[,"Gene.names"])>0)
length(tmp13a)
dim(DEG_ectoderm[tmp13a, ])

## PSM
tmp6 <- which(match(gasBulk.aov[,"Gene.names"], DEG_PSM$gene)>0)
length(tmp6)
tmp6a <- which(match(DEG_PSM$gene, gasBulk.aov[,"Gene.names"])>0)
length(tmp6a)
dim(DEG_PSM[tmp6a, ])

## NMPs
tmp7 <- which(match(gasBulk.aov[,"Gene.names"], DEG_NMP$gene)>0)
length(tmp7)
tmp7a <- which(match(DEG_NMP$gene, gasBulk.aov[,"Gene.names"])>0)
length(tmp7a)
dim(DEG_NMP[tmp7a, ])

#Overlapped genes in each cluster

##With gasBulk
gasBulk.aov.DEG <- rbind(gasBulk.aov[tmp12,],
                   gasBulk.aov[tmp13,],
                   gasBulk.aov[tmp6,],
                   gasBulk.aov[tmp7,])
dim(gasBulk.aov.DEG)

##With van der Brink data
Brink.DEG <- rbind(DEG_endoderm[tmp12a,],
                     DEG_ectoderm[tmp13a,],
                     DEG_PSM[tmp6a,],
                     DEG_NMP[tmp7a,]) 
dim(Brink.DEG)

### 2.2. Heatmap for overlapped gene --------

#Choose dataframe of quantitative values and geneID
gasBulk.aov.Brink <- gasBulk.aov.DEG[,c(1:12,37)]

#Transfer into numeric matrix
gasBulk.aov.Brink.hm <- as.matrix(sapply(gasBulk.aov.Brink[1:12], as.numeric))

#Scale matrix
gasBulk.aov.Brink.hm = t(scale(t(gasBulk.aov.Brink.hm)))

#Choose overlapped genes (for adding into Heatmap)
overlapped.genes <- gasBulk.aov.DEG[,37]

#### Hierarchical clustering -----
#dist(1-cor(t(gasBulk.aov.hm))) 
#"pearson"
#spearman
#kendall
#euclidean
lgChan = Legend(labels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg"), 
                legend_gp = gpar(fill = colour.Chan), 
                title = "Population",
                direction = "horizontal")
dist.method <- "euclidean"
gasBulk.aov.Brink.hm.hc <- Heatmap(gasBulk.aov.Brink.hm, name = "Scale Log Abundance",
                                   clustering_distance_rows = dist.method,
                                   clustering_method_rows = "ward.D",
                                   cluster_columns = dendsort(hclust(dist(t(gasBulk.aov.Brink.hm)))),
                                   #row_split = 4,
                                   #column_split = 4,
                                   col = colour.hm,
                                   heatmap_legend_param = list(direction = "horizontal"),
                                   top_annotation = HeatmapAnnotation(Population = Population,
                                                                      border = TRUE,
                                                                      gap = unit(2, "points"),
                                                                      col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                           "Sox17-RFP+" = "#E03E3E",
                                                                                           "Bra-GFP+" = "#53BB79",
                                                                                           "Triple-Neg" = "#EF8228")),
                                                                      show_legend = FALSE),
                                   right_annotation = rowAnnotation(Brink = Brink.DEG$cluster,
                                                                    border = TRUE,
                                                                    annotation_legend_param = list(
                                                                      Brink = list(Brink = c("Ecto" = "dodgerblue2",
                                                                                            "PSM" = "#B7D468",
                                                                                            "NMP" = "#2EC20A",
                                                                                            "Endo" = "#D14239"),
                                                                                  direction = "horizontal")),
                                                                    col = list(Brink = c("Ecto" = "dodgerblue2",
                                                                                           "PSM" = "#B7D468",
                                                                                           "NMP" = "#2EC20A",
                                                                                           "Endo" = "#D14239"),
                                                                               direction = "horizontal"),
                                                                    show_legend = FALSE),
                                   column_title = paste0("Overlapped ",nrow(gasBulk.aov.Brink)," protein/mRNA \nwith van den Brink et al."),
                                   column_names_gp = grid::gpar(fontsize = 15),
                                   row_names_gp = grid::gpar(fontsize = 6)) +
  rowAnnotation(gene = anno_text(overlapped.genes, which = "row",gp = gpar(fontsize = 5)),
                width = max_text_width(unlist(overlapped.genes)) + unit(0.1, "cm"))
gasBulk.aov.Brink.hm.hc



png(paste0("Figure/IV.2. Heatmap on Overlapped ",nrow(gasBulk.aov.Brink)," protein-mRNA with van den Brink et al - remove duplicate.png"), width = 5.25, height = 8.25, units="in", res = 320)
draw(gasBulk.aov.Brink.hm.hc, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
#draw(lgChan, x = unit(1, "npc"), y = unit(0.715, "npc"), just = c("right", "bottom"))
dev.off()


#### k-means clustering -----
set.seed(120)
n <- 4
kclus <- kmeans(gasBulk.aov.Brink.hm, n)
#Heatmap of 991 significant hits
split2 <- paste0("Cluster\n", kclus$cluster)
gasBulk.aov.Brink.hm.km <- Heatmap(gasBulk.aov.Brink.hm, name = "Scale Log\nAbundance",
                                row_dend_reorder = TRUE,
                                cluster_row_slices = FALSE,
                                split=split2, 
                                clustering_method_rows = "ward.D2",
                                cluster_columns = dendsort(hclust(dist(t(gasBulk.aov.Brink.hm)))),
                                #column_split = 4,
                                col = colour.hm,
                                top_annotation = HeatmapAnnotation(Population = Population,
                                                                   border = TRUE,
                                                                   #gap = unit(2, "points"),
                                                                   col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                        "Sox17-RFP+" = "#E03E3E",
                                                                                        "Bra-GFP+" = "#53BB79",
                                                                                        "Triple-Neg" = "#EF8228")),
                                                                   show_legend = FALSE),
                                right_annotation = rowAnnotation(Brink = Brink.DEG$cluster,
                                                                 border = TRUE,
                                                                 col = list(Brink = c("Ecto" = "dodgerblue2",
                                                                                      "PSM" = "#B7D468",
                                                                                      "NMP" = "#2EC20A",
                                                                                      "Endo" = "#D14239"))),
                                column_title = paste0("Overlapped ",nrow(gasBulk.aov.Brink)," protein/mRNA \nwith van den Brink et al. - k-means"),
                                column_names_gp = grid::gpar(fontsize = 6),
                                row_names_gp = grid::gpar(fontsize = 1)) +
  rowAnnotation(gene = anno_text(overlapped.genes, which = "row",gp = gpar(fontsize = 4)),
                width = max_text_width(unlist(overlapped.genes)) + unit(0.1, "cm"))

gasBulk.aov.Brink.hm.km

png(paste0("Figure/IV.2. Heatmap on Overlapped ",nrow(gasBulk.aov.Brink)," protein-mRNA with van den Brink et al. - k-means.png"), width = 5.25, height = 8.25, units="in", res = 320)
gasBulk.aov.Brink.hm.km
draw(lgChan, x = unit(1, "npc"), y = unit(0.715, "npc"), just = c("right", "top"))
dev.off()
