packages <- c(
  ##Data manipulation and analytics
  "magrittr","dplyr","knitr","tidyverse","scales",
  "scp","scater","scran","sva","SC3",
  #Visualization
  "RColorBrewer","dendsort","cowplot", "ggplot2","ComplexHeatmap","circlize",
  #Downstream analysis
  "enrichplot","ggupset","clusterProfiler",
  #Mouse database
  "org.Mm.eg.db")
  
lapply(packages, require, character.only = TRUE)

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("Gas.100Data/gas.100.ProcessedData.Rda")

# I. DATA PROCESSING ----------

## Input ----------

###Input data

#Process MQ files from a MQ run
gas.100 <- read.delim("Gas.100Data/2&3. evidence - 100 cells - 0386 + 0370.txt")
expdeg <- "100 Cells"

#Cleaning
#Remove >2nd leading proteins
gas.100$Leading.razor.protein <- sub(';.*','',gas.100$Leading.razor.protein)
#Remove >2nd leading genes
gas.100$Gene.names <- sub(';.*','',gas.100$Gene.names)

#Combines 2 CVs into one (No need if set names for "Experiments" columns):
#gas.100$Raw.file <- sub('_FS45','',gas.100$Raw.file)
#gas.100$Raw.file <- sub('_FS65','',gas.100$Raw.file)

#Annotation file
sampleAnnotation <- read.csv("Gas.100Data/2&3. sampleAnnotation - Gas_100.csv")

#Sox17 - Mt - Bra
RePro <- c("Q61473","P02802","P20293")
Rege <- c("Sox17", "Mt1", "Tbxt", "T", "Mt","Bra")


##Colour
#For heatmap
colour.hm <-  c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
#For Channels
colour.Chan <- c("dodgerblue","#E03E3E","#53BB79","#EF8228","#937264","#3043A2","#C25D7B")
#For Datasets
colour.100 <- c("#D88C4C","#E49D41","#FEB447","#06CDE0","#1498BE")


## I.1. PSM level ----------

### 1. Reconstructing the dataset to create SCP/QFeatures object ----------
#In this research, all the experiments were TMT16-plex
#read SCP
scp <- readSCP(featureData = gas.100,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Experiment",
               #removeEmptyCols = TRUE
               )
#Number of assays
length(scp)

### 2. Quality control at PSM level ----------

#### 2.1. Cleaning missing data ----------

#The zeros can be biological zeros or technical zeros.
#Therefore any zero should be replaced by NA to avoid artefacts in downstream analysis.
scp <- zeroIsNA(scp, i = names(scp))

#### 2.2. Filter out failed runs based on PSM content (detected features) ----------

#Plot for checking number of PSM per assays
#nPSM <- dims(scp)[1, which(startsWith(names(scp),"D"))]
nPSM <- dims(scp)[1,]

#Create data frame
nPSM <- data.frame(Datasets = names(nPSM),
                       Count = nPSM)
nPSM$Datasets <- factor(nPSM$Datasets, levels = nPSM$Datasets)

#Plot
nPSMplot <- ggplot(data=nPSM, aes(x=Datasets, y=Count, fill=Datasets)) +
  scale_fill_manual(values= colour.100, 
                    labels= nPSM$Datasets) +
  geom_bar(stat="identity")+
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
  ylim(0,6000)+
  xlab("Datasets") + 
  ylab("Number of detected PSM") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank())
nPSMplot

ggsave(paste0("Figure/2.2. ",expdeg," - Number of detected peptides per assay.png"), 
       plot=nPSMplot, width=4, height=4, dpi=320)


#Select the assays that have sufficient PSMs - dont need to perform in 100-cell
#e.g. the number of rows is greater than 500),
#keepAssay <- dims(scp)[1, ] > 500
#scp <- scp[, , keepAssay]

#### 2.3. Filter PSMs for contaminants and noisy spectra ----------

scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        !is.na(PIF) & PIF > 0.8)
#or
#scp <- filterFeatures(scp,
#                      ~ !grepl("REV|CON", protein)) & 
#                       !is.na(PIF) & PIF > 0.8)

#### 2.4. Filter features to control for high false discovery rate FDR ----------

#pep2qvalue function to convert PEPs to q-values that are directly related to FDR
#Either
scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "PEP",
                  rowDataName = "qvalue_PSMs")

#or compute q-values at peptide or protein level rather than PSM
scp <- pep2qvalue(scp,
                  i = names(scp),
                  PEP = "PEP",
                  groupBy = "Leading.proteins",
                  rowDataName = "qvalue_proteins")


#Plot 
qvalueplot <- rowDataToDF(scp,
            i = names(scp),
            vars = c("qvalue_PSMs", "qvalue_proteins")) %>%
  data.frame %>%
  pivot_longer(cols = c("qvalue_PSMs", "qvalue_proteins"),
               names_to = "measure") %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0.01) +
  scale_x_log10() +
  facet_grid(rows = vars(measure)) +
  theme_bw()
qvalueplot

ggsave(paste0("Figure/2.4. ",expdeg," - q-value plot.png"), 
       plot=qvalueplot, width=7.2,height=4.45, dpi=320)

#Filter the PSM to control, the protein FDR at 1%
scp <- filterFeatures(scp,
                      ~ qvalue_proteins < 0.01)

#or
#scp <- filterFeatures(scp,
#                      ~ qvalue_psm < 0.01 & qvalue_protein < 0.01)

#### 2.5. Filter out PSMs with high sample to carrier ratio (for single-cell) ----------

#This experiment was performs on 100 cells, so there is no filtering features based on SCP metrics (single-cells metrics on Carrier and Reference channels)

#The number of channels per cell population (100 cells per channel) 

#### Number channels per pell population ----------
poputable <- t(as.data.frame(table(colData(scp)[, "Population"])))
colnames(poputable) <- poputable[1,]
poputable <- poputable[-1, ]
poputable <- t(as.data.frame(poputable))
rownames(poputable) <- c("Channels")
knitr::kable(poputable)

#Save
scp2 <- scp
scp <- scp2

## I.2. Peptide level ----------


### 4. Aggregate PSM data to peptide data ----------
scp <- aggregateFeaturesOverAssays(scp,
                                   i = names(scp),
                                   fcol = "Sequence",
                                   name = pepAssays <- paste0("Peptides_", names(scp)),
                                   fun = matrixStats::colMedians, na.rm = TRUE)

### 5. Join all datasets into one assay ----------
scp <- joinAssays(scp,
                  i = pepAssays,
                  name = "Peptides")

#Optional: Filter cell of interest
#scp <- scp[, scp$Population %in% c("Sox17-RFP+", "Mt1-BFP+", "Bra-GFP+", "Triple-Neg", "Empty"), ]

#Save
scp5 <- scp
scp <- scp5

### 6. (Single-cell/100 cells) quality control ----------

#### 6.1. Median Relative Reporter Intensity ----------

#Computing median relative reporter ion intensity for each channels/cells separately 
#and apply a filter based on this statistic. 

#Compute MedianRI
MedianRI <- apply(assay(scp,"Peptides"),
                  MARGIN = 2, 
                  FUN = median, 
                  na.rm = TRUE)
#Median Relative Reporter Intensity will be stored in the colData. 
colData(scp)[names(MedianRI), "MedianRI"] <- MedianRI

#Make Population as factor for plotting
scp$Population <- factor(scp$Population,
                         levels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg","Empty"))
#Plot
MedianRIplot <- scp %>% colData %>% data.frame %>%
  ggplot(aes(x = Population, y = MedianRI, fill = Population)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.025) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Dataset),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.100) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  #geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianRI of Peptide Intensity"))
MedianRIplot

ggsave(paste0("Figure/6.1. ",expdeg," - MedianRI for Peptide Intensity.png"),
       plot=MedianRIplot, width=7.2,height=4.45,units="in",dpi=320)

#Save
scp6.1 <- scp
#### 6.2. Median Coefficient of Variation ----------

#The Median Coefficient of Variation measures the consistency of quantification 
#for a group of peptides that belong to a protein. 

#medianCVperCell: computes CV for each protein in each channel/cell.
#groupBy: name of the rowData field containing the protein information
#nobs: compute CVs for proteins that have at least n peptides 
#norm: #"sum", "max", "center.mean", "center.median", "div.mean", 
##"div.median", "diff.meda", "quantiles", "quantiles.robust"  "vsn"
#SCoPE2

#Median Coefficient of Variation will be stored in the colData. 
#Note: Choose peptide assays (pepAssays) for computing MedianCV
scp <- scp6.1

scp <- medianCVperCell(scp,
                       i = pepAssays, #or (l+1):(l*2)
                       groupBy = "Leading.razor.protein",
                       nobs = 6, 
                       norm = "SCoPE2",
                       colDataName = "MedianCV")

#Threshold for minimum of MedianCV for Empty Channel
ts <- min(scp$MedianCV[scp$Population=="Empty"],na.rm = TRUE) - 0.001 # - 0.001 for drawinga line below the minium value

#Plot
MedianCVplot <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = Population, y = MedianCV, fill = Population)) +
  geom_boxplot(width=0.4,outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Dataset),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.100)+
  #labs(colour = "Datasets")+
  ylim(0,1) +
  geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianCV of Peptide Intensity"))
MedianCVplot

ggsave(paste0("Figure/6.2. ", expdeg," - MedianCV for Peptide Intensity - remove outlier.png"), 
       plot=MedianCVplot, width=7.2,height=4.45,units="in",dpi=320)

#MedianCV in channels containning cell(s) should theoretically much more consistent 
#than for Empty channels.

scp6.2 <- scp

####  6.3. Filter based on the MedianRI and MedianCV ----------

#Based on the distribution of the Empty, ones can decides which MedianCV threshold and 
#remove channels/cells that exhibit high MedianCV over the different proteins. 
scp <- scp[, !is.na(scp$MedianCV) & scp$MedianCV < ts, ]

#Optional: to remove some channels/cells
scp <- scp[, scp$Population != "Empty", ]

#Plot after filtering

#MedianRI
MedianRIplotQC <- scp %>% colData %>% data.frame %>%
  ggplot(aes(x = Population, y = MedianRI, fill = Population)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.025) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Dataset),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.100) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  #geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianRI of Peptide Intensity after QC"))
MedianRIplotQC

ggsave(paste0("Figure/6.3. ", expdeg," - MedianRI for Peptide Intensity after QC.png"), 
       plot=MedianRIplotQC, width=7.2,height=4.45,units="in",dpi=320)


#MedianCV
MedianCVplotQC <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = Population, y = MedianCV, fill = Population)) +
  geom_boxplot(width=0.4) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Dataset),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.100)+
  ylim(0,1) +
  geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianCV of Peptide Intensity after QC"))
MedianCVplotQC

ggsave(paste0("Figure/6.3. ", expdeg," - MedianCV Peptide Intensity after QC.png"), 
       plot=MedianCVplotQC, width=7.2,height=4.45,units="in",dpi=320)


#Return Population to character (required for Batch correction)
scp$Population <- as.character(scp$Population)

#Save
scp6 <- scp

#### Number channels per pell population ----------
poputable <- t(as.data.frame(table(colData(scp)[, "Population"])))
colnames(poputable) <- poputable[1,]
poputable <- poputable[-1, ]
poputable <- t(as.data.frame(poputable))
rownames(poputable) <- c("Channels")
knitr::kable(poputable)


### 7. Normalization ----------

# Divide columns by median
scp <- sweep(scp, 
             i = "Peptides",
             MARGIN = 2,
             FUN = "/",
             STATS = colMedians(assay(scp[["Peptides"]]), na.rm = TRUE),
             name = "Peptides_norm_col")
# Divide rows by mean
scp <- sweep(scp,
             i = "Peptides_norm_col",
             MARGIN = 1,
             FUN = "/",
             STATS = rowMeans(assay(scp[["Peptides_norm_col"]]), na.rm = TRUE),
             name = "Peptides_norm")

### 8. Filter peptides based on missing data ----------
scp <- filterNA(scp,
                i = "Peptides_norm",
                pNA = 0.99)

### 9. Log-transformation ----------
scp <- logTransform(scp,
                    base = 2,
                    i = "Peptides_norm",
                    name = "Peptides_log")

#Save
scp9 <- scp
scp <- scp9

## I.3. Protein level ----------

### 10. Aggregate peptide data to protein data ----------

scp <- aggregateFeatures(scp,
                         i = "Peptides_log",
                         name = "Proteins",
                         fcol = "Leading.razor.protein",
                         fun = matrixStats::colMedians, na.rm = TRUE)

scp10 <- scp
scp <- scp10

### 11. Normalization ----------

scp %>%
  ## Center columns with median
  sweep(i = "Proteins",
        MARGIN = 2,
        FUN = "-",
        STATS = colMedians(assay(scp[["Proteins"]]),
                           na.rm = TRUE),
        name = "Proteins_norm_col") %>%
  ## Center rows with mean
  sweep(i = "Proteins_norm_col",
        MARGIN = 1,
        FUN = "-",
        STATS = rowMeans(assay(.[["Proteins_norm_col"]]),
                         na.rm = TRUE),
        name = "Proteins_norm") ->
  scp

#save
scp11 <- scp
scp <- scp11


### 12. Imputation (using the KNN algorithm) ----------
misval <- longFormat(scp[, , "Proteins_norm"]) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram() +
  xlim(0, 1) +
  theme_bw()
misval

ggsave(paste0("Figure/12. ", expdeg," - missing value.png"), 
       plot=misval, width=7.2,height=4.45,units="in",dpi=320)

#Before imputation, there percentage of missingness: 
assay(scp,"Proteins_norm") %>% is.na %>%  mean

#Imputation k = 4
scp <- impute(scp,
              i = "Proteins_norm",
              method = "knn",
              k = 4, rowmax = 1, colmax= 1,
              maxp = Inf, rng.seed = 1234)

#Check after imputation, missing value is equal to (If 0 then fine):
assay(scp,"Proteins_norm") %>% is.na %>%  mean

#Save
scp12 <- scp
scp <- scp12



### 13. Batch correction (using the ComBat algorithm) ----------

sce <- getWithColData(scp, "Proteins_norm")
batch <- colData(sce)$Raw.file
model <- model.matrix(~ as.character(Population), data = colData(sce))

assay(sce) <- ComBat(dat = assay(sce),
                     batch = batch,
                     mod = model)
addAssay(scp,
         y = sce,
         name = "Proteins_batchC") %>%
  addAssayLinkOneToOne(from = "Proteins_norm",
                       to = "Proteins_batchC") ->
  scp

#Save
scpB <- scp
scp <- scpB


### 14. Reconstructing the dataset for Downstream analysis ----------

#Picking only of Experiments/Assays of interest

gas.DR <- getWithColData(scp, "Proteins_batchC")

#Add a scale-log assays for some heatmap/visualization
assay(gas.DR,"scaleassay") <- t(scale(t(assay(gas.DR,"assay")), 
                                      center = TRUE, scale = TRUE))

#Find proteins that have no correspoding gene
blankgene <- gas.100[match(rownames(gas.DR), gas.100$Leading.razor.protein),][,c("Leading.razor.protein","Gene.names")]

#Search on  Uniprot
blankgene[which(blankgene[,"Gene.names"] == ""),]
#2836 - A0A2I3BQR1 <- Nt5dc2
#6472 - A0A3B2WBC6 <- Polr2m
#3545- Q9D735 <- Trir
#2681 - C0HKD9 <- Mfap1b
#3809 - Q9CQE8 <- RTRAF
#Manual add genes
gas.100[2836,"Gene.names"] <- "Nt5dc2"
gas.100[6472,"Gene.names"] <- "Polr2m"
gas.100[3545,"Gene.names"] <- "Trir"
gas.100[2681,"Gene.names"] <- "Mfap1b"
gas.100[3809,"Gene.names"] <- "RTRAF"

#Change rownames gas.DR from proteins names to correspoding gene names
rownames(gas.DR) <- gas.100[match(rownames(gas.DR), gas.100$Leading.razor.protein),]$Gene.names

#Check order
#View(data.frame(rownames(scp[["Proteins_batchC"]]),
#                rownames(gas.DR),
#                gas.100[match(rownames(scp[["Proteins_batchC"]]), gas.100$Leading.razor.protein),]$Gene.names,d,
#                gas.100[match(rownames(scp[["Proteins_batchC"]]), gas.100$Leading.razor.protein),]$Leading.razor.protein))

save(gas.100,
     expdeg,
     sampleAnnotation,
     scp,
     gas.DR,
     MedianCVplot,
     MedianRIplot,
     colour.hm,
     colour.Chan,
     colour.100,
     file="Gas.100Data/gas.100.ProcessedData.Rda")

#Save
#gas.DR1 <- gas.DR

# II. DIMENSIONALITY REDUCTION AND VISUALIZATION ----------

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("Gas.100Data/gas.100.ProcessedData.Rda")
#Change order for PCA plot
colour.PCA <- c("#53BB79","dodgerblue","#E03E3E","#EF8228")

## 1. Principal component analysis ----------

#Choose n top proteins, scale/ non-Scale may generate different results
set.seed(100)
PCAtop <- 500
gas.DR <- runPCA(gas.DR,
                 ncomponents = 10,
                 ntop = PCAtop,
                 scale = FALSE,
                 exprs_values = "assay",
                 name = "PCA")

#Note: If scale=TRUE, the expression values for each feature are standardized 
#so that their variance is unity. This will also remove features with standard deviations below 1e-8.
#1.2. PCA plot
PCAplot <- plotReducedDim(gas.DR,
               dimred = "PCA",
               colour_by = "Population",
               point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR), "percentVar")[2],1), "%)"))
PCAplot

ggsave(paste0("Figure/II.1.1. ",expdeg," - PCA for top ",PCAtop," PCA determinants 2.png"), 
       plot=PCAplot, width=3, height=2, dpi=320)


#1.2. PCA Check Batch
gas.DR$lcbatch <- as.factor(gas.DR$lcbatch ) 
PCAplotbatch <- plotReducedDim(gas.DR,
                          dimred = "PCA",
                          colour_by = "lcbatch",
                          point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = c("#D88C4C","#1498BE"))+
  labs(colour = "Batch")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR), "percentVar")[2],1), "%)"))
PCAplotbatch

ggsave(paste0("Figure/II.1.2. ",expdeg," - PCA Batch for top ",PCAtop," PCA determinants.png"), 
       plot=PCAplotbatch, width=2.8, height=2.2, dpi=320)


#1.3. No batch 

#Take assay without batch
gas.DR.nB <- getWithColData(scp, "Proteins_norm")

set.seed(100)
PCAtop <- 500
gas.DR.nB <- runPCA(gas.DR.nB,
                 ncomponents = 10,
                 ntop = PCAtop,
                 scale = FALSE,
                 exprs_values = "assay",
                 name = "PCA")

#Note: If scale=TRUE, the expression values for each feature are standardized 
#so that their variance is unity. This will also remove features with standard deviations below 1e-8.
#1.4. PCA plot - No Batch
PCAplotNB <- plotReducedDim(gas.DR.nB,
                          dimred = "PCA",
                          colour_by = "Population",
                          point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  theme(legend.position="none")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[2],1), "%)"))
PCAplotNB

ggsave(paste0("Figure/II.1.3. ",expdeg," - PCA for top ",PCAtop," PCA determinants - No batch correction.png"), 
       plot=PCAplotNB, width=2, height=2, dpi=320)


#1.4. PCA Check Batch - No Batch
gas.DR.nB$lcbatch <- as.factor(gas.DR.nB$lcbatch ) 
PCAplotbatchNB <- plotReducedDim(gas.DR.nB,
                               dimred = "PCA",
                               colour_by = "lcbatch",
                               point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = c("#D88C4C","#1498BE"))+
  labs(colour = "Batch")+
  theme(legend.position="none")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[2],1), "%)"))
PCAplotbatchNB

ggsave(paste0("Figure/II.1.4. ",expdeg," - PCA Batch for top ",PCAtop," PCA determinants - No batch correction.png"), 
       plot=PCAplotbatchNB, width=2, height=2, dpi=320)


## 2. Percentage of variance explained ----------

#Percentage of variance explained by each PCs
pv <- attr(reducedDim(gas.DR), "percentVar")
pvplot <- data.frame(c(seq(1,length(pv))),pv)

#Find elbow percentage of variance explained
elbow <- PCAtools::findElbowPoint(pv)

#Plot 
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

ggsave(paste0("Figure/II.2.1. ",expdeg," - Percentage of variance explained Plot.png"), 
       plot=pvPlot, width=2.8,height=2.2, dpi=320)

#Top PCs (elbows) explain for %
pvplot$pv[1:3] %>% sum
pvplot$pv[1:elbow] %>% sum

#Take only 3 PC
reducedDim(gas.DR, "PCA.elbow") <- reducedDim(gas.DR)[,1:elbow]
#reducedDimNames(gas.DR)

#Check first 3 PC correlation
PCAtoplot <- plotPCA(gas.DR, 
                     ncomponents = elbow, 
                     colour_by = "Population")+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  theme(legend.position="none")
PCAtoplot

ggsave(paste0("Figure/II.2.2. ",expdeg," - PCA with more PCs plot.png"), 
       plot=PCAtoplot, width=4,height=4, dpi=320)

## 3. Uniform Manifold Approximation and Projection ----------

#850
set.seed(100)
gas.DR <- runUMAP(gas.DR,
                  ncomponents = 3,
                  exprs_values = "assay",
                  n_neighbors = 4,
                  dimred = "PCA.elbow",
                  name = "UMAP")
UMAPplot <- plotReducedDim(gas.DR,
               dimred = "UMAP",
               colour_by = "Population",
               point_alpha = 0.5)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")
UMAPplot

ggsave(paste0("Figure/II.3. ",expdeg," - UMAPplot.png"), 
              plot=UMAPplot, width=3,height=2,units="in",dpi=320)

## 4. Detect the PCA determinants ----------

## Detect PCA determinants/gene that contribute to the most variance PC - Based on Rotation 
pcrot <- attr(reducedDim(gas.DR),"rotation")
#Number of PC
npc <- 4
genenum <- 118
pclist <- list()
for (i in 1:npc) {
  pc <- round(pvplot$pv[i]/sum(pvplot$pv[1:npc])*genenum)
  genePC <- rownames(pcrot[order(abs(pcrot[,i]), decreasing = T)[1:pc],])
  pclist[[i]] <- genePC
}
#list of protein contribute to PC1, PC2, PC3
listPC <- c(pclist[[1]],pclist[[2]],pclist[[3]],pclist[[4]])
length(listPC)

#Find the top PCA determinants - remove last one to keep top 100 determinant
PCAdtmn <- which(match(rownames(gas.DR), listPC[-119])>0)
length(PCAdtmn)


#Plot Expression of these PCA determinants
ntop <- 20

#non-Scale
PCAEplot <- plotExpression(gas.DR,
                           rownames(gas.DR)[head(PCAdtmn,20)],
                           exprs_values = "assay",
                           log2_values = FALSE,
                           colour_by = "Population",) +
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ",ntop," PCA determinent"))
PCAEplot

ggsave(paste0("Figure/II.4. ",expdeg," - Top ",ntop," PCA determinent - non-Scale.png"), 
       plot=PCAEplot, width=7,height=6, dpi=320)

#Scale
PCAEscaleplot <- plotExpression(gas.DR, 
                           rownames(gas.DR)[head(PCAdtmn,20)],
                           exprs_values = "scaleassay",
                           log2_values = FALSE,
                           colour_by = "Population",)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ",ntop," PCA determinent - Scale"))
PCAEscaleplot

ggsave(paste0("Figure/II.4. ",expdeg," - Top ",ntop," PCA determinent - Scale.png"), 
       plot=PCAEscaleplot, width=7,height=6, dpi=320)


## 5. Heatmap for PCA determinants/genes ----------

#Take data of n top PC determinants - Scale data (non-Scale is also possible)
gas.HM <- assay(gas.DR,"scaleassay")[PCAdtmn,]

#Shorten col names
colnames(gas.HM) <- sub('D0370_100_1Reporter.intensity.*','D0370_100_1',colnames(gas.HM))
colnames(gas.HM) <- sub('D0370_100_2Reporter.intensity.*','D0370_100_2',colnames(gas.HM))
colnames(gas.HM) <- sub('D0370_100_3Reporter.intensity.*','D0370_100_3',colnames(gas.HM))
colnames(gas.HM) <- sub('D0386_100_1Reporter.intensity.*','D0386_100_1',colnames(gas.HM))
colnames(gas.HM) <- sub('D0386_100_2Reporter.intensity.*','D0386_100_2',colnames(gas.HM))

#Create Heatmap for top PC determinants
#Manually create legend of Cells Population and Datasets
lgChan <- Legend(labels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg"), 
                legend_gp = gpar(fill = colour.Chan), 
                title = "Population")
lgData <- Legend(labels = c("0370_1","0370_2","0370_3","0386_1","0386_2"),
                 legend_gp = gpar(fill = c("#D88C4C","#E49D41","#FEB447","#06CDE0","#1498BE")),
                 title = "Datasets")
pL = packLegend(list = list(lgChan, lgData))

### Hierarchical clustering ----------

#(Wihout spliting: remove "column_split = 4,")
gas.HM.hc <- Heatmap(gas.HM, name = "Scale Log \nAbundance",
                     clustering_distance_rows = "pearson",
                     clustering_method_rows = "ward.D2",
                     cluster_columns = dendsort(hclust(dist(t(gas.HM)))),
                     col = colour.hm,
                     heatmap_legend_param = list(direction = "horizontal"),
                     top_annotation = HeatmapAnnotation(Population =gas.DR$Population,
                                                        Datasets = colnames(gas.HM),
                                                        border = TRUE,gap = unit(0, "points"),
                                                        col = list(Population =c("Mt1-BFP+" = "dodgerblue",
                                                                                 "Sox17-RFP+" = "#E03E3E",
                                                                                 "Bra-GFP+" = "#53BB79",
                                                                                 "Triple-Neg" = "#EF8228"),
                                                                   Datasets = c("D0370_100_1"="#D88C4C",
                                                                                "D0370_100_2"="#E49D41",
                                                                                "D0370_100_3"="#FEB447",
                                                                                "D0386_100_1"="#06CDE0",
                                                                                "D0386_100_2"="#1498BE")),
                                                        show_legend = FALSE),
                     heatmap_width = unit(8, "cm"), 
                     heatmap_height = unit(18, "cm"),
                     column_names_gp = grid::gpar(fontsize = 3),
                     row_names_gp = grid::gpar(fontsize = 5),
                     column_title = paste0("Top ", length(PCAdtmn)," proteins contribute\n to PC1 and PC2"))

gas.HM.hc

png(paste0("Figure/II.5.1. ",expdeg," - Heatmap on Top ",length(PCAdtmn)," PCA determiants - bottom.png"),width=7.25,height=8.25,units="in",res=300)
draw(gas.HM.hc, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
draw(pL, x = unit(0.835, "npc"), y = unit(0.775, "npc"), just = c("right", "top"))
dev.off()

### k-means clustering ----------

n <- 4
kclus <- kmeans(gas.HM, n)
#Heatmap of 991 significant hits
split2 <- paste0("Cluster\n", kclus$cluster)
gas.HM.km <- Heatmap(gas.HM, name = "Scale Log \nAbundance",
                         row_dend_reorder = TRUE,
                         cluster_row_slices = FALSE,
                         split=split2, 
                         clustering_method_rows = "ward.D2",
                         cluster_columns = dendsort(hclust(dist(t(gas.HM)))),
                         col = colour.hm,
                         top_annotation = HeatmapAnnotation(Population =gas.DR$Population,
                                                            Datasets = colnames(gas.HM),
                                                            border = TRUE,gap = unit(0, "points"),
                                                            col = list(Population =c("Mt1-BFP+" = "dodgerblue",
                                                                                 "Sox17-RFP+" = "#E03E3E",
                                                                                 "Bra-GFP+" = "#53BB79",
                                                                                 "Triple-Neg" = "#EF8228"),
                                                                       Datasets = c("D0370_100_1"="#D88C4C",
                                                                                    "D0370_100_2"="#E49D41",
                                                                                    "D0370_100_3"="#FEB447",
                                                                                    "D0386_100_1"="#06CDE0",
                                                                                    "D0386_100_2"="#1498BE")),
                                                            show_legend = FALSE),
                         column_names_gp = grid::gpar(fontsize = 3),
                         row_names_gp = grid::gpar(fontsize = 5),
                         column_title = paste0("Top ", length(PCAdtmn)," proteins contribute\n to PC1 and PC2"))


gas.HM.km
png(paste0("Figure/II.5.2. ",expdeg," - Heatmap on Top ",length(PCAdtmn)," PCA determiants by k-means.png"),
    width=5.25,height=8.25,units="in",res=300)
gas.HM.km
draw(pL, x = unit(1, "npc"), y = unit(0.785, "npc"), just = c("right", "top"))
dev.off()

#order <- row_order(gas.HM.km)
#genelist <- as.numeric(unlist(row_order(gas.HM.km)))


# III. DOWNSTREAM ANALYSES ----------

## 1. Detecting Gene Markers ----------

#Find Gene Marker with scran

#Test <- c("t", "wilcox", "binom")
#pval.type <- c("any", "some", "all")

#t-test
markers <- findMarkers(gas.DR,
                       assay.type = "assay",
                       groups=gas.DR$Population,
                       pval.type="any",
                       lfc=1)


#Wilcoxon rank sum test -> seem to work fine
#markers.WW <- findMarkers(gas.DR,
#                          assay.type = "scaleassay",
#                          groups=gas.DR$Population,
#                          pval.type="any",
#                          test="wilcox",
#                          lfc=0.0125)


### 1.1. Heatmap for Marker proteins/genes ----------

#Layer
layer <-  c("Mt1-BFP+","Bra-GFP+","Triple-Neg","Sox17-RFP+")

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

png(paste0("Figure/III.1.1. ",expdeg," - Heatmap on top ", ntop," significant Log2 Fold Change of Mt1-BFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerMt1-BFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. ",expdeg," - Heatmap on top ", ntop," significant Log2 Fold Change of Sox17-RFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerSox17-RFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. ",expdeg," - Heatmap on top ", ntop," significant Log2 Fold Change of Bra-GFP+ - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerBra-GFP+`, heatmap_legend_side = "bottom")
dev.off()

png(paste0("Figure/III.1.1. ",expdeg," - Heatmap on top ", ntop," significant Log2 Fold Change of Triple-Neg - t-test.png"),width = 5.25,height = 8.25,units="in",res=320)
draw(`MarkerTriple-Neg`, heatmap_legend_side = "bottom")
dev.off()

#Wilcoxon rank sum test -> seem to work fine

#for (i in layer){
#  layer.marker <- markers.WW[[i]]
#  top.marker <- layer.marker[layer.marker$Top <= ntop,]
#  AUCs <- getMarkerEffects(top.marker, prefix="AUC")
#  plot <- Heatmap(AUCs,
#                  col = col.hm,
#                  column_title = paste0("Top ", ntop, " differential expressed proteins of ",i," \nover other populations -WW"),
#                  heatmap_legend_param = list(
#                    #legend_height = unit(4, "cm")))
#                    #title_position = "lefttop-rot"),
#                    direction = "horizontal",
#                    title = "Log2FC"),
#                  show_heatmap_legend = T)
#  assign(paste0("MarkerWW",i),plot)
#}

#MarkerWWMt1-BFP+
#MarkerWWSox17-RFP+
#MarkerWWBra-GFP+
#MarkerWWTriple-Neg


### 1.2. Detect all the Marker proteins/genes ----------

#List all marker gene
genemarker <- list()
for(i in layer){
  layer.marker <- markers[[i]]
  top.marker <- layer.marker[layer.marker$Top <= ntop,]
  genemarker[[i]] <- rownames(top.marker[,])
}

#Marker gene list 
EMK <- which(match(rownames(gas.DR), genemarker %>% unlist)>0)
#Remove duplicated genes
#EMK <- EMK[-which(duplicated(rownames(gas.DR)[EMK]))]

### 1.3. Plot all Marker proteins/genes ----------

#### GeneMarker plot ----------

#non-Scale
GMplot <- plotExpression(gas.DR, 
                         rownames(gas.DR)[EMK],
                         exprs_values = "assay",
                         log2_values = FALSE, #already log2
                         one_facet=TRUE,
                         colour_by = "Population",)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ", ntop," marker per cell population"))
GMplot
ggsave(paste0("Figure/III.1.3.1. ",expdeg," - Expression of top ", ntop, " Marker genes in each cell population - non-Scale.png"),
       plot=GMplot,width=7.25,height=4.25,units="in",dpi=320)

#Scale
GMscaleplot <- plotExpression(gas.DR, 
                         rownames(gas.DR)[EMK],
                         exprs_values = "scaleassay",
                         log2_values = FALSE, #already log2
                         one_facet=TRUE,
                         colour_by = "Population")+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  ggtitle(paste0("Expression of top ", ntop," marker per cell population - Scale"))
GMscaleplot
ggsave(paste0("Figure/III.1.3.2. ",expdeg," - Expression of top ", ntop, " Marker genes in each cell population - Scale.png"), 
       plot=GMscaleplot, width=7.25,height=4.25,units="in",dpi=320)

#### Dot plot ----------

#non-Scale
dotEMK <- plotDots(gas.DR, 
                   features=rownames(gas.DR)[EMK], 
                   exprs_values = "assay",
                   group="Population") + 
  ggtitle(paste0("Top ", ntop," markers per cell population"))+
  theme_bw()
dotEMK

ggsave(paste0("Figure/III.1.3.3. ",expdeg," - Dot expression of top ", ntop, " Marker genes in each cell population - non-Scale.png"),
       plot=dotEMK, width=5,height=7.50,units="in",dpi=320)
       
#Scale
dotscaleEMK <- plotDots(gas.DR, 
         features=rownames(gas.DR)[EMK], 
         exprs_values = "scaleassay",
         group="Population") + 
  ggtitle(paste0("Top ", ntop," markers per cell population - Scale"))+
  theme_bw()
dotscaleEMK

ggsave(paste0("Figure/III.1.3.4. ",expdeg," - Dot expression of top ", ntop, " Marker genes in each cell population - Scale.png"), 
       plot=dotscaleEMK, width=5,height=7.50,units="in",dpi=320)

## 2. Testing clustering with SC3 ----------

gas.SC3 <- gas.DR

#Re-names in to logcoutns and counts for SC3 
#Scale
assay(gas.SC3, "logcounts") <- assay(gas.SC3, "scaleassay") #or non-Scale
assay(gas.SC3, "counts") <- assay(gas.SC3, "aggcounts")
rowData(gas.SC3)$feature_symbol <- rownames(gas.SC3)
gas.SC3 <- gas.SC3[!duplicated(rowData(gas.SC3)$feature_symbol), ]

#Computeing SC3
gas.SC3 <- sc3(gas.SC3, ks = 2:4, biology = TRUE, gene_filter = FALSE)

##Plot Markers
gas.SC3plot <- sc3_plot_markers(gas.SC3, 
                                       k = 3,
                                       auroc = 0.85,
                                       p.val = 0.01,
                                       show_pdata = c("Population",
                                                      "sc3_3_clusters"))
gas.SC3plot

ggsave(paste0("Figure/III.2. ",expdeg," - Marker plot usng SC3 - Scale.png"),
       plot=gas.SC3plot,width=8.25,height=7.25,units="in",dpi=320)

## 3. Enrichment analysis ----------

#LFC > 0.0125 for retaining more inforeation
markers <- findMarkers(gas.DR,
                       assay.type = "assay",
                       groups=gas.DR$Population,
                       pval.type="any",
                       lfc=0.0125)


#Detect DE genes with t-test
#Already compute above
#markers

##layer = ("Mt1-BFP+","Bra-GFP+","Triple-Neg","Sox17-RFP+")

### 3.1. On Mt1-BFP+ population ----------

layer <- "Mt1-BFP+"
DElayer <- as.data.frame(subset(markers[[layer]],FDR <0.05 ))
dim(DElayer)
#Choose LFC
sig_lfc <- 0.5
#Gene up
DEup <- rownames(DElayer[DElayer$summary.logFC>sig_lfc,])
length(DEup)
#Gene down
DEdown <- rownames(DElayer[DElayer$summary.logFC<(-sig_lfc),])
length(DEdown)

###
genes <- keys(org.Mm.eg.db, keytype="SYMBOL")
geneUniverse <- AnnotationDbi::select(org.Mm.eg.db, 
                                      keys = genes,
                                      columns=c("ENTREZID"), 
                                      keytype = "SYMBOL")

#### 1. Gene Ontology ----------
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2

#MF: Molecular Function: molecular activities of gene products 
goDEMF <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)

#CC: Cellular Component where gene products are active 
goDECC <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)

#BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
goDEBP <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)

#Visualization of Functional Enrichment Result

#1.Bar Plot
goDEMFplot <- dotplot(goDEMF, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
goDECCplot <- dotplot(goDECC, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
goDEBPplot <- dotplot(goDEBP, showCategory=30)+ 
  ggtitle(paste0("UpGene in ",layer, "\n Biological Process - FDR < ",pvalueCutoff))
title <- "Over Representation Analysis of Mt1-BFP+"
DElayerGO <- plot_grid(goDEMFplot, 
                       #goDECCplot,
                       goDEBPplot,
                       labels = title, 
                       ncol=2)
DElayerGO
ggsave(paste0("Figure/III.3.1. ",expdeg," - GO term of ",layer,".png"),
       plot = DElayerGO,width=20,height=7.25,units="in",dpi=320)

#### 2. Gene Set Enrichment Analysis ----------

#Gene List, include up and down

#Gene up
DELup <- DElayer[DElayer$summary.logFC>sig_lfc,]

#Gene down
DELdown <- DElayer[DElayer$summary.logFC<(-sig_lfc),]

#Gene List
DElist <- rbind(DELup, DELdown)
dim(DElist)
DElistFC = DElist[,"summary.logFC"]
names(DElistFC) = as.character(rownames(DElist))
DElist <- sort(DElistFC, decreasing = TRUE)

###
pvalueCutoff <- 0.2

#MF: Molecular Function: molecular activities of gene products
gseDEMF <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "MF",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)

#CC: Cellular Component where gene products are active 
gseDECC <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "CC",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)

##BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
gseDEBP <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "BP",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)


#2. Dot plot
gseDEMFplot <- dotplot(gseDEMF, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
gseDECCplot <- dotplot(gseDECC, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
title <- "Gene Set Enrichment Analysis of Mt1-BFP+"
gseDEplot <- plot_grid(gseDEMFplot, 
                       #gseDECCplot,
                       #No BP for Mt1-BFP+
                       labels = title, 
                       ncol=1)

ggsave(paste0("Figure/III.3.2. ",expdeg," - GSEA term of ",layer,".png"),
       plot = gseDEMFplot,width=20,height=7.25,units="in",dpi=320)


#### 3. Gene-Concept Network ----------

#or colorEdge = TRUE
cnegoDEMF <- cnetplot(goDEMF,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
cnegoDEBP <- cnetplot(goDEBP,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
title <- "Gene-Concept Network - MF and BP"
cneplot <- cowplot::plot_grid(cnegoDEMF,
                              cnegoDEBP,
                              labels = title,
                              ncol=2)

ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer,".png"), plot = cneplot,
       width=19,height=7.25,units="in",dpi=320)

ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," MF.png"), 
       plot = cnegoDEMF + ggtitle(paste0(layer," - Molecular Function")),width=12,height=7,units="in",dpi=320)
ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," BP.png"), 
       plot = cnegoDEBP + ggtitle(paste0(layer," - Biological Process")),width=12,height=7,units="in",dpi=320)

#### 4.Heatmap-like functional classification
heatgoDEMFplot <- heatplot(goDEMF, 
                       foldChange=DElist)+
  ggtitle("Heatmap-like functional classification")

png(paste0("Figure/III.3.4. ",expdeg," - Heatmap-like functional classification of ",layer,".png"),width=12,height=7.25,units="in",res=300)
heatgoDEMFplot
dev.off()

#### 5. Enrichment Map

DElayerpw <- pairwise_termsim(goDEMF)
#layout	
#Layout:'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 
#'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
emgoDEMFplot <- emapplot(DElayerpw,
                     pie_scale=0.5,
                     layout="kk", 
                     cex_category=1.5)+
  ggtitle("Enrichment Map")
png(paste0("Figure/III.3.5. ",expdeg," - Enrichment Map of",layer,".png"),width=12,height=7.25,units="in",res=300)
emgoDEMFplot
dev.off()


#### 6. UpSet Plot
upsetgoDEMF <- upsetplot(goDEMF)+
  ggtitle("Over Representation Analysis")
upsetgseDEMF <- upsetplot(gseDEMF)+
  ggtitle("Gene Set Enrichment Analysis")

title <- "UpSet Plot"
upsetgoDEMFplot <- plot_grid(upsetgoDEMF,
                             upsetgseDEMF,
                             labels = title,
                             nrow=2)

png(paste0("Figure/III.3.6. ",expdeg," - UpSet Plot of",layer,".png"),width=12,height=7.25,units="in",res=300)
upsetgoDEMFplot
dev.off()

#### 7. ridgeline plot for expression distribution of GSEA result
ridgegseDEMF <- ridgeplot(gseDEMF)+
  ggtitle("Ridgeline Plot Analysis")

png(paste0("Figure/III.3.7. ",expdeg," - Ridgeline Plot of",layer,".png"),width=12,height=7.25,units="in",res=300)
ridgegseDEMF
dev.off()


### 3.2. On Sox17-RFP+ population ----------

layer <- "Sox17-RFP+"
DElayer <- as.data.frame(subset(markers[[layer]],FDR <0.05 ))
dim(DElayer)
#Choose LFC
sig_lfc <- 0.5
#Gene up
DEup <- rownames(DElayer[DElayer$summary.logFC>sig_lfc,])
length(DEup)
#Gene down
DEdown <- rownames(DElayer[DElayer$summary.logFC<(-sig_lfc),])
length(DEdown)

###
genes <- keys(org.Mm.eg.db, keytype="SYMBOL")
geneUniverse <- AnnotationDbi::select(org.Mm.eg.db, 
                                      keys = genes,
                                      columns=c("ENTREZID"), 
                                      keytype = "SYMBOL")
#### 1. Gene Ontology ----------
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2
#MF: Molecular Function: molecular activities of gene products 
goDEMF <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#CC: Cellular Component where gene products are active 
goDECC <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
goDEBP <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#Visualization of Functional Enrichment Result
#1.Bar Plot
barplot(goDEMF, showCategory=20)
#2. Dot plot
goDEMFplot <- dotplot(goDEMF, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
goDECCplot <- dotplot(goDECC, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
goDEBPplot <- dotplot(goDEBP, showCategory=30)+ 
  ggtitle(paste0("UpGene in ",layer, "\n Biological Process - FDR < ",pvalueCutoff))
title <- "Over Representation Analysis"
DElayerGO <- plot_grid(goDEMFplot, 
                       #goDECCplot,
                       goDEBPplot,
                       labels = title, 
                       ncol=2)

png(paste0("Figure/III.3.1. ",expdeg," - GO term of ",layer,".png"),width=20,height=7.25,units="in",res=300)
DElayerGO
dev.off()

#### 2. Gene Set Enrichment Analysis ----------

#Gene List, include up and down
#Gene up
DELup <- DElayer[DElayer$summary.logFC>sig_lfc,]
#Gene down
DELdown <- DElayer[DElayer$summary.logFC<(-sig_lfc),]
#Gene List
DElist <- rbind(DELup, DELdown)
dim(DElist)
DElistFC = DElist[,"summary.logFC"]
names(DElistFC) = as.character(rownames(DElist))
DElist <- sort(DElistFC, decreasing = TRUE)
#
###
pvalueCutoff <- 0.2
#MF: Molecular Function: molecular activities of gene products
gseDEMF <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "MF",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
#CC: Cellular Component where gene products are active 
gseDECC <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "CC",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
##BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
gseDEBP <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "BP",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
#2. Dot plot
gseDEMFplot <- dotplot(gseDEMF, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
gseDECCplot <- dotplot(gseDECC, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
gseDEBPplot <- dotplot(gseDEBP, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Biological Component - FDR < ",pvalueCutoff))
title <- "Gene Set Enrichment Analysis"
gseDEplot <- plot_grid(gseDEMFplot, 
                       #gseDECCplot,
                       #gseDEBPplot,
                       #No BP for Mt1-BFP+
                       labels = title, 
                       ncol=2)
png(paste0("Figure/III.3.2. ",expdeg," - GSEA term of ",layer,".png"),width=16,height=7.25,units="in",res=300)
gseDEplot
dev.off()

#### 3. Gene-Concept Network ----------
#or colorEdge = TRUE
cnegoDEMF <- cnetplot(goDEMF,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
cnegoDEBP <- cnetplot(goDEBP,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
title <- "Gene-Concept Network"
cneplot <- cowplot::plot_grid(cnegoDEMF,
                              cnegoDEBP,
                              labels = title,
                              ncol=2)

png(paste0("Figure/III.3.3. ",expdeg," - Gene-Concept Network of ",layer,".png"),width=19,height=7.25,units="in",res=300)
cneplot
dev.off()

ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," MF.png"), 
       plot = cnegoDEMF + ggtitle(paste0(layer," - Molecular Function")),width=12,height=7,units="in",dpi=320)
ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," BP.png"), 
       plot = cnegoDEBP + ggtitle(paste0(layer," - Biological Process")),width=12,height=7,units="in",dpi=320)


### 3.3. On Bra-GFP+ population ----------

layer <- "Bra-GFP+"
DElayer <- as.data.frame(subset(markers[[layer]],FDR <0.05 ))
dim(DElayer)
#Choose LFC
sig_lfc <- 0.5
#Gene up
DEup <- rownames(DElayer[DElayer$summary.logFC>sig_lfc,])
length(DEup)
#Gene down
DEdown <- rownames(DElayer[DElayer$summary.logFC<(-sig_lfc),])
length(DEdown)



#### 1. Gene Ontology ----------
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2


#MF: Molecular Function: molecular activities of gene products 
goDEMF <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#CC: Cellular Component where gene products are active 
goDECC <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
goDEBP <- enrichGO(DEup,
                   org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = "fdr",
                   qvalueCutoff = qvalueCutoff)
#Visualization of Functional Enrichment Result
#1.Bar Plot
barplot(goDEMF, showCategory=20)
#2. Dot plot
goDEMFplot <- dotplot(goDEMF, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
goDECCplot <- dotplot(goDECC, showCategory=30) + 
  ggtitle(paste0("UpGene in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
goDEBPplot <- dotplot(goDEBP, showCategory=30)+ 
  ggtitle(paste0("UpGene in ",layer, "\n Biological Process - FDR < ",pvalueCutoff))
title <- "Over Representation Analysis"
DElayerGO <- plot_grid(goDEMFplot, 
                       #goDECCplot,
                       goDEBPplot,
                       labels = title, 
                       ncol=2)

png(paste0("Figure/III.3.1. ",expdeg," - GO term of ",layer,".png"),width=20,height=7.25,units="in",res=300)
DElayerGO
dev.off()

#### 2. Gene Set Enrichment Analysis ----------

##GSEA  
#Gene List, include up and down
#Gene up
DELup <- DElayer[DElayer$summary.logFC>sig_lfc,]
#Gene down
DELdown <- DElayer[DElayer$summary.logFC<(-sig_lfc),]
#Gene List
DElist <- rbind(DELup, DELdown)
dim(DElist)
DElistFC = DElist[,"summary.logFC"]
names(DElistFC) = as.character(rownames(DElist))
DElist <- sort(DElistFC, decreasing = TRUE)
#
###
pvalueCutoff <- 0.2
#MF: Molecular Function: molecular activities of gene products
gseDEMF <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "MF",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
#CC: Cellular Component where gene products are active 
gseDECC <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "CC",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
##BP: Biological Process pathways and larger processes made up of 
#the activities of multiple gene products
gseDEBP <- gseGO(geneList = DElist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont          = "BP",
                 pvalueCutoff = pvalueCutoff,
                 pAdjustMethod = "fdr",
                 verbose      = TRUE)
#2. Dot plot
gseDEMFplot <- dotplot(gseDEMF, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Molecular Function - FDR < ",pvalueCutoff))
gseDECCplot <- dotplot(gseDECC, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Cellular Component - FDR < ",pvalueCutoff))
gseDEBPplot <- dotplot(gseDEBP, showCategory=5) +
  ggtitle(paste0("\n DElayer in ",layer, "\n Biological Process - FDR < ",pvalueCutoff))

title <- "Gene Set \n Enrichment Analysis"
gseDEplot <- plot_grid(gseDEMFplot, 
                       #gseDECCplot,
                       
                       labels = title, 
                       ncol=2)
png(paste0("Figure/III.3.2. ",expdeg," - GSEA term of ",layer,".png"),width=8,height=7.25,units="in",res=300)
gseDEMFplot
dev.off()

#### 3. Gene-Concept Network ----------
#or colorEdge = TRUE
cnegoDEMF <- cnetplot(goDEMF,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
cnegoDEBP <- cnetplot(goDEBP,
                      showCategory = 5,
                      foldChange=DElist,
                      layout = "kk",
                      categorySize="pvalue",
                      node_label="all")
title <- "Gene-Concept Network"
cneplot <- cowplot::plot_grid(cnegoDEMF,
                              cnegoDEBP,
                              labels = title,
                              ncol=2)

png(paste0("Figure/III.3.3. ",expdeg," - Gene-Concept Network of ",layer,".png"),width=19,height=7.25,units="in",res=300)
cneplot
dev.off()

ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," MF.png"), 
       plot = cnegoDEMF + ggtitle(paste0(layer," - Molecular Function")),width=12,height=7,units="in",dpi=320)
ggsave(paste0("Figure/III.3.3. ",expdeg," Gene-Concept Network of ",layer," BP.png"), 
       plot = cnegoDEBP + ggtitle(paste0(layer," - Biological Process")),width=12,height=7,units="in",dpi=320)


# IV. Integrate with scRNA-seq dataset ------

# Load lists from four cluster in 120h mouse gastruloids 
DEG_endoderm <- read.delim("Gas.100Data/DEG_endoderm.txt", header = F)
DEG_ectoderm <- read.delim("Gas.100Data/DEG_ectoderm.txt", header = F)
DEG_PSM <- read.delim("Gas.100Data/DEG_PSM.txt", header = F)
DEG_NMP <- read.delim("Gas.100Data/DEG_NMPs.txt", header = F)

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



## Endoderm
tmp12 <- which(match(rownames(gas.DR), DEG_endoderm$gene)>0)
length(tmp12)
tmp12a <- which(match(DEG_endoderm$gene, rownames(gas.DR))>0)
length(tmp12a)
dim(DEG_endoderm[tmp12a, ])

## Ectoderm
tmp13 <- which(match(rownames(gas.DR), DEG_ectoderm$gene)>0)
length(tmp13)
tmp13a <- which(match(DEG_ectoderm$gene, rownames(gas.DR))>0)
length(tmp13a)
dim(DEG_ectoderm[tmp13a, ])

## PSM
tmp6 <- which(match(rownames(gas.DR), DEG_PSM$gene)>0)
length(tmp6)
tmp6a <- which(match(DEG_PSM$gene, rownames(gas.DR))>0)
length(tmp6a)
dim(DEG_PSM[tmp6a, ])

## NMPs
tmp7 <- which(match(rownames(gas.DR), DEG_NMP$gene)>0)
length(tmp7)
tmp7a <- which(match(DEG_NMP$gene, rownames(gas.DR))>0)
length(tmp7a)
dim(DEG_NMP[tmp7a, ])

#Overlapped genes in each cluster

##With gasBulk
gas.DR.HM.DEG <- rbind(assay(gas.DR,"scaleassay")[tmp12,],
                    assay(gas.DR,"scaleassay")[tmp13,],
                    assay(gas.DR,"scaleassay")[tmp6,],
                    assay(gas.DR,"scaleassay")[tmp7,])
dim(gas.DR.HM.DEG)

##With van der Brink data
Brink.DEG <- rbind(DEG_endoderm[tmp12a,],
                   DEG_ectoderm[tmp13a,],
                   DEG_PSM[tmp6a,],
                   DEG_NMP[tmp7a,]) 
dim(Brink.DEG)

#Shorten col names
colnames(gas.DR.HM.DEG) <- sub('D0370_100_1Reporter.intensity.*','D0370_100_1',colnames(gas.DR.HM.DEG))
colnames(gas.DR.HM.DEG) <- sub('D0370_100_2Reporter.intensity.*','D0370_100_2',colnames(gas.DR.HM.DEG))
colnames(gas.DR.HM.DEG) <- sub('D0370_100_3Reporter.intensity.*','D0370_100_3',colnames(gas.DR.HM.DEG))
colnames(gas.DR.HM.DEG) <- sub('D0386_100_1Reporter.intensity.*','D0386_100_1',colnames(gas.DR.HM.DEG))
colnames(gas.DR.HM.DEG) <- sub('D0386_100_2Reporter.intensity.*','D0386_100_2',colnames(gas.DR.HM.DEG))

#Create Heatmap for top PC determinants
#Manually create legend of Cells Population and Datasets
lgChan <- Legend(labels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg"), 
                 legend_gp = gpar(fill = colour.Chan), 
                 title = "Population")
lgData <- Legend(labels = c("0370_1","0370_2","0370_3","0386_1","0386_2"),
                 legend_gp = gpar(fill = c("#D88C4C","#E49D41","#FEB447","#06CDE0","#1498BE")),
                 title = "Datasets")
pL = packLegend(list = list(lgChan, lgData))


#Choose overlapped genes (for adding into Heatmap)
overlapped.genes <- rownames(gas.DR.HM.DEG)

#### Hierarchical clustering -----
#dist(1-cor(t(gasBulk.aov.hm))) 
#"pearson"
#spearman
#kendall
#euclidean

####

gas.DR.HM.DEG.hc <- Heatmap(gas.DR.HM.DEG, name = "Scale Log Abundance",
                            clustering_distance_rows = "pearson",
                            clustering_method_rows = "ward.D2",
                            cluster_columns = dendsort(hclust(dist(t(gas.DR.HM.DEG)))),
                            col = colour.hm,
                            heatmap_legend_param = list(direction = "horizontal"),
                            top_annotation = HeatmapAnnotation(Population =gas.DR$Population,
                                                               Datasets = colnames(gas.DR.HM.DEG),
                                                               border = TRUE,gap = unit(0, "points"),
                                                               col = list(Population =c("Mt1-BFP+" = "dodgerblue",
                                                                                    "Sox17-RFP+" = "#E03E3E",
                                                                                    "Bra-GFP+" = "#53BB79",
                                                                                    "Triple-Neg" = "#EF8228"),
                                                                          Datasets = c("D0370_100_1"="#D88C4C",
                                                                                       "D0370_100_2"="#E49D41",
                                                                                       "D0370_100_3"="#FEB447",
                                                                                       "D0386_100_1"="#06CDE0",
                                                                                       "D0386_100_2"="#1498BE")),
                                                                          show_legend = FALSE),
                            right_annotation = rowAnnotation(Brink = Brink.DEG$cluster,
                                                             border = TRUE,
                                                             col = list(Brink = c("Ecto" = "dodgerblue2",
                                                                                  "PSM" = "#B7D468",
                                                                                  "NMP" = "#2EC20A",
                                                                                  "Endo" = "#D14239"))),
                            column_title = paste0("Overlapped ",nrow(gas.DR.HM.DEG)," protein/mRNA \nwith van den Brink et al."),
                            column_names_gp = grid::gpar(fontsize = 6),
                            row_names_gp = grid::gpar(fontsize = 1)) +
  rowAnnotation(gene = anno_text(overlapped.genes, which = "row",gp = gpar(fontsize = 4)),
                width = max_text_width(unlist(overlapped.genes)) + unit(0.1, "cm"))
gas.DR.HM.DEG.hc

png(paste0("Figure/IV.2. Heatmap on Overlapped ",nrow(gas.DR.HM.DEG)," protein-mRNA with van den Brink et al - unique genes.png"), width = 5.25, height = 8.25, units="in", res = 320)
draw(gas.DR.HM.DEG.hc, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
#draw(pL, x = unit(1, "npc"), y = unit(0.845, "npc"), just = c("right", "top"))
dev.off()









#### k-means clustering -----
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
                                   top_annotation = HeatmapAnnotation(Population =Population,
                                                                      border = TRUE,
                                                                      #gap = unit(2, "points"),
                                                                      col = list(Population =c("Mt1-BFP+" = "dodgerblue",
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
draw(lgChan, x = unit(0.935, "npc"), y = unit(0.715, "npc"), just = c("right", "top"))
dev.off()




#V. Monitoring data processing ----------

#Some top PCA determinants for plotting
gas.100[PCAdtmn,"Leading.razor.protein"]
###Good test: "Q8K3F7"

#Mt1-BFP+: O70400
#Sox17-RFP+: P10922
#Bra-GFP+: Q9JMG7 

proID <- "Q9JMG7"

proIDplot <- scp %>%
  ## Get the features related to proID
  subsetByFeature(proID) %>%
  ## Format the `QFeatures` to a long format table
  longFormat(colDataCols = c("Raw.file", "Population", "Channel")) %>%
  data.frame %>%
  ## This is used to preserve ordering of the samples and assays in ggplot2
  mutate(assay = factor(assay, levels = names(scp)),
         Channel = sub("Reporter.intensity.", "TMT-", Channel),
         Channel = factor(Channel, levels = unique(Channel))) %>%
  ## Start plotting
  ggplot(aes(x = Channel, y = value, group = rowname, col = Population)) +
  geom_point() +
  scale_colour_manual(values = colour.Chan) +
  ## Plot every assay in a separate facet
  facet_wrap(facets = vars(assay), scales = "free_y", ncol = round(length(scp)/5)) +
  ## Annotate plot
  xlab("Channels") +
  ylab("Intensity (arbitrary units)") +
  ## Improve plot aspect
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        strip.text = element_text(hjust = 0),
        legend.position = "bottom")+
  ggtitle(paste0("Protein: ",proID))
proIDplot

png(paste0("Figure/IV. Testing Protein-",proID,".png"),width=10,height=6,units="in",res=300)
proIDplot
dev.off()