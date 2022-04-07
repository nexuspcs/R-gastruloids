packages <- c(
  ##Data manipulation and analytics
  "magrittr","dplyr","knitr","tidyverse","scales",
  "scp","scater","scran","sva","SC3",
  #Visualization
  "RColorBrewer","dendsort","cowplot", "ggplot2","ComplexHeatmap","circlize")
  #Downstream analysis
  #"enrichplot","ggupset","clusterProfiler",
  #Mouse database
  #"org.Mm.eg.db")
  
lapply(packages, require, character.only = TRUE)

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("Gas.SCData/gas.SC.ProcessedData.Rda")

# I. DATA PROCESSING ----------

## Input ----------

###Input data

#Process MQ files from a MQ run
gas.SC <- read.delim("Gas.SCData/4&5. evidence - Single Cells - LMNOIJK.txt")
expdeg <- "Single Cell"

#Cleaning
#Remove >2nd leading proteins
gas.SC$Leading.razor.protein <- sub(';.*','',gas.SC$Leading.razor.protein)
#Remove >2nd leading genes
gas.SC$Gene.names <- sub(';.*','',gas.SC$Gene.names)

#Combines 2 CVs into one (No need if set names for "Experiments" columns):
#gas.SC$Raw.file <- sub('_FS45','',gas.SC$Raw.file)
#gas.SC$Raw.file <- sub('_FS65','',gas.SC$Raw.file)

#Annotation file
sampleAnnotation <- read.csv("Gas.SCData/4&5. sampleAnnotation - Gas_SC_LMNOJJK.csv")

#Colour

#For heatmap
colour.hm <-  c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
#For Channels
colour.Chan <- c("dodgerblue","#E03E3E","#53BB79","#EF8228","#937264","#3043A2","#C25D7B")
#For Datasets
colour.SC <- c("#776AD6","#A068B0","#E487B8",
               "#BBC385","#7FAD81","#0F7B6C","#0D5B11")


## I.1. PSM level ----------

### 1. Reconstructing the dataset to create SCP/QFeatures object ----------
#read SCP
scp <- readSCP(featureData = gas.SC,
               colData = sampleAnnotation,
               channelCol = "Channel",
               batchCol = "Experiment",
               #removeEmptyCols = TRUE
               )

#Number of assays
l <- length(scp)

### 2. Quality control at PSM level ----------

#### 2.1. Cleaning missing data ----------

#The zeros can be biological zeros or technical zeros.
#Therefore, any zero should be replaced by NA to avoid artefacts in downstream analysis.
scp <- zeroIsNA(scp, i = names(scp))

#### 2.2. Filter out failed runs based on PSM content (detected features) ----------

#Plot for checking number of PSM per assays
#Number of PSM per datasets 
#nPSM <- dims(scp)[1, which(startsWith(names(scp),"sc"))]
nPSM <- dims(scp)[1, ]

#Create data frame
nPSM <- data.frame(Datasets = names(nPSM),
                       Count = nPSM)
nPSM$Datasets <- factor(nPSM$Datasets, levels = nPSM$Datasets)

#Plot
nPSMplot <- ggplot(data=nPSM, aes(x=Datasets, y=Count, fill=Datasets)) +
  scale_fill_manual(values= colour.SC, 
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
       plot=nPSMplot, width=5, height=4, dpi=320)


#Select the assays that have sufficient PSMs
#e.g. the number of rows is greater than 150),
keepAssay <- dims(scp)[1, ] > 500
scp <- scp[, , keepAssay]

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

#pep2qvalue function to convert PEPs to q-values (DART-ID algorithm) that are directly related to FDR
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

#Gilter the PSM to control, the protein FDR at 1%
scp <- filterFeatures(scp,
                      ~ qvalue_proteins < 0.01)

#or
#scp <- filterFeatures(scp,
#                      ~ qvalue_psm < 0.01 & qvalue_protein < 0.01)

#Save
scp2.4 <- scp

#### 2.5. Filter out PSMs with high sample to carrier ratio (for single-cell) ----------

scp <- computeSCR(scp, 
                  i = names(scp),
                  colDataCol = "Population",
                  carrierPattern = "Carrier",
                  samplePattern = "Mt1-BFP+|Sox17-RFP+|Bra-GFP+|Triple-Neg",
                  rowDataName = "MeanSCR")


scp2.5 <- scp
#Mean for MeanSCR
MeanSCR.df <- scp %>%  rowDataToDF(i = names(scp), 
                                 vars = "MeanSCR") %>%
  data.frame 


#Plot
SCRplot <- ggplot(MeanSCR.df, aes(x = MeanSCR)) +
  geom_histogram(binwidth=0.015) +
  geom_vline(xintercept = c(1/100, median(MeanSCR.df$MeanSCR, na.rm = TRUE), 0.1),
             lty = c(2, 2, 1),
             linetype='dashed', 
             col = c('red', 'blue','black')) +
  scale_x_log10()+
  xlab("Mean Single cells to Carrier ratio")+
  ylab("Count")+
  theme_bw()
SCRplot

ggsave(paste0("Figure/2.5.1. ",expdeg," - MeanSCR .png"), 
       plot=SCRplot, width=7.2,height=4.45, dpi=320)

#Filter based on MeanSCR
scp <- filterFeatures(scp,
                      ~ !is.na(MeanSCR) & MeanSCR < 0.1)

#Mean for MeanSCR
MeanSCR.df <- scp %>%  rowDataToDF(i = names(scp), 
                                   vars = "MeanSCR") %>%
  data.frame 

#Plot QC
SCRqcplot <- ggplot(MeanSCR.df, aes(x = MeanSCR)) +
  geom_histogram(binwidth=0.015) +
  geom_vline(xintercept = c(1/100, median(MeanSCR.df$MeanSCR, na.rm = TRUE), 0.1),
             lty = c(2, 2, 1),
             linetype='dashed',
             col = c('red', 'blue','black')) +
  scale_x_log10()+
  xlab("Mean Single cells to Carrier ratio")+
  ylab("Count")+
  theme_bw()
SCRqcplot


ggsave(paste0("Figure/2.5.2. ",expdeg," - MeanSCR after QC .png"), 
       plot=SCRqcplot, width=7.2,height=4.45, dpi=320) 

#### 2.6. Normalize to reference ----------
scp <- divideByReference(scp,
                         i = names(scp),
                         colDataCol = "Population",
                         samplePattern = ".",
                         refPattern = "Reference")

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
                         levels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg",
                                    "Empty","Reference","Carrier"))

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
  scale_colour_manual(values = colour.SC) +
  #scale_y_log10()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^-2, 10^2))+
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
                       nobs = 3, 
                       norm = "SCoPE2",
                       colDataName = "MedianCV")

#Threshold for minimum of MedianCV for Empty Channel
#ts <- min(scp$MedianCV[scp$Population=="Empty"],na.rm = TRUE) - 0.001 # - 0.001 for drawinga line below the minium value
ts <- 0.4
#Plot
MedianCVplot <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = Population, y = MedianCV, fill = Population)) +
  geom_boxplot(width=0.4,outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Dataset),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.SC)+
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
scp <- scp[, scp$Population %in% c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg")]

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
  scale_colour_manual(values = colour.SC) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^-2, 10^2))+
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
  scale_colour_manual(values = colour.SC)+
  ylim(0,1) +
  geom_hline(yintercept = 0.4, linetype='dashed', col = 'red') +
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
gas.DR <- scp[["Proteins_batchC"]]

#Add a scale-log assays for some heatmap/visualization
assay(gas.DR,"scaleassay") <- t(scale(t(assay(gas.DR,"assay")), 
                                      center = TRUE, scale = TRUE))

#Find proteins that have no correspoding gene
blankgene <- gas.SC[match(rownames(gas.DR), gas.SC$Leading.razor.protein),][,c("Leading.razor.protein","Gene.names")]

#Search on  Uniprot
blankgene[which(blankgene[,"Gene.names"] == ""),]

#No blank gene

#Change rownames gas.DR from proteins names to correspoding gene names
rownames(gas.DR) <- gas.SC[match(rownames(gas.DR), gas.SC$Leading.razor.protein),]$Gene.names

#Check order
#View(data.frame(rownames(scp[["Proteins_batchC"]]),
#                rownames(gas.DR),
#                gas.SC[match(rownames(scp[["Proteins_batchC"]]), gas.SC$Leading.razor.protein),]$Gene.names,d,
#                gas.SC[match(rownames(scp[["Proteins_batchC"]]), gas.SC$Leading.razor.protein),]$Leading.razor.protein))

save(gas.SC,
     expdeg,
     sampleAnnotation,
     nPSM,
     scp2.4,
     scp2.5,
     scp2,
     scp5,
     scp6.1,
     scp6.2,
     scp9,
     scp10,
     scp11,
     scp12,
     scpB,
     scp,
     gas.DR,
     MedianCVplot,
     MedianRIplot,
     colour.hm,
     colour.Chan,
     colour.SC,
     file="gas.SCData/gas.SC.ProcessedData.Rda")

#Save
#gas.DR1 <- gas.DR

# II. DIMENSIONALITY REDUCTION AND VISUALIZATION ----------

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("gas.SCData/gas.SC.ProcessedData.Rda")


## 1. Principal component analysis ----------

#Change order for PCA plot
colour.PCA <- c("#53BB79","dodgerblue","#E03E3E","#EF8228")

#Choose n top proteins, scale/ non-Scale may generate different results
set.seed(100)
PCAtop <- 100
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

ggsave(paste0("Figure/II.1.1. ",expdeg," - PCA for top ",PCAtop," PCA determinants.png"), 
       plot=PCAplot, width=3, height=2, dpi=320)


#1.2. PCA Check Batch
gas.DR$lcbatch <- as.factor(sub('0386_sc_','',gas.DR$lcbatch)) 
PCAplotbatch <- plotReducedDim(gas.DR,
                               dimred = "PCA",
                               colour_by = "lcbatch",
                               point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = c("#776AD6","#0D5B11"))+
  labs(colour = "Batch")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR), "percentVar")[2],1), "%)"))
PCAplot

ggsave(paste0("Figure/II.1.2. ",expdeg," - PCA Batch for top ",PCAtop," PCA determinants.png"), 
       plot=PCAplotbatch, width=2.8, height=2.2, dpi=320)

##1.3. Batch

gas.DR.nB <- getWithColData(scp, "Proteins_norm")

#Choose n top proteins, scale/ non-Scale may generate different results
set.seed(100)
PCAtop <- 100
gas.DR.nB <- runPCA(gas.DR.nB,
                 ncomponents = 10,
                 ntop = PCAtop,
                 scale = FALSE,
                 exprs_values = "assay",
                 name = "PCA")

#Note: If scale=TRUE, the expression values for each feature are standardized 
#so that their variance is unity. This will also remove features with standard deviations below 1e-8.
#1.3. PCA plot
PCAplotNB <- plotReducedDim(gas.DR.nB,
                          dimred = "PCA",
                          colour_by = "Population",
                          point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[2],1), "%)"))
PCAplotNB

ggsave(paste0("Figure/II.1.3. ",expdeg," - PCA for top ",PCAtop," PCA determinants - NB.png"), 
       plot=PCAplotNB, width=3, height=2, dpi=320)


#1.4. PCA Check Batch
gas.DR.nB$lcbatch <- as.factor(sub('0386_sc_','',gas.DR.nB$lcbatch)) 
PCAplotbatchNB <- plotReducedDim(gas.DR.nB,
                               dimred = "PCA",
                               colour_by = "lcbatch",
                               point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = c("#776AD6","#0D5B11"))+
  labs(colour = "Batch")+
  xlab(paste0("PC1 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[1],1), "%)"))+
  ylab(paste0("PC2 (", round(attr(reducedDim(gas.DR.nB), "percentVar")[2],1), "%)"))
PCAplotbatchNB

ggsave(paste0("Figure/II.1.4. ",expdeg," - PCA Batch for top ",PCAtop," PCA determinants - NB.png"), 
       plot=PCAplotbatchNB, width=2.8, height=2.2, dpi=320)


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
  labs(colour = "Population")
PCAtoplot

ggsave(paste0("Figure/II.2.2. ",expdeg," - PCA with more PCs plot.png"), 
       plot=PCAtoplot, width=7,height=6, dpi=320)

## 3. Uniform Manifold Approximation and Projection ----------

#850
set.seed(850)
gas.DR <- runUMAP(gas.DR,
                  ncomponents = elbow,
                  exprs_values = "assay",
                  n_neighbors = 4,
                  dimred = "PCA.elbow",
                  name = "UMAP")
UMAPplot <- plotReducedDim(gas.DR,
                           dimred = "UMAP",
                           colour_by = "Population",
                           point_alpha = 1)+
  ggplot2::scale_color_manual(values = colour.PCA)+
  labs(colour = "Population")
UMAPplot

ggsave(paste0("Figure/II.3. ",expdeg," - UMAPplot.png"), 
              plot=UMAPplot, width=4,height=3,units="in",dpi=320)

## 4. Detect the PCA determinants ----------

## Detect PCA determinants/gene that contribute to the most variance PC - Based on Rotation 
pcrot <- attr(reducedDim(gas.DR),"rotation")

#Number of PCs choosen
npc <- 2
pclist <- list()
for (i in 1:npc) {
  pc <- round(pvplot$pv[i]/sum(pvplot$pv[1:npc])*100)
  genePC <- rownames(pcrot[order(abs(pcrot[,i]), decreasing = T)[1:pc],])
  pclist[[i]] <- genePC
}

#List of protein contribute to PC1 and PC2
listPC <- c(pclist[[1]],pclist[[2]])
length(listPC)

#Find the top PCA determinants
PCAdtmn <- which(match(rownames(gas.DR), listPC)>0)
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
colnames(gas.HM) <- sub('sc_IReporter.intensity.*','sc_I',colnames(gas.HM))
colnames(gas.HM) <- sub('sc_JReporter.intensity.*','sc_J',colnames(gas.HM))
colnames(gas.HM) <- sub('sc_KReporter.intensity.*','sc_K',colnames(gas.HM))
colnames(gas.HM) <- sub('sc_NReporter.intensity.*','sc_N',colnames(gas.HM))
colnames(gas.HM) <- sub('sc_OReporter.intensity.*','sc_O',colnames(gas.HM))


#Create Heatmap for top PC determinants
#Manually create legend of Cells Population and Datasets
lgChan <- Legend(labels = c("Mt1-BFP+","Sox17-RFP+","Bra-GFP+","Triple-Neg"), 
                 legend_gp = gpar(fill = colour.Chan), 
                 title = "Population")
lgData <- Legend(labels = c("sc_I","sc_J","sc_K", "sc_N","sc_O"),
                 legend_gp = gpar(fill = c("#776AD6","#A068B0","#E487B8","#0F7B6C","#0D5B11")),
                 title = "Datasets")
pL = packLegend(list = list(lgChan, lgData))

### Hierarchical clustering ----------

#(Wihout spliting: remove "column_split = 4,")
gas.HM.hc <- Heatmap(gas.HM, name = "Scale Log \nAbundance",
                     clustering_distance_rows = "pearson",
                     clustering_method_rows = "ward.D2",
                     cluster_columns = dendsort(hclust(dist(t(gas.HM)))),
                     col = colour.hm,
                     top_annotation = HeatmapAnnotation(Population = gas.DR$Population,
                                                        Datasets = colnames(gas.HM),
                                                        border = TRUE,gap = unit(0, "points"),
                                                        col = list(Population = c("Mt1-BFP+" = "dodgerblue",
                                                                                  "Sox17-RFP+" = "#E03E3E",
                                                                                  "Bra-GFP+" = "#53BB79",
                                                                                  "Triple-Neg" = "#EF8228"),
                                                                   Datasets = c("sc_I"="#776AD6",
                                                                                "sc_J"="#A068B0",
                                                                                "sc_K"="#E487B8",
                                                                                "sc_O"="#0F7B6C",
                                                                                "sc_N"="#0D5B11")),
                                                        show_legend = FALSE),
                     heatmap_width = unit(8, "cm"), 
                     heatmap_height = unit(18, "cm"),
                     column_names_gp = grid::gpar(fontsize = 3),
                     row_names_gp = grid::gpar(fontsize = 5),
                     column_title = paste0("Top ", length(PCAdtmn)," proteins contribute\n to PC1 and PC2"))

gas.HM.hc

png(paste0("Figure/II.5.1. ",expdeg," - Heatmap on Top ",length(PCAdtmn)," PCA determiants.png"),width=7.25,height=8.25,units="in",res=300)
gas.HM.hc
draw(pL, x = unit(0.84, "npc"), y = unit(0.75, "npc"), just = c("right", "top"))
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
                     top_annotation = HeatmapAnnotation(Population = gas.DR$Population,
                                                        Datasets = colnames(gas.HM),
                                                        border = TRUE,gap = unit(0, "points"),
                                                        col = list(Population = c("Mt1-BFP+" = "dodgerblue",
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
#draw(pL, x = unit(0.96, "npc"), y = unit(0.785, "npc"), just = c("right", "top"))
dev.off()
