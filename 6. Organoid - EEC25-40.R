packages <- c(
  ##Data manipulation and analytics
  "magrittr","dplyr","knitr","tidyverse","scales",
  "scp","scater","scran","sva","SC3",
  #Visualization
  "RColorBrewer","dendsort","cowplot", "ggplot2","ComplexHeatmap","circlize",
  #Downstream analysis
  #"enrichplot","ggupset","clusterProfiler",
  #Mouse database
  #"org.Mm.eg.db")
)
  
lapply(packages, require, character.only = TRUE)

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("orga.EECData/orga.EEC.ProcessedData.Rda")

# I. DATA PROCESSING ----------

## Input ----------

###Input data

#Process MQ files from a MQ run
orga.EEC <- read.delim("orga.EECData/6. evidence - Organoid - EEC25-40.txt")
expdeg <- "Organoid - EEC"

#Cleaning
#Remove >2nd leading proteins
orga.EEC$Leading.razor.protein <- sub(';.*','',orga.EEC$Leading.razor.protein)
#Remove >2nd leading genes
orga.EEC$Gene.names <- sub(';.*','',orga.EEC$Gene.names)

#Combines 2 CVs into one (No need if set names for "Experiments" columns):
#orga.EEC$Raw.file <- sub('_FS45','',orga.EEC$Raw.file)
#orga.EEC$Raw.file <- sub('_FS65','',orga.EEC$Raw.file)

#Annotation file
sampleAnnotation <- read.csv("orga.EECData/6. sampleAnnotation Organoid - EEC.csv")

#Colour

#For heatmap
colour.hm <-  c(brewer.pal(9,'Blues')[9:3],brewer.pal(9,'Reds')[3:9])
#For Channels
colour.Chan <- c("dodgerblue","#E03E3E","#937264")
#For Datasets
colour.EEC <- c("#F52E00","#FF6003","#FE7701","#800BC6",
                "#FFAB26","#E5E3D5","#98D3DB","#5CAFC9","#437DBA","#0D98BB",
                "#2AA7BC","#47B6BC","#65C4BD","#FB33DB","#82D3BD","#9FE2BE")


## I.1. PSM level ----------

### 1. Reconstructing the dataset to create SCP/QFeatures object ----------

#Remove Empty columns in case mixed experiments with TMT 11 and TMT16-plex. 
#In this research, all the experiments were TMT16-plex

#read SCP
scp <- readSCP(featureData = orga.EEC,
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
scp <- zeroIsNA(scp, i = 1:l)

#### 2.2. Filter out failed runs based on PSM content (detected features) ----------

#Plot for checking number of PSM per assays
#Number of PSM per datasets 
nPSM <- dims(scp)[1, which(startsWith(names(scp),"E"))]

#Create data frame
nPSM <- data.frame(Datasets = names(nPSM),
                       Count = nPSM)
nPSM$Datasets <- factor(nPSM$Datasets, levels = nPSM$Datasets)

#Plot
nPSMplot <- ggplot(data=nPSM, aes(x=Datasets, y=Count, fill=Datasets)) +
  scale_fill_manual(values= colour.EEC, 
                    labels= nPSM$Datasets) +
  geom_bar(stat="identity")+
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
  ylim(0,10000)+
  xlab("Datasets") + 
  ylab("Number of detected PSM") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank())
nPSMplot

ggsave(paste0("Figure/2.2. ",expdeg," - Number of detected peptides per assay.png"), 
       plot=nPSMplot, width=8, height=5, dpi=320)


#Select the assays that have sufficient PSMs
#e.g. the number of rows is greater than 150),
keepAssay <- dims(scp)[1, ] > 150
scp <- scp[, , keepAssay]

#Suz:
hist(orga.EEC$PIF)
length(which(orga.EEC$Reverse == "+"))
length(which(orga.EEC$PIF<0.8))
length(which(orga.EEC$PIF == "NaN"))
length(which(orga.EEC$Potential.contaminant == "+"))

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

#### 2.5. Filter out PSMs with high sample to carrier ratio (for single-cell) ----------

#This experiment was performs on 100 cells, so there is no filtering features based on SCP metrics (single-cells metrics on Carrier and Reference channels)

#The number of channels per cell population (100 cells per channel) 

#### Number channels per pell population ----------
poputable <- t(as.data.frame(table(colData(scp)[, "SampleType"])))
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
#scp <- scp[, scp$SampleType %in% c("RFP", "BFP", "GFP", "Rest", "Empty"), ]

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

#Make SampleType as factor for plotting

scp$SampleType <- factor(scp$SampleType,
                          levels = c("100 EECs + stimulus","100 EECs - stimulus","Empty"))
#Plot
MedianRIplot <- scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.025) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.EEC) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  geom_hline(yintercept = 10^4.5, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianRI of Peptide Intensity"))
MedianRIplot

ggsave(paste0("Figure/6.1. ",expdeg," - MedianRIplot of Peptide Intensity.png"),
       plot=MedianRIplot, width=7.2,height=5.45,units="in",dpi=320)


#### 6.1.2. Median Relative Reporter Intensity ----------
scp6.1.2 <- scp

scp$RI <- factor(scp$RI,
                 levels = c("126","127N","127C","128N","128C","129N",'129C',
                            "130N","130C","131N","131C","132N","132C","133N","133C","134"))
colour.RI <- c("#2EC506","#2EC506","#F9D43E",
               "#2EC506","#2EC506","#2EC506","#2EC506","#2EC506","#2EC506",
               "#2EC506","#2EC506","#2EC506","#2EC506","#2EC506","#2EC506",
               "#EF6E18")

#Plot
MedianRIplot1 <- scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.025) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = RI),
              alpha = 0.7, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.RI) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  geom_hline(yintercept = 10^4.5, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianRI of Peptide Intensity"))
MedianRIplot1

ggsave(paste0("Figure/6.1.2 ",expdeg," - MedianRIplot of PSM - Channel.png"),
       plot=MedianRIplot1, width=7.2,height=5.45,units="in",dpi=320)

scp <- scp6.1.2
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
ts <- min(scp$MedianCV[scp$SampleType=="Empty"],na.rm = TRUE) - 0.001 # - 0.001 for drawinga line below the minium value

#
colour.EEC <- c("#85A655","#FF6003","#FE7701","#800BC6",
                "#FFAB26","#E5E3D5","#631B24","#5CAFC9","#437DBA","#0D98BB",
                "#2AA7BC","#47B6BC","#65C4BD","#FB33DB","#82D3BD","#9FE2BE")

#Plot
MedianCVplot <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianCV, fill = SampleType)) +
  geom_boxplot(width=0.4,outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.EEC)+
  ylim(0,1) +
  geom_hline(yintercept = 0.3, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianCV of Peptide Intensity"))
MedianCVplot

ggsave(paste0("Figure/6.2. ", expdeg," - MedianCVplot.png"), 
       plot=MedianCVplot, width=7.2,height=5.45,units="in",dpi=320)

#MedianCV in channels containning cell(s) should theoretically much more consistent 
#than for Empty channels.

scp6.2 <- scp

#### 6.2.2. Median Coefficient of Variation ----------
scp$RI <- factor(scp$RI,
                 levels = c("126","127N","127C","128N","128C","129N",'129C',
                            "130N","130C","131N","131C","132N","132C","133N","133C","134"))

colour.RI <- c("#2EC506","#2EC506","#F9D43E",
               "#070FE5","#2EC506","#2EC506","#2EC506","#2EC506","#DE231D",
               "#2EC506","#FB33DB","#2EC506","#2EC506","#2EC506","#2EC506",
               "#EF6E18")
#Plot
MedianCVplot2 <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianCV, fill = SampleType)) +
  geom_boxplot(width=0.4) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = RI),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.RI)+
  ylim(0,1) +
  geom_hline(yintercept = 0.3, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianCV of Peptide Intensity"))
MedianCVplot2

ggsave(paste0("Figure/6.2. ", expdeg," - MedianCVplot2.png"), 
       plot=MedianCVplot2, width=7.2,height=5.45,units="in",dpi=320)

#MedianCV in channels containning cell(s) should theoretically much more consistent 
#than for Empty channels.

scp6.2 <- scp



####  6.3. Filter based on the MedianRI and MedianCV ----------

#Based on the distribution of the Empty, ones can decides which MedianCV threshold and 
#remove channels/cells that exhibit high MedianCV over the different proteins. 
scp <- scp[, !is.na(scp$MedianCV) & scp$MedianCV < 0.3 & scp$MedianRI > 10^4.5, ]

#Optional: to remove some channels/cells
scp <- scp[, scp$SampleType != "Empty", ]

#Plot after filtering

colour.EEC1 <- c("#F52E00","#FF6003","#FE7701",
                "#FFAB26","#E5E3D5","#98D3DB","#5CAFC9","#437DBA","#0D98BB",
                "#2AA7BC","#47B6BC","#65C4BD","#82D3BD","#9FE2BE")

#MedianRI
MedianRIplotQC <- scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.025) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.EEC1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  #geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianRI of Peptide Intensity after QC"))
MedianRIplotQC

ggsave(paste0("Figure/6.3. ", expdeg," - MedianRIplot of PSM after QC.png"), 
              plot=MedianRIplotQC, width=7.2,height=4.45,units="in",dpi=320)


#MedianCV
MedianCVplotQC <- scp %>% getWithColData("Peptides") %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianCV, fill = SampleType)) +
  geom_boxplot(width=0.4) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan) +
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.001))+
  scale_colour_manual(values = colour.EEC1)+
  ylim(0,1) +
  geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle(paste0(expdeg," - MedianCV of Peptide Intensity after QC"))
MedianCVplotQC

ggsave(paste0("Figure/6.3. ", expdeg," - MedianCVplot after QC.png"), 
       plot=MedianCVplotQC, width=7.2,height=4.45,units="in",dpi=320)


#Return SampleType to character (required for Batch correction)
scp$SampleType <- as.character(scp$SampleType)

#Save
scp6 <- scp

#### Number channels per pell population ----------
poputable <- t(as.data.frame(table(colData(scp)[, "SampleType"])))
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
model <- model.matrix(~ as.character(SampleType), data = colData(sce))

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

addAssay(scp,
         y = getWithColData(scp,"Proteins_norm"),
         name = "Proteins_batchC") %>%
  addAssayLinkOneToOne(from = "Proteins_norm",
                       to = "Proteins_batchC") ->
  scp

#Picking only of Experiments/Assays of interest
orga.DR <- getWithColData(scp, "Proteins_batchC")

#Add a scale-log assays for some heatmap/visualization
assay(orga.DR,"scaleassay") <- t(scale(t(assay(orga.DR,"assay")), 
                                      center = TRUE, scale = TRUE))

#Find proteins that have no correspoding gene
blankgene <- orga.EEC[match(rownames(orga.DR), orga.EEC$Leading.razor.protein),][,c("Leading.razor.protein","Gene.names")]

#Search on  Uniprot
blankgene[which(blankgene[,"Gene.names"] == ""),]

#Change rownames orga.DR from proteins names to correspoding gene names
rownames(orga.DR) <- orga.EEC[match(rownames(orga.DR), orga.EEC$Leading.razor.protein),]$Gene.names

#Check order
#View(data.frame(rownames(scp[["Proteins_batchC"]]),
#                rownames(orga.DR),
#                orga.EEC[match(rownames(scp[["Proteins_batchC"]]), orga.EEC$Leading.razor.protein),]$Gene.names,d,
#                orga.EEC[match(rownames(scp[["Proteins_batchC"]]), orga.EEC$Leading.razor.protein),]$Leading.razor.protein))

save(orga.EEC,
     expdeg,
     sampleAnnotation,
     scp,
     orga.DR,
     MedianCVplot,
     MedianRIplot,
     colour.hm,
     colour.Chan,
     colour.EEC,
     file="orga.EECData/orga.EEC.ProcessedData.Rda")

#Save
#orga.DR1 <- orga.DR

# II. DIMENSIONALITY REDUCTION AND VISUALIZATION ----------

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("orga.EECData/orga.EEC.ProcessedData.Rda")


## 1. Principal component analysis ----------

interest <- c("CHGA", "SST", "GAST", "NTS", "MLN", "CCK", "GCG", "PYY")

#Choose n top proteins, scale/ non-Scale may generate different results
set.seed(100)
#n = 500, 1000, 1500, 2000
PCAtop <- 500
orga.DR <- runPCA(orga.DR,
                 ncomponents = 5,
                 ntop = PCAtop,
                 scale = FALSE,
                 exprs_values = "assay",
                 name = "PCA")

#Note: If scale=TRUE, the expression values for each feature are standardized 
#so that their variance is unity. This will also remove features with standard deviations below 1e-8.
#1.2. PCA plot
PCAplot <- plotReducedDim(orga.DR,
                          dimred = "PCA",
                          colour_by = "SampleType",
                          point_alpha = 0.8)+
  ggplot2::scale_color_manual(values = c("dodgerblue","#E03E3E"))+
  labs(colour = "SampleType")+
  ggtitle(paste0("PCA on top ",PCAtop, " most varible proteins"))
PCAplot

ggsave(paste0("Figure/II.1.1. ",expdeg," - PCA for top ",PCAtop," PCA determinants.png"), 
       plot=PCAplot, width=5, height=4, dpi=320)


#1.2. PCA Check Batch
orga.DR$Designs <- as.factor(orga.DR$Designs) 
PCAplotbatch <- plotReducedDim(orga.DR,
                               dimred = "PCA",
                               colour_by = "Designs",
                               point_alpha = 0.8)+
  ggtitle("Batch Correction")
PCAplotbatch

ggsave(paste0("Figure/II.1.2. ",expdeg," - PCA Batch for top ",PCAtop," PCA determinants.png"), 
       plot=PCAplotbatch, width=4, height=4, dpi=320)



## 2. Percentage of variance explained ----------

#Percentage of variance explained by each PCs
pv <- attr(reducedDim(orga.DR), "percentVar")
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
       plot=pvPlot, width=7,height=6, dpi=320)

#Top PCs (elbows) explain for %
pvplot$pv[1:3] %>% sum
pvplot$pv[1:elbow] %>% sum

#Take only 3 PC
reducedDim(orga.DR, "PCA.elbow") <- reducedDim(orga.DR)[,1:elbow]
#reducedDimNames(orga.DR)

#Check first 3 PC correlation
PCAtoplot <- plotPCA(orga.DR, 
                     ncomponents = elbow, 
                     colour_by = "SampleType")+
  ggplot2::scale_color_manual(values = c("dodgerblue","#E03E3E"))+
  labs(colour = "SampleType")
PCAtoplot

ggsave(paste0("Figure/II.2.2. ",expdeg," - PCA with more PCs plot.png"), 
       plot=PCAtoplot, width=7,height=6, dpi=320)

## 3. Uniform Manifold Approximation and Projection ----------

#850
set.seed(850)
orga.DR <- runUMAP(orga.DR,
                  ncomponents = elbow,
                  exprs_values = "assay",
                  n_neighbors = 4,
                  dimred = "PCA.elbow",
                  name = "UMAP")
UMAPplot <- plotReducedDim(orga.DR,
               dimred = "UMAP",
               colour_by = "SampleType",
               point_alpha = 1)+
  ggplot2::scale_color_manual(values = c("dodgerblue","#E03E3E"))+
  labs(colour = "SampleType")
UMAPplot

ggsave(paste0("Figure/II.3. ",expdeg," - UMAPplot.png"), 
              plot=UMAPplot, width=4,height=3,units="in",dpi=320)

## 4. Detect the PCA determinants ----------

## Detect PCA determinants/gene that contribute to the most variance PC - Based on Rotation 
pcrot <- attr(reducedDim(orga.DR),"rotation")

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
PCAdtmn <- which(match(rownames(orga.DR), listPC)>0)
length(PCAdtmn)


#Plot Expression of these PCA determinants
ntop <- 20

#non-Scale
PCAEplot <- plotExpression(orga.DR,
                           rownames(orga.DR)[head(PCAdtmn,20)],
                           exprs_values = "assay",
                           log2_values = FALSE,
                           colour_by = "SampleType",)+
  ggplot2::scale_color_manual(values = c("dodgerblue","#E03E3E"))+
  labs(colour = "SampleType") +
  ggtitle(paste0("Expression of top ",ntop," PCA determinent"))
PCAEplot

ggsave(paste0("Figure/II.4. ",expdeg," - Top ",ntop," PCA determinent - non-Scale.png"), 
       plot=PCAEplot, width=7,height=6, dpi=320)

#Scale
PCAEscaleplot <- plotExpression(orga.DR, 
                           rownames(orga.DR)[head(PCAdtmn,20)],
                           exprs_values = "scaleassay",
                           log2_values = FALSE,
                           colour_by = "SampleType",)+
  ggplot2::scale_color_manual(values = c("dodgerblue","#E03E3E"))+
  labs(colour = "SampleType")+
  ggtitle(paste0("Expression of top ",ntop," PCA determinent - Scale"))
PCAEscaleplot

ggsave(paste0("Figure/II.4. ",expdeg," - Top ",ntop," PCA determinent - Scale.png"), 
       plot=PCAEscaleplot, width=7,height=6, dpi=320)


## 5. Heatmap for PCA determinants/genes ----------

#Take data of n top PC determinants - Scale data (non-Scale is also possible)
orga.HM <- assay(orga.DR,"scaleassay")[PCAdtmn,]

#Shorten col names
colnames(orga.HM) <- sub('*Reporter.intensity.*','',colnames(orga.HM))

#Create Heatmap for top PC determinants

### Hierarchical clustering ----------

#(Wihout spliting: remove "column_split = 4,")
orga.HM.hc <- Heatmap(orga.HM, name = "Scale Log \nAbundance",
                     clustering_distance_rows = "pearson",
                     clustering_method_rows = "ward.D2",
                     cluster_columns = dendsort(hclust(dist(t(orga.HM)))),
                     col = colour.hm,
                     top_annotation = HeatmapAnnotation(Cells = orga.DR$SampleType,
                                                        Datasets = colnames(orga.HM),
                                                        border = TRUE,gap = unit(0, "points"),
                                                        col = list(Cells = c("100 EECs + stimulus" = "dodgerblue",
                                                                             "100 EECs - stimulus" = "#E03E3E"
                                                                             #"GFP" = "#53BB79",
                                                                             #"Rest" = "#B1C97F"),
                                                                   #Datasets = c("D0370_100_1"="#D88C4C",
                                                                    #            "D0370_100_2"="#E49D41",
                                                                    #            "D0370_100_3"="#FEB447",
                                                                    #            "D0386_100_1"="#06CDE0",
                                                                    #            "D0386_100_2"="#1498BE")),
                                                        #show_legend = FALSE
                                                        ))),
                     heatmap_width = unit(8, "cm"), 
                     heatmap_height = unit(18, "cm"),
                     column_names_gp = grid::gpar(fontsize = 3),
                     row_names_gp = grid::gpar(fontsize = 5),
                     column_title = paste0("Top ", length(PCAdtmn)," genes contribute\n to PC1 and PC2"))
                     
orga.HM.hc

png(paste0("Figure/II.5.1. ",expdeg," - Heatmap on Top ",length(PCAdtmn)," PCA determiants.png"),width=7.25,height=8.25,units="in",res=300)
orga.HM.hc
dev.off()

###For gene of interest
interest

orga.HM1 <- assay(orga.DR,"scaleassay")[rownames(orga.DR)[rownames(orga.DR) %in% interest],]
#Shorten col names
colnames(orga.HM1) <- sub('*Reporter.intensity.*','',colnames(orga.HM1))

orga.HM1.hc <- Heatmap(orga.HM1, name = "Scale Log \nAbundance",
                      clustering_distance_rows = "pearson",
                      clustering_method_rows = "ward.D2",
                      cluster_columns = dendsort(hclust(dist(t(orga.HM1)))),
                      col = colour.hm,
                      top_annotation = HeatmapAnnotation(Cells = orga.DR$SampleType,
                                                         Datasets = colnames(orga.HM),
                                                         border = TRUE,gap = unit(0, "points"),
                                                         col = list(Cells = c("100 EECs + stimulus" = "dodgerblue",
                                                                              "100 EECs - stimulus" = "#E03E3E"
                                                                              #"GFP" = "#53BB79",
                                                                              #"Rest" = "#B1C97F"),
                                                                              #Datasets = c("D0370_100_1"="#D88C4C",
                                                                              #            "D0370_100_2"="#E49D41",
                                                                              #            "D0370_100_3"="#FEB447",
                                                                              #            "D0386_100_1"="#06CDE0",
                                                                              #            "D0386_100_2"="#1498BE")),
                                                                              #show_legend = FALSE
                                                         ))),
                      heatmap_width = unit(18, "cm"), 
                      heatmap_height = unit(8, "cm"),
                      column_names_gp = grid::gpar(fontsize = 3),
                      row_names_gp = grid::gpar(fontsize = 5),
                      column_title = paste0("Top gene of interest"))

orga.HM1.hc

png(paste0("Figure/II.5.1. ",expdeg," - Heatmap on Top gene of interest.png"),width=10,height=4,units="in",res=300)
orga.HM1.hc
dev.off()

### Check with proteinGroup
proteG <- read.delim("orga.EECData/6. proteinGroups - Organoid - EEC25-40.txt")

colnames(proteG)

proteG$Majority.protein.IDs <- sub(';.*','',proteG$Majority.protein.IDs)
proteG$Gene.names <- sub(';.*','',proteG$Gene.names)


#Number of protein of interest found

proteG$Gene.names %in% interest %>% sum()
