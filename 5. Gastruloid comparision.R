packages <- c(
  #Filtering and analytics
  "scp","scales","MASS","viridisLite","matrixStats","sva",
  #Visualization
  "scales","cowplot","ggpubr","gridExtra","ggplot2","VennDiagram",
  #Data manipulation
  "plyr","dplyr")
lapply(packages, require, character.only = TRUE)

#Use Saved Processed data to pass the # I. DATA PROCESSING
load("GasData/gas.ProcessedData.Rda")

#File Annotation
#Data = Bulk, 100 (100 cells) or SC (Single Cells)
#gasData: MQ proteinGroups.txt file
#gas100.anno: annotation file for reporter labeling 
#gasData.scp: transfer for SCP manipulation
#gasData.f: data after quality control/filtering
#gasData.f.n: data after normalization
#gasData.f.n.l: data after log transformation
#gasData.corr: data for computing correlation among datasets

# I. Data processing -------------

## 1.1 Bulk with stringent filter - WRONG, NOT USED ANYMORE-----
### 1.1. Input rawfile -----

gasBulk <- read.delim("1. proteinGroups - Bulk - LFQ.txt")

#Filtering out Contaminant and Reverses proteins
gasBulk <- gasBulk[!grepl(("REV|CON"),gasBulk$Protein.IDs),]

#Choose only first proteins
gasBulk$Protein.IDs <- sub(';.*','',gasBulk$Protein.IDs)

### 1.2. Filtering out peptide and razor peptides-------
#Filter all proteins that have less than 2 peptides and
#no razor+unique peptide in the sum of your experiments
#Filtering very strict 
#Pick up columne of LFQ Intensity
meta <- c("BFP1","BFP2","BFP3",
          "GFP1","GFP2","GFP3",
          "Rest1","Rest2","Rest3",
          "RFP1","RFP2","RFP3")

data <- gasBulk
for (i in meta){
  design <- which(endsWith(colnames(data), i)) #Taking only columns belongs to a particular experimental design
  pep <- data[,paste0("Peptides.",i)]#The total number of peptide in the design
  razor <- data[,paste0("Razor...unique.peptides.",i)] #The total number of razor + unique peptides in the design
  data[pep <3 | razor <1, design] <- NA} #Change into NA

#Remove columns contains NA value. These columns were not passed the filter
gasBulk.rb <- data[complete.cases(data[ , which(startsWith(colnames(data), "Peptides."))]),]

#Count number of Proteins for Bulk
nPro.Bulk <- data.frame(Datasets = "BulkLFQ",
                        Count = nrow(gasBulk.rb))

### 1.3. Cleaning matrix --------------

# Choose LFQ intensity row
gasBulk.f <- gasBulk.rb[,which(startsWith(colnames(gasBulk.rb), "LFQ.intensity."))]

# Change names of rows and columns
rownames(gasBulk.f) <- gasBulk.rb$Protein.IDs
colnames(gasBulk.f) <- c("BFP","BFP","BFP",
                         "GFP","GFP","GFP",
                         "Rest","Rest","Rest",
                         "RFP","RFP","RFP")
# Zero into NA
gasBulk.f[gasBulk.f==0] <- NA

# Missing value plot
misval <- longFormat(gasBulk.f) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme_minimal()
misval

### 1.4. Protein Normalization ----

# Divide columns by median
gasBulk.f.n <- data.matrix(gasBulk.f)
gasBulk.f.n <- sweep(gasBulk.f.n, 
                    MARGIN = 2,
                    FUN = "/",
                    STATS = colMedians(gasBulk.f.n, na.rm = TRUE))
# Divide rows by mean
gasBulk.f.n <- sweep(gasBulk.f.n,
                    MARGIN = 1,
                    FUN = "/",
                    STATS = rowMeans(gasBulk.f.n, na.rm = TRUE))

### 1.5. Log transformation ----
gasBulk.f.n.l <- log2(gasBulk.f.n)
dim(gasBulk.f.n.l)

## Missing value plot
misval <- longFormat(gasBulk.f.n.l) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme_minimal()
misval

## 1.2 Bulk using Perseus -----
gasBulk.f.Perseus1 <- read.delim("GasData/1.2. gasBulk Persues.txt")
gasBulk.f.Perseus <- read.delim("GasData/1.2. gasBulk Persues.txt")
gasBulk.f.Perseus <- gasBulk.f.Perseus[-c(1:2),]
gasBulk.f.Perseus <- gasBulk.f.Perseus[,1:12]
gasBulk.f.Perseus <- as.data.frame(apply(gasBulk.f.Perseus, 2, as.numeric))

rownames(gasBulk.f.Perseus) <- sub(';.*','',gasBulk.f.Perseus1[-c(1:2),]$Protein.IDs)

colnames(gasBulk.f.Perseus) <- c("BFP","BFP","BFP",
                                 "GFP","GFP","GFP",
                                 "Rest","Rest","Rest",
                                 "RFP","RFP","RFP")






## 2. Gastruloid 100 cells -----


### 2.1. Input rawfile -----
gas100 <- read.delim("GasData/2&3. proteinGroups - 100 cells - 0386 + 0370.txt")

#Filtering out Contaminant and Reverses proteins
gas100 <- gas100[!grepl(("REV|CON"),gas100$Protein.IDs),]

#Choose only first proteins
gas100$Protein.IDs <- sub(';.*','',gas100$Protein.IDs)

#Change row names for corresponding proteins
rownames(gas100) <- gas100$Protein.IDs

#Annotation file for SCP data manipulation (later step)
gas100.anno <- read.csv("GasData/2&3. sampleAnnotation - Gas_100.csv")

#Colour
colour.100 <- c("#D88C4C","#E49D41","#FEB447","#06CDE0","#1498BE")

### 2.2. Filtering out peptides and razor peptides-------

#Create each dataframe (For testing values of rbind dataframe below)
#for (i in unique(gas100.anno$Designs)){
#  design <- which(endsWith(colnames(gas100), i)) #Taking only columns belongs to a particular experimental design
#  data <- gas100[, c(which(colnames(gas100)=="Protein.IDs"), design)] #create dataframe of the design
#  pep <- data[,paste0("Peptides.",i)] #The total number of peptide in the design
#  razor <- data[,paste0("Razor...unique.peptides.",i)] #The total number of razor + unique peptides in the design
#  data <- filter(data, 
#                 pep>2 & #Filter out all proteins that have less than 2 peptides 
#                   razor>0) #Filter out all proteins that have  no razor+unique peptides
#  assign(paste0(i), data) # Creating the dataframe for the design
#}


#sum(nrow(D0370_100_1),
#    nrow(D0370_100_2),
#    nrow(D0370_100_3),
#    nrow(D0386_100_1),
#    nrow(D0386_100_2))

gas100.rb = list()
for (i in unique(gas100.anno$Designs)){
  design <- which(endsWith(colnames(gas100), i)) #Taking only columns belongs to a particular experimental design
  data <- gas100[, c(which(colnames(gas100)=="Protein.IDs"), design)] #create dataframe of the design
  pep <- data[,paste0("Peptides.",i)] #The total number of peptide in the design
  razor <- data[,paste0("Razor...unique.peptides.",i)] #The total number of razor + unique peptides in the design
  data <- filter(data, 
                 pep>2 & #Filter out all proteins that have less than 2 peptides 
                   razor>0) #Filter out all proteins that have  no razor+unique peptides
  data$Designs <- i #Assign name of design for its protein (row)
  colnames(data) <- sub(paste0('.',i),'',colnames(data))
  gas100.rb[[i]] <- data
}


### 2.3. Transfer into SCP for manipulation-----
# Create rbind dataframe for SCP package
gas100.scp <- do.call(rbind, gas100.rb)

#### 3.3.1. Transfer into SCP -----
gas100.scp <- readSCP(featureData = gas100.scp,
                      colData = gas100.anno,
                      channelCol = "Channel",
                      batchCol = "Designs",
                      removeEmptyCols = TRUE)


#Order factor (For plot) and corespldong colour
gas100.scp$SampleType <- factor(gas100.scp$SampleType,
                               levels = c("BFP","RFP","GFP","Rest","Empty"))
colour.Chan <- c("dodgerblue","#E03E3E","#53BB79","#B1C97F","#937264","#3043A2","#C25D7B")

# Change zero into NA
gas100.scp <- zeroIsNA(gas100.scp, i = names(gas100.scp))

#Total number of channels
length(colData(gas100.scp)[, "SampleType"])

# Number of channels per cell populations
poputable <- t(as.data.frame(table(colData(gas100.scp)[, "SampleType"])))
colnames(poputable) <- poputable[1,]
poputable <- poputable[-1, ]
poputable <- t(as.data.frame(poputable))
rownames(poputable) <- c("Channels")
knitr::kable(poputable)

#### 3.3.2. Join in one assay---------
gas100.scp <- joinAssays(gas100.scp,
                         i = names(gas100.scp),
                         name = "Proteins")


#### 3.3.3. Quality control------

##### a. Filter based on the median relative intensity-------------

#Find MedianRI
gas100.scp[["Proteins"]] %>%
  assay %>%
  apply(MARGIN = 2, 
        FUN = median, 
        na.rm = TRUE) ->
  medians
colData(gas100.scp)[names(medians), "MedianRI"] <- medians

#Maxium of empty channels
ts <- max(gas100.scp$MedianRI[gas100.scp$SampleType=="Empty"])+ 0.01 #+0.01 for preventing zero value

#Plot
MedianRIplot <- gas100.scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.05) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.15, height = 0.2))+
  scale_colour_manual(values = colour.100) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  geom_hline(yintercept = ts, linetype='dashed', col = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("100 Cells - MedianRI of Protein Intensity")
MedianRIplot

ggsave("Figure/1.3. 100 Cells - MedianRIplot.png", plot=MedianRIplot, 
       width=7.2,height=4.45,units="in",dpi=300)

gas100.scp1 <- gas100.scp


##### b. Filter columne with MedianRI < empty------------
gas100.scp <- gas100.scp[, gas100.scp$MedianRI > ts, ]

gas100.scp2 <- gas100.scp

##### c. Re-check QC -------------
MedianRIplotQC <- gas100.scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.05) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = colour.100) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("100 Cells - QCMedianRI of Protein Intensity")
MedianRIplotQC

ggsave("Figure/1.3. 100 Cell - MedianRIplotQC .png", plot=MedianRIplotQC, 
       width=7.2,height=4.45,units="in",dpi = 300)

gas100.scp3 <- gas100.scp


##### d. Filter high missing rate -----
# Remove proteins per assays with high missing rate
gas100.scp <- filterNA(gas100.scp,
                       i = "Proteins",
                       pNA = 0.99)
#Save
gas100.scp4 <- gas100.scp

#### 3.3.4. Number of proteins per assays -----

#Pick single datasets 
nPro.100 <- dims(gas100.scp)[1, which(startsWith(names(gas100.scp),"D"))]

#Create data frame
nPro.100 <- data.frame(Datasets = names(nPro.100),
                       Count = nPro.100)

nPro.100$Datasets <- factor(nPro.100$Datasets, levels = nPro.100$Datasets)

#Plot
nPro.100plot <- ggplot(data=nPro.100, 
                       aes(x=Datasets, 
                           y=Count, 
                           fill=Datasets)) +
  scale_fill_manual(values= colour.100,
                    labels= nPro.100$Datasets) +
  geom_bar(stat="identity")+
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
  ylim(0,700)+
  xlab("Datasets") + 
  ylab("Number of detected proteins") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank())
nPro.100plot

ggsave("Figure/1.3. 100 cells - Number of Detected Proteins per Datasets2.png", 
       plot=nPro.100plot, width=4, height=4, dpi = 320)

#### 3.3.5. Convert SCP back into data frame-------

#Pick Proteins Assays
gas100.f <- gas100.scp[["Proteins"]]@assays@data@listData[[1]] %>% data.frame

#Create ProteinID column
gas100.f$Proteins <- gas100.scp[["Proteins"]]@rowRanges@elementMetadata@listData[["Protein.IDs"]]

#Aggregate columns based in ProteinID 
gas100.f <- ddply(gas100.f, .(Proteins), 
                  function(x) colSums(x[,grep("Reporter.intensity.", colnames(gas100.f))],na.rm = T)  )

#Change zero into NA
gas100.f[gas100.f==0] <- NA

##Create final data

#Row
rownames(gas100.f) <- gas100.f$Proteins

#Remove the Protein.ID
gas100.f <- gas100.f[,-which(names(gas100.f) =="Proteins")]

#Annotate colnames with SampleType
colnames(gas100.f) <- gas100.scp$SampleType


### 2.4. Protein Normalization ----

# Divide columns by median
gas100.f.n <- data.matrix(gas100.f)
gas100.f.n <- sweep(gas100.f.n, 
                    MARGIN = 2,
                    FUN = "/",
                    STATS = colMedians(gas100.f.n, na.rm = TRUE))

# Divide rows by mean
gas100.f.n <- sweep(gas100.f.n,
                    MARGIN = 1,
                    FUN = "/",
                    STATS = rowMeans(gas100.f.n, na.rm = TRUE))

### 2.5. Log transformation ----
gas100.f.n.l <- log2(gas100.f.n)

# Missing value plot
misval <- longFormat(gas100.f.n.l) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme_minimal()
misval


### 2.6.Imput or not? - NOT IMPUTE NOW-----
gas100.f.n %>%
  assay %>%
  is.na %>%
  mean

gas100.f.n.l.i <- impute::impute.knn(gas100.f.n.l,
                                     k = 4, rowmax = 1, colmax= 1,
                                     maxp = Inf, rng.seed = 1234)

gas100.f.n.l.i <- gas100.f.n.l.i$data
#Check (If 0 then fine)
gas100.f.n.l.i %>%
  assay %>%
  is.na %>%
  mean

misval <- longFormat(gas100.f.n.l.i) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme_minimal()
misval


### 2.7. Batch correction - NOT PERFORMED ----------
sce <- getWithColData(gas100.scp, "Proteins")
batch <- gas100.scp$lcbatch
model <- model.matrix(~ SampleType, data = colData(sce))

gas100.f.n.l.i.b <- ComBat(dat = gas100.f.n.l.i,
                           batch = batch,
                           mod = model)





## 3. Gastruloid Single Cells -----

## 3. Gastruloid Single Cells -----

### 3.1. Input rawfile -----
gasSC <- read.delim("GasData/4&5. proteinGroups - scLMNOIJK.txt")

#Filtering out Contaminant and Reverses proteins
gasSC <- gasSC[!grepl(("REV|CON"),gasSC$Protein.IDs),]

#Choose only first proteins
gasSC$Protein.IDs <- sub(';.*','',gasSC$Protein.IDs)

#Change row names for corresponding proteins
rownames(gasSC) <- gasSC$Protein.IDs

#Annotation file for SCP data manipulation (later step)
gasSC.anno <- read.csv("GasData/4&5. sampleAnnotation - Gas_sc.csv")

#Colour
colour.SC <- c("#776AD6","#A068B0","#E487B8",
               "#BBC385","#7FAD81","#0F7B6C","#0D5B11")

### 3.2. Filtering out peptides and razor peptides-------

##Create each dataframe (For testing values of rbind dataframe below)
#for (i in unique(gasSC.anno$Designs)){
#  design <- which(endsWith(colnames(gasSC), i)) #Taking only columns belongs to a particular experimental design
#  data <- gasSC[, c(which(colnames(gasSC)=="Protein.IDs"), design)] #create dataframe of the design
#  pep <- data[,paste0("Peptides.",i)] #The total number of peptide in the design
#  razor <- data[,paste0("Razor...unique.peptides.",i)] #The total number of razor + unique peptides in the design
#  data <- filter(data, 
#                 pep>2 & #Filter out all proteins that have less than 2 peptides 
#                   razor>0) #Filter out all proteins that have  no razor+unique peptides
#  assign(paste0(i), data) # Creating the dataframe for the design
#}

#sum(nrow(sc_I),
#    nrow(sc_J),
#    nrow(sc_K),
#    nrow(sc_L),
#    nrow(sc_M),
#    nrow(sc_N),
#    nrow(sc_O))

#Filtering out peptides and razor peptides
gasSC.rb = list()
for (i in unique(gasSC.anno$Designs)){
  design <- which(endsWith(colnames(gasSC), i)) #Taking only columns belongs to a particular experimental design
  data <- gasSC[, c(which(colnames(gasSC)=="Protein.IDs"), design)] #create dataframe of the design
  pep <- data[,paste0("Peptides.",i)] #The total number of peptide in the design
  razor <- data[,paste0("Razor...unique.peptides.",i)] #The total number of razor + unique peptides in the design
  data <- filter(data, 
                 pep>2 & #Filter out all proteins that have less than 2 peptides 
                   razor>0) #Filter out all proteins that have  no razor+unique peptides
  data$Designs <- i #Assign name of design for its protein (row)
  colnames(data) <- sub(paste0('.',i),'',colnames(data))
  gasSC.rb[[i]] <- data
}


### 3.3. Transfer into SCP for manipulation-----
# Create rbind dataframe for SCP package
gasSC.scp <- do.call(rbind, gasSC.rb)

#### 3.3.1. Transfer into SCP -----
gasSC.scp <- readSCP(featureData = gasSC.scp,
                     colData = gasSC.anno,
                     channelCol = "Channel",
                     batchCol = "Designs",
                     removeEmptyCols = TRUE)

#Order factor (For plot) and corespldong colour
gasSC.scp$SampleType <- factor(gasSC.scp$SampleType,
                               levels = c("BFP","RFP","GFP","Rest",
                                          "Empty","Reference","Carrier"))
colour.Chan <- c("dodgerblue","#E03E3E","#53BB79","#B1C97F","#937264","#3043A2","#C25D7B")

#Change zero into NA
gasSC.scp <- zeroIsNA(gasSC.scp, i = names(gasSC.scp))

#Total number of channels
length(colData(gasSC.scp)[, "SampleType"])

# Number of channels per cell populations
poputable <- t(as.data.frame(table(colData(gasSC.scp)[, "SampleType"])))
colnames(poputable) <- poputable[1,]
poputable <- poputable[-1, ]
poputable <- t(as.data.frame(poputable))
rownames(poputable) <- c("Channels")
knitr::kable(poputable)

gasSC.scp1 <- gasSC.scp
gasSC.scp <- gasSC.scp1

#### 3.3.2. Join in one assay---------
gasSC.scp <- joinAssays(gasSC.scp,
                        i = names(gasSC.scp),
                        name = "Proteins")
gasSC.scp2 <- gasSC.scp
gasSC.scp <- gasSC.scp2


#### 3.3.3. Quality control------

##### a. Filter based on the median relative intensity-------------

#Find MedianRI
gasSC.scp[["Proteins"]] %>%
  assay %>%
  apply(MARGIN = 2, 
        FUN = median, 
        na.rm = TRUE) ->
  medians
colData(gasSC.scp)[names(medians), "MedianRI"] <- medians

#Maxium of empty channels
ts <- max(gasSC.scp$MedianRI[gasSC.scp$SampleType=="Empty"])+ 0.01 #+0.01 for preventing zero value

#Plot 
MedianRIplot <- gasSC.scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.05) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = colour.SC) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("Single cells - MedianRI of Protein Intensity")
MedianRIplot

ggsave("Figure/1.3. Single Cell - MedianRIplot.png", plot=MedianRIplot, 
       width=7,height=4.25,units="in",dpi=300)

gasSC.scp1 <- gasSC.scp

##### b. Filter columns with MedianRI belong to single cell------------
#gasSC.scp <- gasSC.scp[, gasSC.scp$MedianRI > ts, ]
#Not using this filtering for Single Cells
#Choose only for channels of interests:
gasSC.scp <- gasSC.scp[, gasSC.scp$SampleType %in% c("BFP","RFP","GFP","Rest"), ]

##### c. Re-check QC -------------

MedianRIplotQC <- gasSC.scp %>% colData %>% data.frame %>%
  ggplot(aes(x = SampleType, y = MedianRI, fill = SampleType)) +
  #geom_boxplot(width=0.5) +
  geom_violin(adjust=0.75) + geom_boxplot(width=.05) +
  stat_summary(fun=mean, geom="point", fill="#D92721", shape=21, size=1)+
  scale_fill_manual(values = colour.Chan)+
  geom_jitter(aes(colour = Designs),
              alpha = 0.9, shape=16, size = 2,
              position = position_jitter(width = 0.2, height = 0.2))+
  scale_colour_manual(values = colour.SC) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)),
                limit = c(10^2, 10^6))+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("Single cells - MedianRI of Protein Intensity after QC")
MedianRIplotQC

ggsave("Figure/1.3. Single Cell - MedianRIplot QC.png", plot=MedianRIplotQC, 
       width=7.2,height=4.45,units="in",dpi = 300)

gasSC.scp2 <- gasSC.scp

##### d. Filter high missing rate -----
# Remove proteins per assays with high missing rate
gasSC.scp <- filterNA(gasSC.scp,
                      i = names(gasSC.scp),
                      pNA = 0.99)
#Save
gasSC.scp4 <- gasSC.scp

#### 3.3.4. Number of proteins per assays -----

#Pick single datasets 
nPro.SC <- dims(gasSC.scp)[1, which(startsWith(names(gasSC.scp),"sc_"))]

#Create data frame
nPro.SC <- data.frame(Datasets = names(nPro.SC),
                      Count = nPro.SC)

nPro.SC$Datasets <- factor(nPro.SC$Datasets, levels = nPro.SC$Datasets)

#Plot
nPro.SCplot <- ggplot(data=nPro.SC, aes(x=Datasets, y=Count, fill=Datasets)) +
  scale_fill_manual(values= colour.SC, 
                    labels= nPro.SC$Datasets) +
  geom_bar(stat="identity")+
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
  ylim(0,700)+
  xlab("Datasets") + 
  ylab("Number of detected proteins") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank())
nPro.SCplot

ggsave("Figure./1.3. Single Cell - Number of detected proteins.png", 
       plot=nPro.SCplot, width=4, height=4, dpi = 320)

gasSC.scp3 <- gasSC.scp

#### 3.3.5. Convert SCP back into data frame-------

#Pick Proteins Assays
gasSC.f <- gasSC.scp[["Proteins"]]@assays@data@listData[[1]] %>% data.frame

#Create ProteinID column
gasSC.f$Proteins <- gasSC.scp[["Proteins"]]@rowRanges@elementMetadata@listData[["Protein.IDs"]]

#Aggregate columns based in ProteinID 
gasSC.f <- ddply(gasSC.f, .(Proteins), 
                 function(x) colSums(x[,grep("Reporter.intensity.", colnames(gasSC.f))],na.rm = T)  )

#Change zero into NA
gasSC.f[gasSC.f==0] <- NA

##Create final data

#Row
rownames(gasSC.f) <- gasSC.f$Proteins

#Remove the Protein.ID
gasSC.f <- gasSC.f[,-which(names(gasSC.f) =="Proteins")]

#Annotate colnames with SampleType
colnames(gasSC.f) <- gasSC.scp$SampleType


### 3.4. Protein Normalization ----

# Divide columns by median
gasSC.f.n <- data.matrix(gasSC.f)
gasSC.f.n <- sweep(gasSC.f.n, 
                   MARGIN = 2,
                   FUN = "/",
                   STATS = colMedians(gasSC.f.n, na.rm = TRUE))

# Divide rows by mean
gasSC.f.n <- sweep(gasSC.f.n,
                   MARGIN = 1,
                   FUN = "/",
                   STATS = rowMeans(gasSC.f.n, na.rm = TRUE))


### 3.5. Log transformation ----
gasSC.f.n.l <- log2(gasSC.f.n)


# Missing value plot
misval <- longFormat(gasSC.f.n.l) %>%
  data.frame %>%
  group_by(colname) %>%
  summarize(missingness = mean(is.na(value))) %>%
  ggplot(aes(x = missingness)) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme_minimal()
misval

#Save Big data
save(gasBulk,
     gasBulk.rb,
     gasBulk.f,
     gasBulk.f.n,
     gasBulk.f.n.l,
     gasBulk.f.Perseus,
     gas100,
     gas100.anno,
     gas100.rb,
     gas100.scp,
     gas100.f,
     gas100.f.n,
     gas100.f.n.l,
     gasSC,
     gasSC.anno,
     gasSC.rb,
     gasSC.scp,
     gasSC.f,
     gasSC.f.n,
     gasSC.f.n.l,
     nPro.Bulk,
     nPro.100,
     nPro.SC,
     colour.Chan,
     colour.100,
     colour.SC,
     file="GasData/gas.ProcessedData.Rda")
#
load("GasData/gas.ProcessedData.Rda")

# II. Overlap ------
#
load("GasData/gas.ProcessedData.Rda")

## Venn Diagram -----
colour.Com <- c("#FD8289","lightslategrey","#0962A6")

#Use from Perseus

gasBulk.f <-  gasBulk.f.Perseus

#Number of each experiments

nrow(gasBulk.f)
nrow(gas100.f)
nrow(gasSC.f)

#Venn diagram
venn.diagram(list(rownames(gas100.f),
                  rownames(gasBulk.f),
                  rownames(gasSC.f)),
             category.names = c("100 Cells","BulkLFQ","Single Cells"),
             filename = 'Figure/2.1. Overlapped Protein among Datasets with Bulk Perseus.png',
             fill = colour.Com,
             # print.mode=c("raw","percent"),
             # Output features
             imagetype="png" ,
             height = 1000 , 
             width = 1000 , 
             resolution = 300,
             compression = "lzw",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 30, 0),
             cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "sans",
             rotation = 3
)




## Plot Number of Proteins ----

nPro.Bulk[,2] <- nrow(gasBulk.f.Perseus)

colour.All <- c("lightslategrey",colour.100,colour.SC)
nPRO <- rbind(nPro.Bulk,nPro.100,nPro.SC)
nPRO$Datasets <- as.factor(nPRO$Datasets)

#Mean & SD for 100 cells

mean(nPRO[2:6,2])
sd(nPRO[2:6,2])

#Mean & SD for SC cells

mean(nPRO[7:13,2])
sd(nPRO[7:13,2])

#Plot
nPro.plot <- ggplot(data=nPRO, 
                      aes(x=Datasets, 
                          y=Count, 
                          fill=Datasets)) +
  scale_fill_manual(values= colour.All,
                    labels= nPRO$Datasets) +
  geom_bar(stat="identity")+
  ylim(0,5000)+
  geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Number of detected proteins") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank())
nPro.plot

ggsave("Figure/2.2. Number of detected proteins.png", 
       plot=nPro.plot, width=8, height=4, dpi = 320)


# III.Protein abundance estimates ----
is.character(gasBulk.f)
gasBulk.f[is.numeric(as.matrix(gasBulk.f))]

gasBulk.f[is.nan(as.matrix(gasBulk.f))]<-NA

gasBulk.f[ gasBulk.f == "NaN" ] <- NA

bulkPE <- gasBulk.f
bulkPEm <-rowMedians(as.matrix(gasBulk.f), na.rm=T)

#Find replicate protein in Gastruloid 100 cells
bulkPE$Gas100 <- 0
for(X in rownames(gas100.f)){
  bulkPE$Gas100[grep(X, rownames(bulkPE))]<-1
}

#Find replicate protein in Gastruloid single cells
bulkPE$GasSC <- 0
for(X in rownames(gasSC.f)){
  bulkPE$GasSC[grep(X, rownames(bulkPE))]<-1
}

#
sum(bulkPE$Gas100==1)
sum(bulkPE$GasSC==1)

#Check
#bulkPEm2<-bulkPEm/median(bulkPEm, na.rm=T)*50000
#bulkPEm2<-bulkPEm

#median(bulkPEm2, na.rm=T)
#median(bulkPEm2, na.rm=T)

bulkPEm100<-bulkPEm[bulkPE$Gas100==1]
bulkPEmSC<-bulkPEm[bulkPE$GasSC==1]

bulkPEdf<-data.frame(Quant = c(bulkPEm, bulkPEm100,bulkPEmSC),
                     Datasets = c(rep("BulkLFQ", length(bulkPEm)), 
                                  rep("100 Cells", length(bulkPEm100)),
                                  rep("Single Cells", length(bulkPEmSC))))

#Remove NA
bulkPEdf<-bulkPEdf[!is.na(bulkPEdf$Quant), ]
bulkPEdf$Datasets<-as.factor(bulkPEdf$Datasets)
#bulkPEdf$Quant<-log10(bulkPEdf$Quant)
#bulkPEdf$q2<-10^bulkPEdf$Quant
colour.Com
#Plot
pro_es <- ggplot(bulkPEdf,aes(x=Quant,group=Datasets,fill=Datasets))+
  geom_histogram(position="dodge",binwidth=0.125,alpha = 1)+
  scale_fill_manual(values=colour.Com) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size=18) ) + 
  xlab("Protein Abundance") + 
  ylab("Number of proteins") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)) ) + 
  theme_minimal()
pro_es

ggsave("Figure/3. Subsets Protein Abundance.png", plot=pro_es, width=6, height=4, dpi = 320)

#Percetage of subsets 100 cells over Bulk
print(paste0(round(380/3274*100,2),"%"))

#Percetage of subsets single cells over Bulk
print(paste0(round(74/3274*100,2),"%"))

#Percetage of remaining subsets 100 cells
print(paste0(round((827-380)/827*100,2),"%"))

#Percetage of remaining subsets single cells
print(paste0(round((205-74)/205*100,2),"f%"))

#Totol proteins that are not detected in Bulk
print(323+124+7)

# IV. Absolute Intensity ----

## Absolute Intensity Per Channel ----
df.absI <- list()

#Bulk
K <- gasBulk.f
dataI <- K %>% unlist %>% data.frame() # only focus on intensity 
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "BulkLFQ"
df.absI[["BulkLFQ"]] <- dataI

#100 Cells
K <- gas100.f
dataI <- K %>% unlist %>% data.frame() # only focus on intensity 
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "100 Cells"
df.absI[["100 Cells"]] <- dataI

#Single Cells
K <- gasSC.f
dataI <- K %>% unlist %>% data.frame() # only focus on intensity 
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "Single Cells"
df.absI[["Single Cells"]] <- dataI

#Combine
gas.absI <- do.call(rbind, df.absI)
#
colnames(gas.absI) <- c("Intensity", "Datasets")
gas.absI$Datasets <- as.factor(gas.absI$Datasets)
gas.absI$log10 <- log10(gas.absI$Intensity)


#FD8289 - orange

#Plot

intsplot <- ggplot(gas.absI, aes(Intensity,fill = Datasets)) + 
  geom_histogram(position="dodge",binwidth=0.125,alpha = 1)+
  scale_fill_manual(values=colour.Com)+
  xlab("Protein Abundance") + 
  ylab("Number of proteins") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)) ) + 
  theme_minimal()
intsplot

ggsave("Figure/4.1. Absolute Intensity.png", plot=intsplot, width=6, height=4, dpi = 320)


## Absolute Intensity Per Proteins ----

df.absIM <- list()

#Bulk
K <- gasBulk.f
dataI <-rowMedians(as.matrix(K), na.rm=T) %>% data.frame() 
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "BulkLFQ"
df.absIM[["BulkLFQ"]] <- dataI

#100 Cells
K <- gas100.f
dataI <-rowMedians(as.matrix(K), na.rm=T) %>% data.frame()
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "100 Cells"
df.absIM[["100 Cells"]] <- dataI

#Single Cells
K <- gasSC.f
dataI <-rowMedians(as.matrix(K), na.rm=T) %>% data.frame() 
dataI[dataI==0]<-NA # Change Zero into NA
dataI <- dataI[!is.na(dataI)] %>% data.frame() # remove NA value
dataI$Datasets <- "Single Cells"
df.absIM[["Single Cells"]] <- dataI

#Combine
gas.absIM <- do.call(rbind, df.absIM)
#
colnames(gas.absIM) <- c("Intensity", "Datasets")
gas.absIM$Datasets <- as.factor(gas.absIM$Datasets)
gas.absIM$log10 <- log10(as.numeric(gas.absIM$Intensity))

#Plot

intsplotP <- ggplot(gas.absIM, aes(Intensity,fill = Datasets)) + 
  geom_histogram(position="dodge",binwidth=0.125,alpha = 1)+
  scale_fill_manual(values=colour.Com)+
  xlab("Protein Abundance") + 
  ylab("Number of proteins") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) {10^x}),
                labels = trans_format("log10", math_format(10^.x)) ) + 
  theme_minimal()
intsplotP

ggsave("Figure/4.2. Absolute Intensity.png", plot=intsplotP, width=8, height=4, dpi = 320)


# V. Bulk / 100 Cells/ Single Cells-----------

#Caution: Still need improve for automation; therefore loading processed data provided
load("GasData/gas.corr.Rda")
## 1. Data manipulation-----
#df.cor for creating the correlation table of the loop
df.cor <- list()
##cp: Comparison: 
#5.1 for Bulk vs 100 Cells, 
#5.2 for Bulk vs Single Cells,
#5.3 for 100 Cells vs Single Cells,
##i1|i2: #"BFP"|"RFP"|"GFP"|"Rest"
##data1|data2: Bulk | 100 Cells | Single Cells


#Manual replace
cp <- "5.3 " 
i1 <- "GFP"
i2 <- "Rest"
data1 <- "100 Cells"
data1.RA <- gas100.f.n
data2 <- "Single Cells"
data2.RA <- gasSC.f.n


#Loop should start here
##Data1
#SCoPE2: mat.input.bulk<-cr_norm(filt.mat.rc(dataI, 0.6, 0.8))
#Compute ratios
data1.RA.i1 <- rowMeans(data1.RA[,which(colnames(data1.RA)==i1)], na.rm=T)
data1.RA.i2 <- rowMeans(data1.RA[,which(colnames(data1.RA)==i2)], na.rm=T)
data1.RA.i <- data1.RA.i1/data1.RA.i2

##Data2
#SCoPE2: mat.input.bulk<-cr_norm(filt.mat.rc(dataI, 0.6, 0.8))
#Compute ratios
data2.RA.i1 <- rowMeans(data2.RA[,which(colnames(data2.RA)==i1)], na.rm=T )
data2.RA.i2 <- rowMeans(data2.RA[,which(colnames(data2.RA)==i2)], na.rm=T )
data2.RA.i <- data2.RA.i1/data2.RA.i2

##Compute ratio 

# Compile into a data frame
data.RA <- data.frame(Data1<-data1.RA.i)
data.RA <- cbind(data.RA, Data2<-data2.RA.i[ match(names(data1.RA.i),names(data2.RA.i)) ] )
colnames(data.RA) <- c(data1,data2)

#Log and NA replace
data.RA.l <- log2(data.RA[, 1:2])
data.RA.l[data.RA.l == -Inf] <- NA
data.RA.l[data.RA.l == Inf] <- NA

#Using code of SCoPE2 analysis
k <- 60
x <- data.RA.l[,data2]
y <- data.RA.l[,data1]
temp.df<-data.frame(x, y)
# Record the positions of the NA values in either x or y
na.v<-c()
for(i in 1:nrow(temp.df)){
  na.v<-c( na.v, any(is.na(temp.df[i,])) )
  
}
# Remove those points
temp.df.na<-temp.df[!na.v, ]

x<-temp.df.na$x
y<-temp.df.na$y

contour_cols <- viridis(k, alpha = 0.5)
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dens <- get_density(x, y, k)
## 2. Generating plots -----
def.par<-par()
#png(paste0(cp ,data1," vs ",data2," — ",i1," - ",i2, "Protein ratios.png"), width=4,height=4,units="in",res=300)
par(mar=c(5,6,2,2))
plot(x, y, col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], pch = 16, xlab="", ylab="",
     xlim=c(-2.3,3.5), ylim=c(-1.3,2.5),
     cex.axis=2)
text(-2.5, 2.30, paste0(i1,"/",i2, " (",sum(!is.na(data.RA.l[,2])),")"),pos = 4, cex=1.75)
text(-2.5, 1.80, paste0("ρ = ", round(cor(x,y),2)), pos = 4, cex=2)
mtext(paste0(data1,", log2"), side=2, padj = -3, cex=1.75)
mtext(paste0(data2,", log2"), side=1, padj = 3, cex=1.75)
#abline(a=0, b=1)
abline(lm(y~x), col="black", lty=2)
#print(lm(y~x))
#dev.off()

df.cor[[paste0(cp ,data1," vs ", data2," — ",i1," - ",i2)]] <- 
  data.frame(Data = paste0(data1," vs ",data2),
                               Data1 = data1,
                               Data2= data2,
                               Pop1 = i1,
                               Pop2 = i2,
                               Overlapped = sum(!is.na(data.RA.l[,2])),
                               CoverageD1 = sum(!is.na(data.RA.l[,2]))*100/nrow(data1.RA),
                               CoverageD2 = sum(!is.na(data.RA.l[,2]))*100/nrow(data2.RA),
                               Pearson = round(cor(x,y, method = "pearson"),2),
                               Spearman = round(cor(x,y, method = "spearman"),2),
                               Kendall = round(cor(x,y, method = "kendall"),2))


## Loop ends here

gas.corr <- do.call(rbind, df.cor)
gas.corr$Data <- as.factor(gas.corr$Data)



save(df.cor,
     gas.corr,
     file="GasData/gas.corr.Rda")
## 3. Load processed data -----------
load("GasData/gas.corr.Rda")

# VI. Plot Correlation -----------

#Colour 
colour.corr <- c("#E05700","#1266C0","#349CA5")

## Pearman -----
Pearsonplot <- ggplot(gas.corr, 
             aes(x = Data, y = Pearson, colour = Data, fill = Data)) +
  geom_violin(trim=T, alpha = 0.1) +
  geom_jitter(aes(colour = Data),
              shape=16,
              size = 4,
              position = position_jitter(width = 0.15, height = 0.05),
              alpha = 0.8) + 
  scale_fill_manual(values=colour.corr) +
  scale_colour_manual(values=colour.corr) +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y =  element_blank()) +
  theme(legend.position = "none") +
  ggtitle("Pearson")
Pearsonplot

ggsave("Figure/6.1. Pearsonplot.png", 
       plot=Pearsonplot, width=3, height=3, dpi = 320)

## Spearman-----
Spearmanplot <- ggplot(gas.corr, 
             aes(x = Data, y = Spearman, colour = Data, fill = Data)) +
  geom_violin(trim=T, alpha = 0.1) +
  geom_jitter(aes(colour = Data),
              shape=16,
              size = 4,
              position = position_jitter(width = 0.15, height = 0.05),
              alpha = 0.8) + 
  scale_fill_manual(values=colour.corr)+
  scale_colour_manual(values=colour.corr)+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y =  element_blank()) +
  theme(legend.position = "none") +
  ggtitle("Spearman")
Spearmanplot

ggsave("Figure/6.1. Spearmanplot.png", 
       plot=Spearmanplot, width=3, height=3, dpi = 320)

## Kendall-----
Kendallplot <- ggplot(gas.corr, 
             aes(x = Data, y = Kendall, colour = Data, fill = Data)) +
  geom_violin(trim=T, alpha = 0.1) +
  geom_jitter(aes(colour = Data),
              shape=16,
              size = 4,
              position = position_jitter(width = 0.15, height = 0.05),
              alpha = 0.8) + 
  scale_fill_manual(values=colour.corr)+
  scale_colour_manual(values=colour.corr)+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y =  element_blank()) +
  ggtitle("Kendall") 
Kendallplot

## Overlapped proteins -----
Overplot <- ggplot(gas.corr, aes(x = Data,
                                 y = Overlapped)) +
  geom_jitter(aes(colour = Data),
              size = 4,
              shape=16,
              position = position_jitter(width = 0.2, height = 1),
              alpha = 0.7) +
  scale_colour_manual(values=colour.corr)+ 
  ylim(0,400)+
  ylab("Overlapped") + 
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y =  element_blank()) +
  theme(legend.position = "none") +
  ggtitle("Number of Overlapped Proteins")   
Overplot

ggsave("Figure/6.1. Overplot.png", plot=Overplot, width=3, height=3, dpi = 320)

## Legend
cor.leg <- as_ggplot(leg <-get_legend(Kendallplot))

## Combined plot -----
gas.corr.plot <- plot_grid(Pearsonplot, 
                           Spearmanplot, 
                           Overplot,
                           cor.leg, ncol = 2, nrow = 2)
gas.corr.plot

ggsave("Figure/6.1. Datasets correlation comparison.png", 
       plot=gas.corr.plot, width=10, height=8, dpi = 320)

## Extracted Table -----

gas.table <- gas.corr

colnames(gas.table)[7] <- c("%Coverage for Data 1")
colnames(gas.table)[8] <- c("%Coverage for Data 2")

png("Figure/6.2. GasdataTable2.png",width=14,height=6,units="in",res=300)
grid.table(gas.table[,-1])
dev.off()


## Percent plot---------------
#Take out corresponding columns 
d1 <- gas.corr[,c(1,2,7)]
d2 <- gas.corr[,c(1,3,8)]
#Change names
colnames(d1) <- c("Comparison","Data","Coverage")
colnames(d2) <- c("Comparison","Data","Coverage")
#Combine into long data
gas.corr.c <- rbind(d1,d2)
gas.corr.c$Data <- as.factor(gas.corr.c$Data)
#Plot
coverplot <- ggplot(gas.corr.c, aes(x = Comparison,
                     y = Coverage)) +
  geom_jitter(aes(colour = Data),
              size = 4,
              shape=16,
              position = position_jitter(width = 0.2, height = 1),
              alpha = 0.8) +
  scale_colour_manual(values=colour.Com) +
  ylim(0,100) +
  ylab("% Coverage") +
  theme_minimal() +
  theme(#axis.text.x= element_blank()),
        #axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
  )
coverplot 

ggsave("Figure/6.3. Percentage coverage for data comparison.png", plot = coverplot, width=5.75, height=4.25, dpi = 320)

