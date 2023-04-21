#################Analysis of DDA Controlled Mixture for all 8 methods###########################

##Tool:Skyline
## This file is with native pre-processing of data

##Significant proteins in data
proteins<-c("P02701","P00711","Q29443","Q29550","P0CG53",
            "P68082","P00432","P02754","P24627","P80025",
            "P00915","P02787","P02663","P01008","P00921",
            "P05307","P61769","P02662","P01012","P02666",
            "Q3SX14","P00563","P02769","Q58D62","P00698",
            "P00004","P00711","P00442","P01133","P02753")

#Libraries Required

library(DEP)
library(proDA)
library(pmartR)
library(msmsTests)
library(reshape2)
library(MSstats)
library(dplyr)
library(tidyverse)
library(stringr)
library(MSnbase)
library(msqrob2)
library(MSstats)
library(matrixStats)
library(limma)
library(DEqMS)
library(data.table)
library(prolfqua)

# Reading data
raw <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/DDA_CM_Skyline/ControlMixture_DDA_Skyline_input.csv")
annotation <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/DDA_CM_Skyline/ControlMixture_DDA_Skyline_annotation.csv")

###########################MSstats#########################################################################

##############################
## Extra filtering for this dataset :
## 1: remove 'precursor -64'
unique(raw$Fragment.Ion)
raw <- raw[which(raw$Fragment.Ion %in% c( "precursor", "precursor [M+1]","precursor [M+2]")), ]


##############################
## Make MSstats required format
##############################
quant <- SkylinetoMSstatsFormat(raw,
                                annotation = annotation,
                                fewMeasurements="remove", ## same as default
                                removeProtein_with1Feature = TRUE)

##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################

processed.quant <- dataProcess(quant,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               censoredInt="0",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

##############################
## Model-based comparison + adjust p-value
##############################

comparison1<-matrix(c(1,-1,0,0,0),nrow=1)
comparison2<-matrix(c(1,0,-1,0,0),nrow=1)
comparison3<-matrix(c(1,0,0,-1,0),nrow=1)
comparison4<-matrix(c(1,0,0,0,-1),nrow=1)
comparison5<-matrix(c(0,1,-1,0,0),nrow=1)
comparison6<-matrix(c(0,1,0,-1,0),nrow=1)
comparison7<-matrix(c(0,1,0,0,-1),nrow=1)
comparison8<-matrix(c(0,0,1,-1,0),nrow=1)
comparison9<-matrix(c(0,0,1,0,-1),nrow=1)
comparison10<-matrix(c(0,0,0,1,-1),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3, comparison4, comparison5, 
                  comparison6, comparison7, comparison8, comparison9, comparison10)
row.names(comparison)<-c("M1-M2", "M1-M3", "M1-M4", "M1-M5", "M2-M3", 
                         "M2-M4", "M2-M5", "M3-M4", "M3-M5", "M4-M5")


colnames(comparison)<-c('Condition1','Condition2','Condition3','Condition4','Condition5')
test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

# Calculating TP and FP
TP<-test.MSstats %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #173
FP<-test.MSstats %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& adj.pvalue < .05) %>% nrow() #92

##################################### DeqMS#######################################################################

input_data = raw %>% select(Protein.Name, Peptide.Modified.Sequence, Precursor.Charge, File.Name, Area)


input_data$PSM = paste(input_data$Peptide.Modified.Sequence, input_data$Precursor.Charge, sep="_")
input_data = input_data %>% group_by(Protein.Name, File.Name, PSM) %>% summarise(Intensity=max(Area))
input_data = input_data %>% filter(Protein.Name != "")
input_data = input_data %>% dplyr::rename(ProteinName=Protein.Name, Run = File.Name)


## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ$Intensity = as.numeric(df.LFQ$Intensity)
#run_median = df.LFQ %>% group_by(Run) %>% summarize(run_med = median(Intensity, na.rm=TRUE))
#overall_med = median(run_median$run_med)
#df.LFQ$Intensity = df.LFQ$Intensity + (overall_med - df.LFQ$Intensity)

df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)

## Convert into DEqMS format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

# df.LFQ = medianSummary(df.LFQ, group_col = 1, ref_col=2)

# df.LFQ[2:13] = equalMedianNormalization(df.LFQ[2:13])
df.LFQ[sapply(df.LFQ, is.infinite)] <- NA

#########################################################################################################
#DEqms

# Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:4])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[5:7])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[8:10])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[11:13])))
df.LFQ$na_count_5 = apply(df.LFQ,1,function(x) sum(is.na(x[14:16])))


#rownames(df.LFQ)<-df.LFQ$ProteinName
#df.LFQ %>% filter(row.names(df.LFQ) %in% iprg_spike_prot)

# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                         df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2 |df.LFQ$na_count_5<2,]
# df.LFQ.filter = df.LFQ

row.names(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)

# df.LFQ.filter = equalMedianNormalization(df.LFQ.filter)

# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count <- as.data.frame(raw )%>% select(Protein.Name, Peptide.Modified.Sequence)
pep.count<-pep.count%>%group_by(Protein.Name)%>%
  summarise(count=n_distinct(Peptide.Modified.Sequence),.groups = 'drop')%>%
  as.data.frame()
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count$count = pep.count$count+1
rownames(pep.count) = pep.count$Protein.Name
pep.count = pep.count %>% select(-Protein.Name)

#protein.matrix[sapply(protein.matrix, is.infinite)] <- NA

## DEqMS analysis on LFQ data
protein.matrix <- as.matrix(df.LFQ.filter[1:15])

class = as.factor(rep(1:5, each=3))
design = model.matrix(~0+class) # fitting without intercept

fit1 <-lmFit(protein.matrix,design = design)
cont1 <- makeContrasts(class1-class2, levels = design)
cont2 <- makeContrasts(class1-class3, levels = design)
cont3 <- makeContrasts(class1-class4, levels = design)
cont4 <- makeContrasts(class1-class5, levels = design)
cont5 <- makeContrasts(class2-class3, levels = design)
cont6 <- makeContrasts(class2-class4, levels = design)
cont7 <- makeContrasts(class2-class5, levels = design)
cont8 <- makeContrasts(class3-class4, levels = design)
cont9 <- makeContrasts(class3-class5, levels = design)
cont10 <- makeContrasts(class4-class5, levels = design)
cont = cbind(cont1, cont2, cont3, cont4, cont5, cont6,cont7,cont8,cont9,cont10)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count[rownames(fit3$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error first
min(fit3$count)

fit4 = DEqMS::spectraCounteBayes(fit3)

## Analyze results
results_list = list()
for (i in seq_along(colnames(fit4$coefficients))){
  iprg_deqms_results = DEqMS::outputResult(fit4,coef_col = i)
  iprg_deqms_results$Label = colnames(fit4$coefficients)[[i]]
  results_list[[i]] = iprg_deqms_results
}
iprg_deqms_results = rbindlist(results_list)


TP<-iprg_deqms_results %>% filter(grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #136
FP<-iprg_deqms_results %>% filter(!grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #100

##########################################MSqRob2####################################################################

input_data = raw %>% select(Protein.Name, Peptide.Modified.Sequence, Precursor.Charge, File.Name, Area)
input_data$PSM = paste(input_data$Peptide.Modified.Sequence, input_data$Precursor.Charge, sep="_")
input_data = input_data %>% group_by(Protein.Name, File.Name, PSM) %>% summarise(Intensity=max(Area))
input_data = input_data %>% filter(Protein.Name != "")
input_data = input_data %>% dplyr::rename(ProteinName=Protein.Name, Run = File.Name)
input_data$Intensity = as.numeric(input_data$Intensity)

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, PSM, Run, Intensity)
df.LFQ = as.data.table(df.LFQ)

## Convert into MSqRob format
df.LFQ = dcast(df.LFQ, ProteinName + PSM~Run, 
               value.var = "Intensity", fun.aggregate = max, fill=NA)



colnames(df.LFQ)[3:17] = paste("Intensity", colnames(df.LFQ)[3:17], sep="_")

# Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:4])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[5:7])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[8:10])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[11:13])))
df.LFQ$na_count_5 = apply(df.LFQ,1,function(x) sum(is.na(x[14:16])))

# Filter protein table. require minimum two values for each group.
df.LFQ = df.LFQ[df.LFQ$na_count_1<3 | df.LFQ$na_count_2<3 |
                  df.LFQ$na_count_3<3 | df.LFQ$na_count_4<3| df.LFQ$na_count_5<3,]
ecols <- grep("Intensity", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 1, ecol = ecols,  name = "peptideRaw")

# cond <- which((colnames(pe)[[1]][1], split = "_")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- rep(c("1", "2", "3", "4","5"), each=3) %>% as.factor()
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

## Data preprocessing
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")

Protein_filter <- rowData(pe[["peptideLog"]])$ProteinName %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$ProteinName)
pe <- pe[Protein_filter,,]

pe <- filterFeatures(pe, ~ nNonZero >= 2)

pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)

pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "ProteinName",
                        name = "protein"
)

pe <- msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L1 <- makeContrast("condition2=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4","condition5"))
L2 <- makeContrast("condition3=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4","condition5"))
L3 <- makeContrast("condition4=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4","condition5"))
L4 <- makeContrast("condition5=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4","condition5"))
L5 <- makeContrast("condition3-condition2=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4","condition5"))
L6 <- makeContrast("condition4-condition2=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4","condition5"))
L7 <- makeContrast("condition5-condition2=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4","condition5"))
L8 <- makeContrast("condition4-condition3=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4","condition5"))
L9 <- makeContrast("condition5-condition3=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4","condition5"))
L10 <- makeContrast("condition5-condition3=0", parameterNames = c("(Intercept)",
                                                                  "condition2",
                                                                  "condition3",
                                                                  "condition4","condition5"))
L=cbind(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10)


pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

results_list = list()
for (i in 4:12){
  df = rowData(pe[["protein"]])[[i]]
  iprg_msqrob_results = tibble::rownames_to_column(df, "Protein")
  iprg_msqrob_results$Label = colnames(rowData(pe[["protein"]]))[[i]]
  results_list[[i]] = iprg_msqrob_results
}
iprg_msqrob_results = rbindlist(results_list)

TP<-iprg_msqrob_results %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #137
FP<-iprg_msqrob_results %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #51

########################################### PMartR #############################################################

annot<-annotation

e_mdata<-raw[,c('Peptide.Modified.Sequence','File.Name','Area')]
e_mdata$Area=as.numeric(e_mdata$Area)
e_mdata<-e_mdata%>%group_by(Peptide.Modified.Sequence,File.Name)%>% summarise(Intensity=sum(Area))

e_mdata = e_mdata %>% dplyr::rename(PeptideSequence=Peptide.Modified.Sequence, Run = File.Name)
e_mdata<-e_mdata %>% pivot_wider(names_from = Run, values_from = Intensity,values_fill = NA)

f<-data.frame(i=c('sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12','sample_13','sample_14','sample_15'))

vec<-data.frame(c=colnames(e_mdata))
f_mdata<-annot[match(vec[2:16,],annot$Run),]
f_mdata<-cbind(f,f_mdata)
test<-raw[,c('Peptide.Modified.Sequence','Protein.Name')]
test = test %>% dplyr::rename(PeptideSequence=Peptide.Modified.Sequence, ProteinName = Protein.Name)
test %>%group_by(PeptideSequence) %>%
  fill(ProteinName, .direction = "downup") %>%
  ungroup %>%
  distinct
reorder_idx <- match(e_mdata$PeptideSequence,test$PeptideSequence)
test <- test[match(e_mdata$PeptideSequence,test$PeptideSequence),]  

colnames(test)<-c("PeptideSequence","Protein")

colnames(e_mdata)<-c("PeptideSequence",'sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12','sample_13','sample_14','sample_15')

mypepData <- as.pepData(e_data=e_mdata, f_data = f_mdata, e_meta=test,edata_cname = "PeptideSequence", emeta_cname="Protein",fdata_cname = "i", data_scale = "abundance",check.names=FALSE)
########Created pep data object########################################################################

### Format Data ###
# replace any 0's with NA's #
mypepData <- edata_replace(omicsData = mypepData, x = 0, y = NA) # there are no 0's
# log2 transform the data #
mypepData <- edata_transform(omicsData = mypepData, data_scale = "log2")
# identify the grouping variable, Status #
mypepData <- group_designation(omicsData = mypepData, main_effects = "Condition",
                               covariates = NULL)

### Filter Peptides ###
# Molecule Filter #
myfilter <- molecule_filter(omicsData = mypepData)
mypepData <- applyFilt(filter_object = myfilter, omicsData = mypepData, min_num = 2)


# Proteomics Filter #
myfilter <- proteomics_filter(omicsData = mypepData)
mypepData <- applyFilt(filter_object = myfilter, omicsData = mypepData,min_num_peps=2)



# IMD-ANOVA Filter #
myfilter <- imdanova_filter(omicsData = mypepData)

mypepData <- applyFilt(filter_object = myfilter, omicsData = mypepData,
                       min_nonmiss_anova = 2, min_nonmiss_gtest = 3)


##Normalize
mypepData <-  normalize_global(omicsData = mypepData,
                               subset_fn = "all",
                               norm_fn = "median",
                               apply_norm = TRUE,
                               backtransform = TRUE)


### Protein Quantification ###
myprodata <- protein_quant(pepData = mypepData, method = "rollup", combine_fn =
                             "median")
### Statistical Analysis ###
comparisons <- create_comparisonDF(comp_type="pairwise",myprodata)

mystats <- imd_anova(omicsData = myprodata, comparisons = comparisons,
                     test_method = "anova",pval_adjust_a="bonferroni")



#Results
m<-mystats %>% select(Protein,starts_with('P_Value'))

m<-m %>% pivot_longer(cols=starts_with('P_Value'),names_to='label',values_to='p-value')


TP<-m %>% filter(grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #131
FP<-m %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& `p-value` < .05) %>% nrow() #142


########################################### DEP ####################################################################



df.LFQ = raw %>% select(Protein.Name, File.Name, Area)

df.LFQ = df.LFQ %>% group_by(Protein.Name, File.Name) %>% summarise(Intensity = max(Area, na.rm = TRUE))
df.LFQ$Intensity=as.numeric(df.LFQ$Intensity)
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "Protein.Name", timevar = "File.Name", direction = "wide")

data_unique <- make_unique(df.LFQ, "Protein.Name", "Protein.Name", delim = ";")

expdesign<-annotation
colnames(expdesign)=c('label','condition',"replicate")   
expdesign$replicate=str_split_fixed(expdesign$label,".raw",2)[,1]
expdesign$replicate=str_split_fixed(expdesign$replicate,"_",11)[,11]

columns <- grep("Intensity.", colnames(data_unique))
data<-make_se(data_unique,columns,expdesign)



# Filter, normalize and impute missing values
filt <- filter_missval(data, thr = 1)


norm <- normalize_vsn(data)
imputed <- DEP::impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- DEP::test_diff(data, type="all")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Generate a results table
data_results <- get_results(dep)

###############results
df<-select(data_results,c("ID",(matches("p.adj"))))
df<-df %>% pivot_longer(cols=c((matches("p.adj"))),
                        names_to='label',
                        values_to='p.adj')

TP<-df %>% filter(grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #31
FP<-df %>% filter(!grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #151



######################################proDA###############################################################


########################################################################################################

df.LFQ = raw %>% select(Protein.Name, File.Name, Area)

df.LFQ = df.LFQ %>% group_by(Protein.Name, File.Name) %>% summarise(Intensity = max(Area, na.rm = TRUE))
df.LFQ$Intensity=as.numeric(df.LFQ$Intensity)
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "Protein.Name", timevar = "File.Name", direction = "wide")

protein_table<-df.LFQ
colnames(protein_table)[1] ="Protein.IDs"
intensity_colnames <- grep("^Intensity\\.", colnames(protein_table), value=TRUE)
# Create the intensity matrix
abundance_matrix <- as.matrix(protein_table[, intensity_colnames])
# Adapt column and row maxquant_protein_table
colnames(abundance_matrix) <- sub("^Intensity\\.", "", intensity_colnames)
rownames(abundance_matrix) <- protein_table$Protein.IDs


## Replace the zeros with NAs and take the log2() of the data

abundance_matrix[abundance_matrix == 0] <- NA
abundance_matrix <- log2(abundance_matrix)

## Normalize the data using median_normalization()
normalized_abundance_matrix <- median_normalization(abundance_matrix)

## Fit the probabilistic dropout model
sample_info_df <- data.frame(name = colnames(normalized_abundance_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- str_split_fixed(sample_info_df$name,'_',11)[,10]
sample_info_df$replicate <- str_split_fixed(sample_info_df$name,'_',11)[,11]


fit <- proDA(normalized_abundance_matrix, design = sample_info_df$condition, 
             col_data = sample_info_df)
result_names(fit)

# Test which proteins differ between condition 
test_res1 <- proDA::test_diff(fit,contrast ="`1`-`2`")
test_res2 <- proDA::test_diff(fit,contrast ="`1`-`3`")
test_res3 <- proDA::test_diff(fit,contrast ="`1`-`4`")
test_res4 <- proDA::test_diff(fit,contrast ="`1`-`5`")
test_res5 <- proDA::test_diff(fit,contrast ="`2`-`3`")
test_res6 <- proDA::test_diff(fit,contrast ="`2`-`4`")
test_res7 <- proDA::test_diff(fit,contrast ="`2`-`5`")
test_res8 <- proDA::test_diff(fit,contrast ="`3`-`4`")
test_res9 <- proDA::test_diff(fit,contrast ="`3`-`5`")
test_res10 <- proDA::test_diff(fit,contrast ="`4`-`5`")

test<-rbind(test_res1,test_res2,test_res3,test_res4,test_res5,test_res6,test_res7,
            test_res8,test_res9,test_res10)

TP<-test %>% filter(grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #8
FP<-test %>% filter(!grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #11



######################################prolfqua###############################################################


########################################################################################################


#Preprocessing
annot<-annotation
colnames(annot)<-c("File.Name","Condition","BioReplicate")
raw$Area<-as.numeric(raw$Area)
startdata<-raw%>%group_by(Protein.Name,File.Name)%>%summarise(Area=sum(Area,na.rm = TRUE))

startdata <- dplyr::inner_join(annot, startdata, by = "File.Name")


atable <- AnalysisTableAnnotation$new()

atable$fileName = "File.Name"
atable$hierarchy[["protein_Id"]] <- c("Protein.Name")
atable$factors[["dilution."]] = "Condition"
atable$set_response("Area")
config <- AnalysisConfiguration$new(atable)
adata <- setup_analysis(startdata, config)


lfqpep<- prolfqua::LFQData$new(adata, config)

#Filtering
lfqpep$remove_small_intensities()$filter_proteins_by_peptide_count()

lfqpep<-lfqpep$get_Transformer()$log2()$robscale()$lfq

#ModelFit
lfqpep$config$table$get_response()
formula_Condition <-  strategy_lm("transformedIntensity ~ dilution.")

# specify model definition
modelName  <- "Model"
unique(lfqpep$data$dilution.)


contrasts<-c("1_vs_2"="dilution.Condition2-dilution.Condition1",
             "1_vs_3"="dilution.Condition3-dilution.Condition1",
             "1_vs_4"="dilution.Condition4-dilution.Condition1",
             "1_vs_5"="dilution.Condition5-dilution.Condition1",
             "2_vs_3"="dilution.Condition3-dilution.Condition2",
             "2_vs_4"="dilution.Condition4-dilution.Condition2",
             "2_vs_5"="dilution.Condition5-dilution.Condition2",
             "3_vs_4"="dilution.Condition4-dilution.Condition3",
             "3_vs_5"="dilution.Condition5-dilution.Condition3",
             "4_vs_5"="dilution.Condition5-dilution.Condition4")

mod <- prolfqua::build_model(
  lfqpep$data,
  formula_Condition,
  subject_Id = lfqpep$config$table$hierarchy_keys() )

aovtable <- mod$get_anova()
head(aovtable)



contr<-prolfqua::Contrasts$new(mod,contrasts)

conI<-prolfqua::ContrastsMissing$new(lfqpep,contrasts)


contrasts<-prolfqua::merge_contrasts_results(prefer=contr,add=conI)


prolfqua::ContrastsModerated$undebug("get_contrasts")
contrasts<-contrasts$merged |> prolfqua::ContrastsModerated$new()

results<-contrasts$to_wide()


#Fetching the results

m<-results %>% select(protein_Id,starts_with('FDR'))

m<-m %>% pivot_longer(cols=starts_with('FDR'),names_to='label',values_to='p-value')

TP<- m%>%filter(grepl(paste(proteins, collapse="|"),protein_Id) & `p-value` < .05) %>% nrow()#190
FP<-m%>%filter(!grepl(paste(proteins, collapse="|"),protein_Id) & `p-value` < .05) %>% nrow()#159


######################################Limma###############################################################


########################################################################################################




input_data = raw %>% select(Protein.Name, Peptide.Modified.Sequence, Precursor.Charge, File.Name, Area)


input_data$PSM = paste(input_data$Peptide.Modified.Sequence, input_data$Precursor.Charge, sep="_")
input_data = input_data %>% group_by(Protein.Name, File.Name, PSM) %>% summarise(Intensity=max(Area))
input_data = input_data %>% filter(Protein.Name != "")
input_data = input_data %>% dplyr::rename(ProteinName=Protein.Name, Run = File.Name)


## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ$Intensity = as.numeric(df.LFQ$Intensity)
#run_median = df.LFQ %>% group_by(Run) %>% summarize(run_med = median(Intensity, na.rm=TRUE))
#overall_med = median(run_median$run_med)
#df.LFQ$Intensity = df.LFQ$Intensity + (overall_med - df.LFQ$Intensity)

df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)

## Convert into DEqMS format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

# df.LFQ = medianSummary(df.LFQ, group_col = 1, ref_col=2)

# df.LFQ[2:13] = equalMedianNormalization(df.LFQ[2:13])
df.LFQ[sapply(df.LFQ, is.infinite)] <- NA


row.names(df.LFQ) = df.LFQ$ProteinName
df.LFQ = df.LFQ %>% select(-ProteinName)


protein.matrix <- as.matrix(df.LFQ[1:15])

class = as.factor(rep(1:5, each=3))
design = model.matrix(~0+class) # fitting without intercept

fit1 <-lmFit(protein.matrix,design = design)
cont1 <- makeContrasts(class1-class2, levels = design)
cont2 <- makeContrasts(class1-class3, levels = design)
cont3 <- makeContrasts(class1-class4, levels = design)
cont4 <- makeContrasts(class1-class5, levels = design)
cont5 <- makeContrasts(class2-class3, levels = design)
cont6 <- makeContrasts(class2-class4, levels = design)
cont7 <- makeContrasts(class2-class5, levels = design)
cont8 <- makeContrasts(class3-class4, levels = design)
cont9 <- makeContrasts(class3-class5, levels = design)
cont10 <- makeContrasts(class4-class5, levels = design)
cont = cbind(cont1, cont2, cont3, cont4, cont5, cont6,cont7,cont8,cont9,cont10)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)


m<-as.data.frame(fit3$p.value)
m$protein<-rownames(fit3$p.value)
m<-m %>% pivot_longer(cols=starts_with('class'),names_to='label',values_to='p-value')
m$adj.pvalue = p.adjust(m$`p-value`,method = 'BH')


TP<- m%>%filter(grepl(paste(proteins, collapse="|"),protein) & adj.pvalue < .05) %>% nrow()#140
FP<-m%>%filter(!grepl(paste(proteins, collapse="|"),protein) & adj.pvalue < .05) %>% nrow()#93
