#################Analysis of DDA Controlled Mixture for all 6 methods###########################
##Tool: PD

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

# Reading data
##############################
raw <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/DDA_CM_PD/ControlMixture_DDA_ProteomeDiscoverer_input.csv", stringsAsFactors=F) # the data file

annot <- read.csv('C:/Users/tushi/OneDrive/Desktop/Data/DDA_CM_PD/ControlMixture_DDA_ProteomeDiscoverer_annotation.csv')


###########################MSstats#########################################################################

##############################
## Make MSstats required format
##############################
quant <- PDtoMSstatsFormat(raw, 
                           annotation=annot,
                           removeProtein_with1Peptide=TRUE)

##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################

processed.quant <- dataProcess(quant,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               censoredInt="NA",
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


#Results
TP<-test.MSstats %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #199
FP<-test.MSstats %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #269


##################################### DeqMS###################################################################

input_data=quant[quant$ProteinName!="",]
## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(log2Intensity = log2(sum(Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

df.LFQ[sapply(df.LFQ, is.infinite)] <- NA

#########################################################################################################
#DEqms

# Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:4])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[5:7])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[8:10])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[11:13])))
df.LFQ$na_count_5 = apply(df.LFQ,1,function(x) sum(is.na(x[14:16])))



# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                         df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2 |df.LFQ$na_count_5<2,]
# df.LFQ.filter = df.LFQ

row.names(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)


# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count <- as.data.frame(quant )%>% select(ProteinName, PeptideModifiedSequence)
pep.count<-pep.count%>%group_by(ProteinName)%>%
  summarise(count=n_distinct(PeptideModifiedSequence),.groups = 'drop')%>%
  as.data.frame()
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count$count = pep.count$count+1
rownames(pep.count) = pep.count$ProteinName
pep.count = pep.count %>% select(-ProteinName)

#protein.matrix[sapply(protein.matrix, is.infinite)] <- NA

## DEqMS analysis on LFQ data
protein.matrix <- na.omit(log2(as.matrix(df.LFQ.filter[1:15])))

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


TP<-iprg_deqms_results %>% filter(grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #186
FP<-iprg_deqms_results %>% filter(!grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #58


###################################MSqRob2########################################################

input_data=quant[quant$ProteinName!="",]
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ$PSM = paste(input_data$PeptideModifiedSequence, input_data$PrecursorCharge, sep="_")
df.LFQ = as.data.table(df.LFQ)

## Convert into MSqRob2 format
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

TP<-iprg_msqrob_results %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #188
FP<-iprg_msqrob_results %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #127


########################################## PMartR #################################################################



e_mdata<-quant[,c('PeptideModifiedSequence','Run','Intensity')]
e_mdata<-e_mdata%>%group_by(PeptideModifiedSequence,Run)%>% summarise(Intensity=sum(Intensity))
e_mdata<-e_mdata %>% pivot_wider(names_from = Run, values_from = Intensity,values_fill = NA)

f<-data.frame(i=c('sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12','sample_13','sample_14','sample_15'))

vec<-data.frame(c=colnames(e_mdata))
annot$Run=gsub('[.]', '', annot$Run)
f_mdata<-annot[match(vec[2:16,],annot$Run),]
f_mdata<-cbind(f,f_mdata)
test<-quant[,c('PeptideModifiedSequence','ProteinName')]
test %>%group_by(PeptideModifiedSequence) %>%
  fill(ProteinName, .direction = "downup") %>%
  ungroup %>%
  distinct
reorder_idx <- match(e_mdata$PeptideModifiedSequence,test$PeptideModifiedSequence)
test <- test[match(e_mdata$PeptideModifiedSequence,test$PeptideModifiedSequence),]  

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


TP<-m %>% filter(grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow()  #99
FP<-m %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #165


################################## DEP ##########################################################################
input_data=quant[quant$ProteinName!="",]

df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")

data_unique <- make_unique(df.LFQ, "ProteinName", "ProteinName", delim = ";")

expdesign<-select(annot,c('Run',"Condition","BioReplicate"))

#expdesign$Run=str_split_fixed(expdesign$Run, "_01_", 2)[,2]
colnames(expdesign)=c('label','condition','replicate')  

#expdesign$label=paste0(expdesign$condition,"_",expdesign$replicate)
columns <- grep("Intensity", colnames(data_unique))
data<-make_se(data_unique,columns,expdesign)



# Filter, normalize and impute missing values
filt <- filter_missval(data, thr = 0)


norm <- normalize_vsn(filt)
imputed <- DEP::impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- DEP::test_diff(imputed, type="all")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Generate a results table
data_results <- get_results(dep)
#data_results

###############results
df<-select(data_results,c("ID",(matches("p.adj"))))
df<-df %>% pivot_longer(cols=c((matches("p.adj"))),
                        names_to='label',
                        values_to='p.adj')



TP<-df %>% filter(grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #201
FP<-df %>% filter(!grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #124



######################################################################################################

input_data=quant

df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
protein_table<-as.data.frame(df.LFQ)

intensity_colnames <- grep("^Intensity\\.", colnames(protein_table), value=TRUE)
# Create the intensity matrix
abundance_matrix <- as.matrix(protein_table[, intensity_colnames])
# Adapt column and row maxquant_protein_table
colnames(abundance_matrix) <- sub("^Intensity\\.", "", intensity_colnames)
rownames(abundance_matrix) <- protein_table$ProteinName


## Replace the zeros with NAs and take the log2() of the data

abundance_matrix[abundance_matrix == 0] <- NA
abundance_matrix <- log2(abundance_matrix)

## Normalize the data using median_normalization()
normalized_abundance_matrix <- median_normalization(abundance_matrix)

## Fit the probabilistic dropout model
sample_info_df <- data.frame(name = colnames(normalized_abundance_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- annot$Condition
sample_info_df$replicate <- annot$BioReplicate


fit <- proDA(normalized_abundance_matrix, design = sample_info_df$condition, 
             col_data = sample_info_df)
result_names(fit)

# Test which proteins differ between condition 
test_res1 <- proDA::test_diff(fit,contrast ="Condition1-Condition2")
test_res2 <- proDA::test_diff(fit,contrast ="Condition1-Condition3")
test_res3 <- proDA::test_diff(fit,contrast ="Condition1-Condition4")
test_res4 <- proDA::test_diff(fit,contrast ="Condition1-Condition5")
test_res5 <- proDA::test_diff(fit,contrast ="Condition2-Condition3")
test_res6 <- proDA::test_diff(fit,contrast ="Condition2-Condition4")
test_res7 <- proDA::test_diff(fit,contrast ="Condition2-Condition5")
test_res8 <- proDA::test_diff(fit,contrast ="Condition3-Condition4")
test_res9 <- proDA::test_diff(fit,contrast ="Condition3-Condition5")
test_res10 <- proDA::test_diff(fit,contrast ="Condition4-Condition5")

test<-rbind(test_res1,test_res2,test_res3,test_res4,test_res5,test_res6,test_res7,
            test_res8,test_res9,test_res10)

TP<-test %>% filter(grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #165
FP<-test %>% filter(!grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #86

######################################################################################################