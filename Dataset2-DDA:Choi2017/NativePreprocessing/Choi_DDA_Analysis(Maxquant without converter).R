#################Analysis of DDA choi for all 6 methods###########################
##Tool: Maxquant

##Significant proteins in data
proteins<-c("P44015","P44752","P44374","P44983","P44683","P55249")

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
##############################
# read in proteinGroups file, in order to use protein ids
proteinGroups<-read.table("C:/Users/tushi/OneDrive/Desktop/Data/DDA_choi_maxquant/Choi2017_DDA_MaxQuant_proteinGroups.txt", sep="\t", header=TRUE)

# Read in MaxQuant file: evidence.txt
infile <- read.table("C:/Users/tushi/OneDrive/Desktop/Data/DDA_choi_maxquant/Choi2017_DDA_MaxQuant_evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates
annot <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/DDA_choi_maxquant/Choi2017_DDA_MaxQuant_annotation.csv", header=TRUE)


##################################MSstats########################################################

##############################
## Make MSstats required format
##############################
quant <- MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups,
                             useUniquePeptide = TRUE,
                             summaryforMultipleRows = max,
                             removeProtein_with1Peptide=TRUE)

##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################

## censoredInt='NA' for MaxQuant
processed.quant <- dataProcess(quant,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               censoredInt="NA",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

##############################
## Model-based comparison + adjust p-value
##############################

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")
colnames(comparison)<-c('Condition1','Condition2','Condition3','Condition4')

test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

# Calculating TP and FP
TP<-test.MSstats %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #21
FP<-test.MSstats %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& adj.pvalue < .05) %>% nrow() #27


##################################### DeqMS###################################################################




input_data = infile %>% select(Proteins, Sequence, Charge, Reverse,
                               Potential.contaminant, Raw.file, Intensity)
input_data = input_data[!input_data$Reverse=="+",]
input_data = input_data[!input_data$Potential.contaminant=="+",]

input_data$PSM = paste(input_data$Sequence, input_data$Charge, sep="_")
input_data = input_data %>% group_by(Proteins, Raw.file, PSM) %>% summarise(Intensity=max(Intensity))
input_data = input_data %>% filter(Proteins != "")
input_data = input_data %>% dplyr::rename(ProteinName=Proteins, Run = Raw.file)

## Just use MSstats converter to make life easier
# input_data = MaxQtoMSstatsFormat(evidence, annotation, df.prot, removeFewMeasurements = FALSE)

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)


df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = (sum(Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

## Convert into DEqMS format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

# df.LFQ = medianSummary(df.LFQ, group_col = 1, ref_col=2)

# df.LFQ[2:13] = equalMedianNormalization(df.LFQ[2:13])
df.LFQ[sapply(df.LFQ, is.infinite)] <- NA


# Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:4])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[5:7])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[8:10])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[11:13])))


df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                         df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2,]


row.names(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)


# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count <- proteinGroups %>% select(Majority.protein.IDs, Unique.peptides)

# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count$Unique.peptides = pep.count$Unique.peptides+1
rownames(pep.count) = pep.count$Majority.protein.IDs
pep.count = pep.count %>% select(-Majority.protein.IDs)

## DEqMS analysis on LFQ data
protein.matrix = na.omit(as.matrix(df.LFQ.filter[1:12]))

class = as.factor(rep(1:4, each=3))
design = model.matrix(~0+class) # fitting without intercept

fit1 <-lmFit(protein.matrix,design = design)
cont1 <- makeContrasts(class1-class2, levels = design)
cont2 <- makeContrasts(class1-class3, levels = design)
cont3 <- makeContrasts(class1-class4, levels = design)
cont4 <- makeContrasts(class2-class3, levels = design)
cont5 <- makeContrasts(class2-class4, levels = design)
cont6 <- makeContrasts(class3-class4, levels = design)
cont = cbind(cont1, cont2, cont3, cont4, cont5, cont6)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count[rownames(fit3$coefficients),"Unique.peptides"]

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


TP<-iprg_deqms_results %>% filter(grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #8
FP<-iprg_deqms_results %>% filter(!grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #3











##################################MSqRob2############################################################


input_data = infile %>% filter(Potential.contaminant != "+") %>% select(Proteins, Sequence, Charge, Raw.file, Intensity)
input_data$PSM = paste(input_data$Sequence, input_data$Charge, sep="_")
input_data = input_data %>% group_by(Proteins, Raw.file, PSM) %>% summarise(Intensity=max(Intensity))
input_data = input_data %>% filter(Proteins != "")
input_data = input_data %>% dplyr::rename(ProteinName=Proteins, Run = Raw.file)

# input_data$PSM = paste(input_data$PeptideSequence, input_data$PrecursorCharge, sep="_")

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, PSM, Run, Intensity)
df.LFQ = as.data.table(df.LFQ)

## Convert into DEqMS format
df.LFQ = dcast(df.LFQ, ProteinName + PSM~Run, 
               value.var = "Intensity", fun.aggregate = max, fill=NA)
colnames(df.LFQ)[3:14] = paste("Intensity", colnames(df.LFQ)[3:14], sep="_")

## Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[3:5])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[6:8])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[9:11])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[12:14])))

# Filter protein table. require minimum two values for each group.
df.LFQ = df.LFQ[df.LFQ$na_count_1<3 | df.LFQ$na_count_2<3 |
                  df.LFQ$na_count_3<3 | df.LFQ$na_count_4<3,]
ecols <- grep("Intensity", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 2, ecol = ecols,  name = "peptideRaw")

# cond <- which((colnames(pe)[[1]][1], split = "_")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- rep(c("1", "2", "3", "4"), each=3) %>% as.factor()
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
                                                      "condition4"))
L2 <- makeContrast("condition3=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4"))
L3 <- makeContrast("condition4=0", parameterNames = c("(Intercept)",
                                                      "condition2",
                                                      "condition3",
                                                      "condition4"))
L4 <- makeContrast("condition3-condition2=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4"))
L5 <- makeContrast("condition4-condition2=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4"))
L6 <- makeContrast("condition4-condition3=0", parameterNames = c("(Intercept)",
                                                                 "condition2",
                                                                 "condition3",
                                                                 "condition4"))
L=cbind(L1,L2,L3,L4,L5,L6)

pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

results_list = list()
for (i in 4:9){
  df = rowData(pe[["protein"]])[[i]]
  iprg_msqrob_results = tibble::rownames_to_column(df, "Protein")
  iprg_msqrob_results$Label = colnames(rowData(pe[["protein"]]))[[i]]
  results_list[[i]] = iprg_msqrob_results
}
iprg_msqrob_results = rbindlist(results_list)


TP<-iprg_msqrob_results %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #19
FP<-iprg_msqrob_results %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #9


################################# PMartR ###################################################################



e_mdata<-infile %>% filter(Potential.contaminant != "+") %>% select(Sequence, Raw.file, Intensity)
e_mdata<-e_mdata%>%group_by(Sequence,Raw.file)%>% summarise(Intensity=sum(Intensity))
e_mdata<-e_mdata %>% pivot_wider(names_from = Raw.file, values_from = Intensity,values_fill = NA)

vec<-data.frame(c=colnames(e_mdata))
f_mdata<-annot[match(vec[2:13,],annot$Raw.file),]
test<-infile[,c('Sequence','Proteins')]
test %>%group_by(Sequence) %>%
  fill(Proteins, .direction = "downup") %>%
  ungroup %>%
  distinct
reorder_idx <- match(e_mdata$Sequence,test$Sequence)
test <- test[match(e_mdata$Sequence,test$Sequence),]  

colnames(test)<-c("Sequence","Protein")

mypepData <- as.pepData(e_data=e_mdata, f_data = f_mdata, e_meta=test,edata_cname = "Sequence", emeta_cname="Protein",fdata_cname = "Raw.file", data_scale = "abundance",check.names=FALSE)
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


TP<-m %>% filter(grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow()  #20
FP<-m %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #227

############################################## DEP #################################################################




input_data<-infile %>% filter(Potential.contaminant != "+") %>% select(Proteins,Gene.names, Raw.file, Intensity)

df.LFQ = as.data.frame(input_data) %>% select(Proteins, Gene.names,Raw.file, Intensity)
df.LFQ = df.LFQ %>% group_by(Proteins,Gene.names,Raw.file) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
df.LFQ = df.LFQ %>% filter(Proteins != "")
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = c("Proteins","Gene.names"), timevar = "Raw.file", direction = "wide")



data_unique <- make_unique(df.LFQ, "Gene.names", "Proteins", delim = ";")
colnames(data_unique)<-str_replace(colnames(data_unique),"-","_")
expdesign<-annot%>%select(Raw.file,Condition,BioReplicate)
colnames(expdesign)=c('label','condition','replicate')   
expdesign$label<-str_replace(expdesign$label,"-","_")
expdesign$replicate<-str_split_fixed(expdesign$label,"_",4)[,4]
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

TP<-df %>% filter(grepl(paste(proteins, collapse="|"),ID) & `p.adj` < .05) %>% nrow() #21
FP<-df %>% filter(!grepl(paste(proteins, collapse="|"),ID) & `p.adj` < .05) %>% nrow() #216


#####################################################################################################

####################proDA##########################################################################


input_data<-infile %>% filter(Potential.contaminant != "+") %>% select(Proteins, Raw.file, Intensity)

df.LFQ = as.data.frame(input_data) %>% select(Proteins,Raw.file, Intensity)
df.LFQ = df.LFQ %>% group_by(Proteins,Raw.file) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
df.LFQ = df.LFQ %>% filter(Proteins != "")
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "Proteins", timevar = "Raw.file", direction = "wide")


maxquant_protein_table<-df.LFQ

intensity_colnames <- grep("Intensity\\.", colnames(maxquant_protein_table), value=TRUE)

# Create the intensity matrix
abundance_matrix <- as.matrix(maxquant_protein_table[, intensity_colnames])
# Adapt column and row maxquant_protein_table
colnames(abundance_matrix) <- sub("Intensity\\.", "", intensity_colnames)
rownames(abundance_matrix) <- maxquant_protein_table$Proteins
colnames(abundance_matrix)<-str_replace(colnames(abundance_matrix),"-","_")

## Replace the zeros with NAs and take the log2() of the data

abundance_matrix[abundance_matrix == 0] <- NA
abundance_matrix <- log2(abundance_matrix)

## Normalize the data using median_normalization()
normalized_abundance_matrix <- median_normalization(abundance_matrix)

## Fit the probabilistic dropout model
sample_info_df <- data.frame(name = colnames(normalized_abundance_matrix),
                             stringsAsFactors = FALSE)
sample_info_df$condition <- str_split_fixed(sample_info_df$name,"_",4)[,3]
sample_info_df$replicate <- str_split_fixed(sample_info_df$name,"_",4)[,4]

normalized_abundance_matrix<-na.omit(normalized_abundance_matrix)
fit <- proDA(normalized_abundance_matrix, design = sample_info_df$condition, 
             col_data = sample_info_df)
result_names(fit)
# Test which proteins differ between condition 
test_res1 <- proDA::test_diff(fit,contrast ="sample1-sample2")
test_res2 <- proDA::test_diff(fit,contrast ="sample1-sample3")
test_res3 <- proDA::test_diff(fit,contrast ="sample1-sample4")
test_res4 <- proDA::test_diff(fit,contrast ="sample2-sample3")
test_res5 <- proDA::test_diff(fit,contrast ="sample2-sample4")
test_res6 <- proDA::test_diff(fit,contrast ="sample3-sample4")


test<-rbind(test_res1,test_res2,test_res3,test_res4,test_res5,test_res6)

TP<-test %>% filter(grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #20
FP<-test %>% filter(!grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #30