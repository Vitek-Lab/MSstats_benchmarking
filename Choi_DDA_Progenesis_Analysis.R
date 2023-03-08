#################Analysis of DDA choi for all 6 methods###########################
##Tool: Skyline

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
raw <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/DDA_choi_progenesis/Choi2017_DDA_Progenesis_input.csv", stringsAsFactors=F) # the data file

annot <- read.csv('C:/Users/tushi/OneDrive/Desktop/Data/DDA_choi_progenesis/Choi2017_DDA_Progenesis_annotation.csv')



##################################MSstats########################################################
##############################
## Make MSstats required format
##############################
quant <- ProgenesistoMSstatsFormat(raw, annotation=annot, 
                                   removeProtein_with1Peptide = TRUE)

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
TP<-test.MSstats %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #22
FP<-test.MSstats %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& adj.pvalue < .05) %>% nrow() #110


##################################### DeqMS###################################################################

quant$Run=gsub('[-]','_',quant$Run)
annot$Run=gsub('[-]','_',annot$Run)
input_data=quant
input_data$ProteinName=data.frame(do.call("rbind", strsplit(quant$ProteinName,split="|",fixed=3)))$X2



## Just use MSstats converter to make life easier
# input_data = MaxQtoMSstatsFormat(evidence, annotation, df.prot, removeFewMeasurements = FALSE)

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(log2Intensity = log2(sum(Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

## Convert into Intensity format
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



# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                         df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2,]


row.names(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)



# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count <- as.data.frame(input_data )%>% select(ProteinName, PeptideModifiedSequence)
pep.count<-pep.count%>%group_by(ProteinName)%>%
  summarise(count=n_distinct(PeptideModifiedSequence),.groups = 'drop')%>%
  as.data.frame()
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count$count = pep.count$count+1
rownames(pep.count) = pep.count$ProteinName
pep.count = pep.count %>% select(-ProteinName)



## DEqMS analysis on LFQ data
protein.matrix <- na.omit(as.matrix(df.LFQ.filter[1:12]))

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
# iprg_deqms_results$adj.P.Val = p.adjust(iprg_deqms_results$P.Value, method="BH")
# iprg_deqms_results$sca.adj.pval = p.adjust(iprg_deqms_results$sca.P.Value, method="BH")

TP<-iprg_deqms_results %>% filter(grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #24
FP<-iprg_deqms_results %>% filter(!grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow() #204

############################################################################################

input_data=quant
input_data$ProteinName=data.frame(do.call("rbind", strsplit(quant$ProteinName,split="|",fixed=3)))$X2
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ$PSM = paste(quant$PeptideSequence, quant$PrecursorCharge, sep="_")
df.LFQ = as.data.table(df.LFQ)

## Convert into MSqRob2 format
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

pe = readQFeatures(df.LFQ, fnames = 1, ecol = ecols,  name = "peptideRaw")

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
for (i in 8:13){
  df = rowData(pe[["protein"]])[[i]]
  iprg_msqrob_results = tibble::rownames_to_column(df, "Protein")
  iprg_msqrob_results$Label = colnames(rowData(pe[["protein"]]))[[i]]
  results_list[[i]] = iprg_msqrob_results
}
iprg_msqrob_results = rbindlist(results_list)


TP<-iprg_msqrob_results %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #24
FP<-iprg_msqrob_results %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #254



########################################## PMartR #################################################################



e_mdata<-quant[,c('PeptideModifiedSequence','Run','Intensity')]
e_mdata<-e_mdata%>%group_by(PeptideModifiedSequence,Run)%>% summarise(Intensity=sum(Intensity))
e_mdata<-e_mdata %>% pivot_wider(names_from = Run, values_from = Intensity,values_fill = NA)


vec<-data.frame(c=colnames(e_mdata))
annot$Run=gsub('[.]', '', annot$Run)
f_mdata<-annot[match(vec[2:13,],annot$Run),]

test<-quant[,c('PeptideModifiedSequence','ProteinName')]
test %>%group_by(PeptideModifiedSequence) %>%
  fill(ProteinName, .direction = "downup") %>%
  ungroup %>%
  distinct
reorder_idx <- match(e_mdata$PeptideModifiedSequence,test$PeptideModifiedSequence)
test <- test[match(e_mdata$PeptideModifiedSequence,test$PeptideModifiedSequence),]  

colnames(test)<-c("PeptideModifiedSequence","Protein")

mypepData <- as.pepData(e_data=e_mdata, f_data = f_mdata, e_meta=test,edata_cname = "PeptideModifiedSequence", emeta_cname="Protein",fdata_cname = "Run", data_scale = "abundance",check.names=FALSE)
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


TP<-m %>% filter(grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #28
FP<-m %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #413

#####################DEP######################################################################


input_data=quant

df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")



data_unique <- make_unique(df.LFQ, "ProteinName", "ProteinName", delim = ";")


expdesign<-select(annot,c('Run',"Condition","BioReplicate"))


colnames(expdesign)=c('label','condition','replicate')  
expdesign<-expdesign[order(expdesign$condition),]
expdesign$replicate<-rep(1:3,4)

columns <- grep("Intensity", colnames(data_unique))
data<-make_se(data_unique,columns,expdesign)



# Filter, normalize and impute missing values
filt <- filter_missval(data, thr = 0)


norm <- normalize_vsn(filt)
imputed <- DEP::impute(norm, fun = "MinProb")

# Test for differentially expressed proteins
diff <- DEP::test_diff(imputed, "all")


dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Generate a results table
data_results <- get_results(dep)
#data_results

###############results
df<-select(data_results,c("ID",(matches("p.adj"))))
df<-df %>% pivot_longer(cols=c((matches("p.adj"))),
                        names_to='label',
                        values_to='p.adj')

TP<-df %>% filter(grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow()  #28
FP<-df %>% filter(!grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow()  #413


##########################################  ProDA  ############################################################

input_data=quant

df.LFQ = as.data.frame(input_data) %>% select(ProteinName, Run, Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))
df.LFQ = as.data.frame(df.LFQ)
## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")

protein_table<-df.LFQ

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

TP<-test %>% filter(grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #25
FP<-test %>% filter(!grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #219

######################################################################################################