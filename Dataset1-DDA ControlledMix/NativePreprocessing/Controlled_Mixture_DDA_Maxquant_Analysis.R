#################Analysis of DDA Controlled Mixture for all 8 methods###########################
###Tool: MaxQuant
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

# read in proteinGroups file, in order to use protein ids
proteinGroups <- read.table("C:/Users/tushi/OneDrive/Desktop/Data/ControlMixture_DDA_MaxQuant_proteinGroups.txt", sep="\t", header=TRUE)

# Read in MaxQuant file: evidence.txt
infile <- read.table("C:/Users/tushi/OneDrive/Desktop/Data/ControlMixture_DDA_MaxQuant_evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates
annot <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/ControlMixture_DDA_MaxQuant_annotation.csv", header=TRUE)

###########################MSstats#########################################################################

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

colnames(comparison)<-c('1','2','3','4','5')
test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

# Calculating TP and FP
TP<-test.MSstats %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adj.pvalue < .05) %>% nrow() #191
FP<-test.MSstats %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& adj.pvalue < .05) %>% nrow() #77


#######################DeqMS##################################################################################


df.prot<-proteinGroups
# remove decoy matches and matches to contaminant
df.prot = df.prot[!df.prot$Reverse=="+",]
# Extract columns of LFQ intensites
df.LFQ = df.prot[,128:142]
df.LFQ[df.LFQ==0] <- NA

rownames(df.LFQ) = df.prot$Majority.protein.IDs

df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[1:3])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[4:6])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[7:9])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[10:12])))
df.LFQ$na_count_5 = apply(df.LFQ,1,function(x) sum(is.na(x[13:15])))

# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                         df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2|df.LFQ$na_count_5<2,1:15]



# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count.table = data.frame(count = rowMins(as.matrix(df.prot[,28:42])),
                             row.names = df.prot$Majority.protein.IDs)
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count.table$count = pep.count.table$count+1


protein.matrix = na.omit(log2(as.matrix(df.LFQ.filter)))

class = as.factor(rep(1:5, each=3))
class=as.character(class)
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
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
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]

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

TP<-iprg_deqms_results %>% filter(grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow()#196
FP<-iprg_deqms_results %>% filter(!grepl(paste(proteins, collapse="|"),gene) & adj.P.Val < .05) %>% nrow()#80

###########################################MSqRob2#########################################################

evidence=infile
input_data = evidence %>% filter(Reverse != "+") %>% select(Proteins, Sequence, Charge, Raw.file, Intensity)
input_data$PSM = paste(input_data$Sequence, input_data$Charge, sep="_")
input_data = input_data %>% group_by(Proteins, Raw.file, PSM) %>% summarise(Intensity=max(Intensity))
input_data = input_data %>% filter(Proteins != "")
input_data = input_data %>% dplyr::rename(ProteinName=Proteins, Run = Raw.file)

## Take logsum of proteins
df.LFQ = as.data.frame(input_data) %>% select(ProteinName, PSM, Run, Intensity)
df.LFQ = as.data.table(df.LFQ)

## Convert into MSqRob2 format
df.LFQ = dcast(df.LFQ, ProteinName + PSM~Run, 
               value.var = "Intensity", fun.aggregate = max, fill=NA)
colnames(df.LFQ)[3:17] = paste("Intensity", colnames(df.LFQ)[3:17], sep="_")

## Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[3:5])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[6:8])))
df.LFQ$na_count_3 = apply(df.LFQ,1,function(x) sum(is.na(x[9:11])))
df.LFQ$na_count_4 = apply(df.LFQ,1,function(x) sum(is.na(x[12:14])))
df.LFQ$na_count_5 = apply(df.LFQ,1,function(x) sum(is.na(x[15:17])))
# Filter protein table. require minimum two values for each group.
df.LFQ = df.LFQ[df.LFQ$na_count_1<2 | df.LFQ$na_count_2<2 |
                  df.LFQ$na_count_3<2 | df.LFQ$na_count_4<2| df.LFQ$na_count_5<2,]
ecols <- grep("Intensity", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 2, ecol = ecols,  name = "peptideRaw")

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
L10 <- makeContrast("condition5-condition3=4", parameterNames = c("(Intercept)",
                                                                  "condition2",
                                                                  "condition3",
                                                                  "condition4","condition5"))
L=cbind(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10)

pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

results_list = list()
for (i in 4:9){
  df = rowData(pe[["protein"]])[[i]]
  print(df)
  iprg_msqrob_results = tibble::rownames_to_column(df, "Protein")
  iprg_msqrob_results$Label = colnames(rowData(pe[["protein"]]))[[i]]
  results_list[[i]] = iprg_msqrob_results
}
iprg_msqrob_results = rbindlist(results_list)

TP<-iprg_msqrob_results %>% filter(grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #123
FP<-iprg_msqrob_results %>% filter(!grepl(paste(proteins, collapse="|"),Protein) & adjPval < .05) %>% nrow() #85

################################################PMartR######################################################


e_mdata<-infile %>% select(Sequence, Raw.file, Intensity)
e_mdata<-e_mdata%>%group_by(Sequence,Raw.file)%>% summarise(Intensity=sum(Intensity))
e_mdata<-e_mdata %>% pivot_wider(names_from = Raw.file, values_from = Intensity,values_fill = NA)


# Processing to convert the data into pep data object
f<-data.frame(i=c('sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12','sample_13','sample_14','sample_15'))

vec<-data.frame(c=colnames(e_mdata))
f_mdata<-annot[match(vec[2:16,],annot$Raw.file),]
f_mdata<-cbind(f,f_mdata)
test<-infile[,c('Sequence','Proteins')]
test %>%group_by(Sequence) %>%
  fill(Proteins, .direction = "downup") %>%
  ungroup %>%
  distinct
reorder_idx <- match(e_mdata$Sequence,test$Sequence)
test <- test[match(e_mdata$Sequence,test$Sequence),]  

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



##Storing result of PmartR in a tv=ble
m<-mystats %>% select(Protein,starts_with('P_Value'))
m<-m %>% pivot_longer(cols=starts_with('P_Value'),names_to='label',values_to='p-value')


TP<-m %>% filter(grepl(paste(proteins, collapse="|"),Protein) & `p-value` < .05) %>% nrow() #117
FP<-m %>% filter(!grepl(paste(proteins, collapse="|"),Protein)& `p-value` < .05) %>% nrow() #197


###########################################DEP####################################################################


data<-select(proteinGroups,c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names",
                             "Fasta.headers","Peptides","Razor...unique.peptides","Unique.peptides",
                             "Only.identified.by.site","Reverse","Potential.contaminant",starts_with("LFQ")))

expdesign<-select(annot,c('Raw.file',"Condition","BioReplicate"))

expdesign$Raw.file=str_split_fixed(expdesign$Raw.file, "_01_", 2)[,2]
colnames(expdesign)=c('label','condition','replicate')               

data<-import_MaxQuant(data,expdesign)



# Filter, normalize and impute missing values
filt <- filter_missval(data, thr = 0)


norm <- normalize_vsn(filt)
imputed <- DEP::impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- DEP::test_diff(imputed, type="all")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Generate a results table
data_results <- get_results(dep)
data_results

###############results
df<-select(data_results,c("ID",(matches("p.adj"))))
df<-df %>% pivot_longer(cols=c((matches("p.adj"))),
                        names_to='label',
                        values_to='p.adj')

TP<-df %>% filter(grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #199
FP<-df %>% filter(!grepl(paste(proteins, collapse="|"),ID) & p.adj < .05) %>% nrow() #92
#####################################################################################################

####################proDA#################################################################################



maxquant_protein_table<-proteinGroups

intensity_colnames <- grep("^LFQ\\.intensity\\.", colnames(maxquant_protein_table), value=TRUE)
# Create the intensity matrix
abundance_matrix <- as.matrix(maxquant_protein_table[, intensity_colnames])
# Adapt column and row maxquant_protein_table
colnames(abundance_matrix) <- sub("^LFQ\\.intensity\\.", "", intensity_colnames)
rownames(abundance_matrix) <- maxquant_protein_table$Protein.IDs


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

TP<-test %>% filter(grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #160
FP<-test %>% filter(!grepl(paste(proteins, collapse="|"),name) & adj_pval < .05) %>% nrow() #72

#####################################################################################################

####################prolfqua#################################################################################



startdata <- prolfqua::tidyMQ_ProteinGroups(proteinGroups)
annot$raw.file<-tolower(str_split_fixed(annot$Raw.file,"_",5)[,5])
annot$Condition<-as.character(annot$Condition)
startdata <- dplyr::inner_join(annot, startdata, by = "raw.file")

startdata <- dplyr::filter(startdata, nr.peptides > 1)
atable <- AnalysisTableAnnotation$new()

atable$fileName = "raw.file"
atable$hierarchy[["protein_Id"]] <- c("proteinID")
atable$factors[["dilution."]] = "Condition"
atable$set_response("mq.protein.intensity")
config <- AnalysisConfiguration$new(atable)
adata <- setup_analysis(startdata, config)


lfqpep<- prolfqua::LFQData$new(adata, config)

lfqpep$remove_small_intensities()$filter_proteins_by_peptide_count()

lfqpep<-lfqpep$get_Transformer()$log2()$robscale()$lfq



formula_Condition <-  strategy_lm("transformedIntensity ~ dilution.")

contrasts<-c("1_vs_2"="dilution.2-dilution.1",
             "1_vs_3"="dilution.3-dilution.1",
             "1_vs_4"="dilution.4-dilution.1",
             "1_vs_5"="dilution.5-dilution.1",
             "2_vs_3"="dilution.3-dilution.2",
             "2_vs_4"="dilution.4-dilution.2",
             "2_vs_5"="dilution.5-dilution.2",
             "3_vs_4"="dilution.4-dilution.3",
             "3_vs_5"="dilution.5-dilution.3",
             "4_vs_5"="dilution.5-dilution.4")

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


m<-results %>% select(protein_Id,starts_with('FDR'))

m<-m %>% pivot_longer(cols=starts_with('FDR'),names_to='label',values_to='p-value')

TP<- m%>%filter(grepl(paste(proteins, collapse="|"),protein_Id) & `p-value` < .05) %>% nrow()#185
FP<-m%>%filter(!grepl(paste(proteins, collapse="|"),protein_Id) & `p-value` < .05) %>% nrow()#124




#####################################################################################################

####################limma#################################################################################

df.prot<-proteinGroups
# remove decoy matches and matches to contaminant
df.prot = df.prot[!df.prot$Reverse=="+",]
# Extract columns of LFQ intensites
df.LFQ = df.prot[,128:142]
df.LFQ[df.LFQ==0] <- NA

rownames(df.LFQ) = df.prot$Majority.protein.IDs

protein.matrix = na.omit(log2(as.matrix(df.LFQ)))

class = as.factor(rep(1:5, each=3))
class=as.character(class)
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
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
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

proteins<-c("P02701","P00711","Q29443","Q29550","P0CG53",
            "P68082","P00432","P02754","P24627","P80025",
            "P00915","P02787","P02663","P01008","P00921",
            "P05307","P61769","P02662","P01012","P02666",
            "Q3SX14","P00563","P02769","Q58D62","P00698",
            "P00004","P00711","P00442","P01133","P02753")

m<-as.data.frame(fit3$p.value)
m$protein<-rownames(fit3$p.value)
m<-m %>% pivot_longer(cols=starts_with('class'),names_to='label',values_to='p-value')
m$adj.pvalue = p.adjust(m$`p-value`,method = 'BH')


TP<- m%>%filter(grepl(paste(proteins, collapse="|"),protein) & adj.pvalue < .05) %>% nrow()#196
FP<-m%>%filter(!grepl(paste(proteins, collapse="|"),protein) & adj.pvalue < .05) %>% nrow()#83
