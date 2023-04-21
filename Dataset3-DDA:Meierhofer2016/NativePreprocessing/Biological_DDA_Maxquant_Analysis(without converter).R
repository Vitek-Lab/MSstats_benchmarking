#################Analysis of Biological(DDA:Meierhofer2016) for all 8 methods###########################
##Tool: Maxquant

## Using Native preprocessing

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



##############################
## Read MaxQuant report
##############################
# read in proteinGroups file, in order to use protein ids
proteinGroups<-read.table("C:/Users/tushi/OneDrive/Desktop/Data/Biological/proteinGroups.txt", sep="\t", header=TRUE)

# Read in MaxQuant file: evidence.txt
infile <- read.table("C:/Users/tushi/OneDrive/Desktop/Data/Biological/evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates
annot <- read.csv("C:/Users/tushi/OneDrive/Desktop/Data/Biological/experimentalDesign_quant_annotation.csv", header=TRUE)


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
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

##############################
## Model-based comparison + adjust p-value
##############################

comparison <- matrix(c(-1,1), nrow=1)
row.names(comparison) <- c("Liver-Cerebellum")
colnames(comparison)<-c('Liver','Cerebellum')
test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

# Results
Significant<-test.MSstats %>% filter(adj.pvalue < .05) %>% nrow() #2410



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
df.LFQ$Intensity = log2(df.LFQ$Intensity)
# run_median = df.LFQ %>% group_by(Run) %>% summarize(run_med = median(Intensity, na.rm=TRUE))
# overall_med = median(run_median$run_med)
# df.LFQ$Intensity = df.LFQ$Intensity + (overall_med - df.LFQ$Intensity)

df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(log2Intensity = log2(sum(2^Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

## Convert into DEqMS format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

# df.LFQ = medianSummary(df.LFQ, group_col = 1, ref_col=2)

# df.LFQ[2:13] = equalMedianNormalization(df.LFQ[2:13])
df.LFQ[sapply(df.LFQ, is.infinite)] <- NA


df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:13])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[14:25])))

# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 & df.LFQ$na_count_2<2, ]

# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count <- proteinGroups %>% select(Majority.protein.IDs, Unique.peptides)

# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
rownames(pep.count) = pep.count$Majority.protein.IDs
pep.count = pep.count %>% select(-Majority.protein.IDs)

## DEqMS analysis on LFQ data
rownames(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)

protein.matrix = as.matrix(df.LFQ.filter[1:24])

class = as.factor(rep(1:2, each=12))
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(class1-class2, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count[rownames(fit3$coefficients),]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error first
min(fit3$count)

fit4 = spectraCounteBayes(fit3)

## Analyze results
results_list = list()
for (i in seq_along(colnames(fit4$coefficients))){
  meier_deqms_results = outputResult(fit4,coef_col = i)
  meier_deqms_results$Label = colnames(fit4$coefficients)[[i]]
  results_list[[i]] = meier_deqms_results
}
meier_deqms_results = rbindlist(results_list)


Signficant<-meier_deqms_results %>% filter(adj.P.Val < .05) %>% nrow() #444












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
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

## Convert into DEqMS format
df.LFQ = dcast(df.LFQ, ProteinName + PSM~Run, 
               value.var = "Intensity", fun.aggregate = max, fill=NA)
colnames(df.LFQ)[3:26] = paste("Intensity", colnames(df.LFQ)[3:26], sep="_")
setcolorder(df.LFQ, c("ProteinName", "PSM", "Intensity_WT1_CB_rep1", "Intensity_WT1_CB_rep2", 
                      "Intensity_WT2_CB_rep1", "Intensity_WT2_CB_rep2", 
                      "Intensity_WT3_CB_rep1", "Intensity_WT3_CB_rep2", 
                      "Intensity_WT4_CB_rep1", "Intensity_WT4_CB_rep2", 
                      "Intensity_WT5_CB_rep1", "Intensity_WT5_CB_rep2", 
                      "Intensity_WT6_CB_rep1", "Intensity_WT6_CB_rep2", 
                      "Intensity_WT1_liver_rep1", "Intensity_WT1_liver_rep2", 
                      "Intensity_WT2_liver_rep1", "Intensity_WT2_liver_rep2", 
                      "Intensity_WT3_liver_rep1", "Intensity_WT3_liver_rep2", 
                      "Intensity_WT4_liver_rep1", "Intensity_WT4_liver_rep2", 
                      "Intensity_WT5_liver_rep1", "Intensity_WT5_liver_rep2", 
                      "Intensity_WT6_liver_rep1", "Intensity_WT6_liver_rep2"))

## Extract columns of LFQ intensites
df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[3:14])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[15:26])))

# Filter protein table. require minimum two values for each group.
df.LFQ = df.LFQ[df.LFQ$na_count_1<11 & df.LFQ$na_count_2<11,]
ecols <- grep("Intensity", colnames(df.LFQ))

pe = readQFeatures(df.LFQ, fnames = 2, ecol = ecols,  name = "peptideRaw")

# cond <- which((colnames(pe)[[1]][1], split = "_")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- rep(c("1", "2"), each=12) %>% as.factor()
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

L <- makeContrast("condition2=0", parameterNames = c("(Intercept)",
                                                     "condition2"))

pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

msqrob_meier_results = rowData(pe[["protein"]])$condition2
msqrob_meier_results$Protein=rownames(msqrob_meier_results)

Significant<-msqrob_meier_results %>% filter(adjPval < .05) %>% nrow() #854


################################# PMartR ###################################################################



e_mdata<-infile %>% filter(Potential.contaminant != "+") %>% select(Sequence, Raw.file, Intensity)
e_mdata<-e_mdata%>%group_by(Sequence,Raw.file)%>% summarise(Intensity=sum(Intensity))
e_mdata<-e_mdata %>% pivot_wider(names_from = Raw.file, values_from = Intensity,values_fill = NA)

vec<-data.frame(c=colnames(e_mdata))
f_mdata<-annot[match(vec[2:25,],annot$Raw.file),]
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


Significant<-m %>% filter( `p-value` < .05) %>% nrow()  #1067

############################################## DEP #################################################################




input_data<-infile %>% filter(Potential.contaminant != "+") %>% select(Proteins,Gene.names, Raw.file, Intensity)

df.LFQ = as.data.frame(input_data) %>% select(Proteins, Gene.names,Raw.file, Intensity)
df.LFQ = df.LFQ %>% group_by(Proteins,Gene.names,Raw.file) %>% summarise(Intensity = sum(Intensity, na.rm = TRUE))

df.LFQ = df.LFQ %>% filter(Proteins != "")
df.LFQ = df.LFQ %>% dplyr::rename(ProteinName=Proteins, Run = Raw.file)


## Take logsum of proteins
df.LFQ = as.data.frame(df.LFQ) %>% select(ProteinName,Gene.names, Run, Intensity)
df.LFQ$Intensity = log2(df.LFQ$Intensity)
df.LFQ = df.LFQ %>% group_by(ProteinName, Gene.names,Run) %>% summarise(log2Intensity = log2(sum(2^Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

df.LFQ[sapply(df.LFQ, is.infinite)] <- NA


## Convert into Intensity format
df.LFQ = reshape(df.LFQ, idvar = c("ProteinName","Gene.names"), timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))


data_unique <- make_unique(df.LFQ, "Gene.names", "ProteinName", delim = ";")

expdesign<-annot%>%select(Raw.file,Condition,BioReplicate)
colnames(expdesign)=c('label','condition','replicate')   

expdesign$replicate<-paste0(expdesign$replicate,"_",str_split_fixed(expdesign$label,"_",4)[,3])
columns <- grep("Intensity", colnames(data_unique))
data<-make_se(data_unique,columns,expdesign)




# Filter, normalize and impute missing values
filt <- filter_missval(data, thr = 0)


norm <- normalize_vsn(filt)
imputed <- DEP::impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- DEP::test_diff(imputed, type="manual",test=c("Liver_vs_Cerebellum"))
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)

# Generate a results table
data_results <- get_results(dep)
#data_results

###############results
df<-select(data_results,c("ID",(matches("p.adj"))))
df<-df %>% pivot_longer(cols=c((matches("p.adj"))),
                        names_to='label',
                        values_to='p.adj')

Significant<-df %>% filter( `p.adj` < .05) %>% nrow() #0


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
sample_info_df$condition <- str_split_fixed(sample_info_df$name,"_",3)[,2]
sample_info_df$replicate <- paste0(str_split_fixed(sample_info_df$name,"_",3)[,1],"_",str_split_fixed(sample_info_df$name,"_",3)[,3])

normalized_abundance_matrix<-na.omit(normalized_abundance_matrix)
fit <- proDA(normalized_abundance_matrix, design = sample_info_df$condition, 
             col_data = sample_info_df)
result_names(fit)
# Test which proteins differ between condition 
test <- proDA::test_diff(fit,contrast ="liver-CB")


Significant<-test %>% filter(adj_pval < .05) %>% nrow() #601

##############################################Prolfqua #########################################################


startdata<-infile%>%group_by(Proteins,Raw.file)%>%summarise(Intensity=sum(Intensity,na.rm = TRUE))

startdata <- dplyr::inner_join(annot, startdata, by = "Raw.file")


atable <- AnalysisTableAnnotation$new()

atable$fileName = "Raw.file"
atable$hierarchy[["protein_Id"]] <- c("Proteins")
atable$factors[["dilution."]] = "Condition"
atable$set_response("Intensity")
config <- AnalysisConfiguration$new(atable)
adata <- setup_analysis(startdata, config)




lfqpep<- prolfqua::LFQData$new(adata, config)

lfqpep$remove_small_intensities()$filter_proteins_by_peptide_count()

lfqpep<-lfqpep$get_Transformer()$log2()$robscale()$lfq



formula_Condition <-  strategy_lm("transformedIntensity ~ dilution.")

contrasts<-c("Liver-Cerebellum"="dilution.Liver-dilution.Cerebellum")

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

Significant<- m%>%filter(`p-value` < .05) %>% nrow()#3016

#####################################Limma#############################################

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
df.LFQ$Intensity = log2(df.LFQ$Intensity)
# run_median = df.LFQ %>% group_by(Run) %>% summarize(run_med = median(Intensity, na.rm=TRUE))
# overall_med = median(run_median$run_med)
# df.LFQ$Intensity = df.LFQ$Intensity + (overall_med - df.LFQ$Intensity)

df.LFQ = df.LFQ %>% group_by(ProteinName, Run) %>% summarise(log2Intensity = log2(sum(2^Intensity, na.rm = TRUE)))
df.LFQ = as.data.frame(df.LFQ)

## Convert into DEqMS format
df.LFQ = reshape(df.LFQ, idvar = "ProteinName", timevar = "Run", direction = "wide")
df.LFQ = df.LFQ %>% filter(!grepl(";", ProteinName))

# df.LFQ = medianSummary(df.LFQ, group_col = 1, ref_col=2)

# df.LFQ[2:13] = equalMedianNormalization(df.LFQ[2:13])
df.LFQ[sapply(df.LFQ, is.infinite)] <- NA


df.LFQ$na_count_1 = apply(df.LFQ,1,function(x) sum(is.na(x[2:13])))
df.LFQ$na_count_2 = apply(df.LFQ,1,function(x) sum(is.na(x[14:25])))

# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_1<2 & df.LFQ$na_count_2<2, ]



## DEqMS analysis on LFQ data
rownames(df.LFQ.filter) = df.LFQ.filter$ProteinName
df.LFQ.filter = df.LFQ.filter %>% select(-ProteinName)

protein.matrix = as.matrix(df.LFQ.filter[1:24])

class = as.factor(rep(1:2, each=12))
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(class1-class2, levels = design)
fit2 = contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2)

m<-as.data.frame(fit3$p.value)
m$protein<-rownames(fit3$p.value)
m<-m %>% pivot_longer(cols=starts_with('class'),names_to='label',values_to='p-value')
m$adj.pvalue = p.adjust(m$`p-value`,method = 'BH')


Significant<- m%>%filter(adj.pvalue < .05) %>% nrow()#444


##################################################################################################
# Graphical reprentation of all methods with their number of significant proteins

data<-data.frame(methods=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),
                 values=c(2410,444,854,1067,0,601,3016,444))


## Method Comp
data %>% mutate(Version = factor(methods, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=values, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  labs(y = "Number of Significant Proteins Identified",
       x="", title="Dataset 3: Biological Mixture Meierhofer2016 - DDA - MaxQuant") + theme(axis.text.y=element_blank(),
                                                                             axis.title.y=element_text(size=10, face='bold'),
                                                                             axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  
  annotate("text", x = c(1), y=550, label = c("2410"), fontface =2, size=6) +
  annotate("text", x = c(2), y=550, label = c("444"), fontface =2, size=6) +
  annotate("text", x = c(3), y=550, label = c("854"), fontface =2, size=6)+
  annotate("text", x = c(4), y=550, label = c("1067"), fontface =2, size=6) +
  annotate("text", x = c(5), y=90, label = c("0"), fontface =2, size=6) +
  annotate("text", x = c(6), y=550, label = c("601"), fontface =2, size=6) +
  annotate("text", x = c(7), y=550, label = c("3016"), fontface =2, size=6)  +
  annotate("text", x = c(8), y=550, label = c("444"), fontface =2, size=6)  