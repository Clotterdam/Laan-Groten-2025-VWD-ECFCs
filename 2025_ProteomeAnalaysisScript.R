
#Scipt for analysis proteome data corresponding to manuscript:
#"Von willebrand disease specific defects and proteomic signatures in endothelial colony forming cells"
#script author: Stijn Groten, date: 07-01-2025


# libraries
{
library(tidyverse)
library(limma)
library(sva)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
library(splines)
library(WGCNA) 
library(plyr)
library(dplyr)
library(scales)
library(EnhancedVolcano)
library(UpSetR)
library(eulerr)
library(circlize)
library(umap)
library(corrplot)
library(ggvenn)
library(cluster)
}


# setseed
set.seed(1221)

### setworking directory
getwd()
setwd("C:\\")
setwd("C:\\Users\\grote\\OneDrive - Sanquin\\Experiments\\20231113_VWF_Bas/240417_VWFbas_output/ScriptGithub/")

#load raw data
RAW = read.csv('Proteomedata.csv', header=T, stringsAsFactors = F, dec=',')

# clean up
prot_clean = RAW #subset(RAW, Potential.contaminant != '+' & Only.identified.by.site !='+' & Reverse !='+')

#pick columns with intensities
prot_clean_small <- RAW %>% dplyr::select(contains(c("Desktop")))

###rename_columns
colnames(prot_clean_small) <- gsub("C..Users.massspecuser.Desktop.Stijn.240415_VWF_fin.RawFiles.231128_SG_BL_VWF_"  ,"", colnames(prot_clean_small))
prot_clean_small <- prot_clean_small %>% dplyr::rename_with(~stringr::str_remove(., '_S1.*$'))

#create protein info file
prot_info = prot_clean %>% dplyr::select(-contains(c("Desktop", "X")))

#LFQ only file
prot_clean_small<- as.data.frame(sapply(prot_clean_small, as.numeric))

# replace 0 with NA and log2
prot_clean_small[prot_clean_small==0] = NA
prot_clean_small = log2(prot_clean_small)
rownames(prot_clean_small) <- prot_info$Protein.Group

#check LFQ distribution across samples
boxplot(prot_clean_small)

#create pheno with sample information
pheno <- read.csv("pheno.csv", header = T)

#heatmap sample correlations
prot_clean_small_NAomit <-  na.omit(prot_clean_small)
M<- cor(prot_clean_small_NAomit)
col_fun = colorRamp2(c(0.85, 1), c("white", "red"))
Heatmap(M, cluster_columns = T, cluster_rows = T, col = col_fun)


# Calculate PCA
pci <- prot_clean_small %>% 
  na.omit() %>%
  t() %>%
  prcomp(rank=10)


# Plot PCA
## Scores plot
pci$x %>%
  as.data.frame() %>%
  rownames_to_column( "sample") %>%
  left_join(pheno, "sample") %>%
  ggplot(aes(x = PC1, y = PC2, col = TYPE, label = sample)) +
  geom_point(size = 5) +
  theme_bw() + 
  geom_text(nudge_y = 0, col = "black", size = 2) +
  scale_shape_manual(values = c("no VWf" = 21, "no DDAVP" = 22, "Control" = 23)) +
  labs(
    title = "Principal component analysis",
    subtitle = "Scores plot",
    x=paste0("PC1: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,1]*100,1),"%"),
    y=paste0("PC2: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,3]*100,1),"%"))

## Loading plot
pci$rotation %>%
  as.data.frame() %>%
  rownames_to_column("Protein.Group") %>%
  left_join(prot_info, by = "Protein.Group") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(data = . %>% slice_max(order_by = c(abs(PC1)), n = 15), aes(label = Genes), show.legend = F, max.overlaps = Inf) +
  geom_point(data = . %>% slice_max(order_by = c(abs(PC1)), n = 15), shape = 21,size = 3, fill = 'red4', show.legend = F) +
  geom_text_repel(data = . %>% slice_max(order_by = c(abs(PC2)), n = 15), aes(label = Genes), show.legend = F, max.overlaps = Inf) +
  geom_point(data = . %>% slice_max(order_by = c(abs(PC2)), n = 15), shape = 21,size = 3, fill = 'red4', show.legend = F) +
  theme_bw() +
  labs(
    title = "Principal component analysis",
    subtitle = "Loadings plot",
    x=paste0("PC1: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,1]*100,1),"%"),
    y=paste0("PC2: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,2]*100,1),"%"))


## VALID VALUES PER CONDITION
vvs <- prot_clean_small %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  t() %>%
  as.data.frame() %>%
  add_column(sample = pheno$sample)  %>%
  replace(is.na(.),0) %>%
  group_by(sample) %>%
  summarise_all(sum) %>%
  t() %>%
  as.data.frame()

#create filter column
names(vvs) <- lapply(vvs[1, ], as.character)
vvs_1 <- vvs[-1,] 
vvs_1 <- as.data.frame(sapply(vvs_1, as.numeric))
rownames(vvs_1) <- prot_info$Protein.Group
vvs_1 <- vvs_1 %>% dplyr::mutate(No_pass = rowSums(.)) %>% # 
  rownames_to_column(var = "Protein.Group")

# valid value filtering
df_prot_valid = prot_clean_small %>% 
  rownames_to_column('Protein.Group') %>%
  dplyr::select(pheno$sample, Protein.Group) %>%
  subset(Protein.Group %in% vvs_1$Protein.Group[vvs_1$No_pass >=6]) # valid treshold is 6 pep ids over all samples

#name non-imputed columns
colnames(df_prot_valid) <- paste(colnames(df_prot_valid),"_NONIMP")
colnames(df_prot_valid) <- gsub(" "  ,"", colnames(df_prot_valid))
colnames(df_prot_valid) <- gsub("Protein.Group_NONIMP"  ,"Protein.Group", colnames(df_prot_valid))

#save non imputed protein data
write.csv(df_prot_valid, file = "prot_valid.csv")


#calculate quantified rot numbers per sample
ProtNum <- prot_clean_small %>% select(-contains(c("Protein"))) %>% summarise_all(funs(sum(!is.na(.)))) %>%
  melt()
colnames(ProtNum) <- c("sample", "value")
ProtNum <- left_join(ProtNum, pheno, by = "sample")

ggplot(ProtNum, aes(x=sample, y=value, fill = CLUSTER_PROT)) + 
  geom_bar(stat='identity', position='identity') +
  coord_cartesian(ylim=c(0,9000))+ 
  scale_y_continuous(limits=c(0,9000), breaks = seq(0,9000,by=1000), expand=c(0,0))+
  xlab('')+ ylab('# Proteins')+
  theme_classic()+theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5))


#imputation of data
dat= NULL
dat = df_prot_valid
rownames(dat) = dat$Protein.Group
dat <- subset(dat, select = -c(Protein.Group))

#distrubtion check before imputation
boxplot(dat)

# script to get random numbers from a normal distribution for imputation in analysis
mann_downshift = 1.8
mann_width = 0.3

# make a loop for each column
for(i in 1:ncol(dat)) {    
  
  m = mean(na.omit(dat[,i])) - mann_downshift*sd(na.omit(dat[,i]))
  s= sd(na.omit(dat[,i]))*mann_width
  
  imputed_values = rnorm(dim(dat)[1],m,s)
  
  mask = is.na(dat[,i])
  dat[mask,i] = imputed_values[mask]
}

# plot boxplots to see if normalization is OK
boxplot(dat)

# replace data in DF with imputed data      
prot_valid_imp = dat
colnames(prot_valid_imp) <- gsub("_NONIMP","_IMP", colnames(prot_valid_imp))

#save imputed protein data
write.csv(prot_valid_imp, file = "prot_valid_imp.csv")


# check how sampels are distributed by PCA
colnames(prot_valid_imp)
#  # PCA on imputed data
pca = prcomp(t(prot_valid_imp))
# print the results
summary(pca) 
# make a df with bath info, samples info and PCA results
pci = data.frame(pca$x, CLUSTER_PROT = pheno$CLUSTER_PROT, TYPE = pheno$TYPE
                 )
# plot PC1 vs PC2
ggplot(pci, aes(x=PC1, y=PC2, color=CLUSTER_PROT, fill = CLUSTER_PROT, shape = TYPE)) + geom_point(size=4) + 
    theme(panel.grid.major = element_blank())  



#DETERMINE clustering on imputed data
df_cluster <- prot_valid_imp

#FUN =  kmeans / pam , K MAX = number of tests, B = number of iterations
gap=NULL
gap = clusGap(t(df_cluster), FUN= pam, K.max = 20, B=100)
plot(gap)
k = maxSE(gap$Tab[,'gap'], gap$Tab[,'SE.sim'], method ='Tibs2001SEmax') 

#set k based on outcome above
k=3
#add clustering to pheno (this case already done)
#pheno$CLUSTER_PROT = pam(t(df_cluster),3)$cluster


#Statistical analysis
#set groups to be tested as factors

#cluster_prot -> for teting different proteome clusters
pheno$CLUSTER_PROT <- as.factor(pheno$CLUSTER_PROT)
#nocluster -> for testing disease type without clustering
pheno$TYPE <- as.factor(pheno$TYPE)
#cluster_RNA_TYPE -> for testing disease type with qPCR based clustering
pheno$ClusterRNA_TYPE <- as.factor(pheno$ClusterRNA_TYPE)
#cluster_prot_type -> for testing disease type with protome based clustering
pheno$CLUSTER_PROT_TYPE <- as.factor(pheno$CLUSTER_PROT_TYPE)
#individual -> for testing individual changes
pheno$ID_TYPE <- as.factor(pheno$ID_TYPE)


#make a model
ct  = as.numeric(pheno$CLUSTER_PROT_TYPE) #change pheno$XXXXX to testing choice
design = model.matrix(~0+as.factor(ct))
colnames(design) = levels(pheno$CLUSTER_PROT_TYPE) #change pheno$XXXXX to testing choice

# fit in limma
fit = lmFit(prot_valid_imp, design)

# define comparisons with a contrast matrix
# define all possible combinations
contrast_list = combn(colnames(design), 2)
contrasts = NULL
for(i in seq(1, length(contrast_list), 2)){
  contrasts = c(contrasts, sprintf("%s - %s", contrast_list[i], contrast_list[i + 1]))
}

#select combinations of choice based on tests
contrasts
contrasts_select <- contrasts[c(grep(contrasts, pattern = "A_Control - A_no_DDAVP|A_Control - A_no_VWF|B_Control - B_no_DDAVP|B_Control - B_no_VWF|C_Control - C_no_DDAVP"))]
contrasts_select

#make contrasts matrix
contrast_matrix = makeContrasts(contrasts = contrasts_select,
                                levels = design)
# fit contrasts
fit2 = contrasts.fit(fit, contrast_matrix)
# final ebayes calc
ebayes = eBayes(fit2)
#all results
results = decideTests(ebayes, p.value = 0.05, adjust.method = 'BH', lfc = 1) #pvalue and LFC can be adjusted
summary(results)

#save significant testing overview
write.csv(summary(results), file = 'XXXX.csv')

# export final dataframe combining both imputed, non-impute and protein info data
export_df_genenames  = cbind(results, df_prot_valid, prot_valid_imp)
export_df_genenames <- left_join(export_df_genenames, prot_info, by = "Protein.Group")

#save export file
write.csv(export_df_genenames, file ="XXXX.csv")

