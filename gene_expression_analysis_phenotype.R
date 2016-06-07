icbp<-as.matrix(read.table("~/Dropbox/bild_signatures/Datasets/icbp_Rsubread_tpmlog.txt", sep='\t', stringsAsFactors=FALSE, header=1, row.names=1))
exp<-t(icbp[c("BAX","BAK1","BID","BCL2L11","BAD","BIK","PMAIP1","HRK","BBC3","BMF","BCL2","BCL2L1","BCL2L2","MCL1","BCL2A1"),])
rownames(exp)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")
heatmap(exp)
source('~/Dropbox/bild_signatures/bild_signatures/Rmarkdowns_scripts//Key_ASSIGN_functions.Rmd', echo=TRUE)
single_pathway_best_icbp<-gatherFile("~/Desktop/20160201_ICBP_single_anchorGenes/")[,c(
  "akt_20_gene_list/adap_adap_single/pathway_activity_testset.csv/akt"
  ,"bad_250_gene_list/adap_adap_single/pathway_activity_testset.csv/bad"
  ,"egfr_25_gene_list/adap_adap_single/pathway_activity_testset.csv/egfr"
  ,"her2_10_gene_list/adap_adap_single/pathway_activity_testset.csv/her2"
  ,"igf1r_100_gene_list/adap_adap_single/pathway_activity_testset.csv/igf1r"
  ,"krasgv_175_gene_list/adap_adap_single/pathway_activity_testset.csv/krasgv"
  #,"krasqh_300_gene_list/adap_adap_single/pathway_activity_testset.csv/krasqh"
  ,"raf_350_gene_list/adap_adap_single/pathway_activity_testset.csv/raf"
)]
colnames(single_pathway_best_icbp) <-
  toupper(gsub(pattern = "_.*",replacement = "",colnames(single_pathway_best_icbp)))
colnames(single_pathway_best_icbp)

rownames(single_pathway_best_icbp)[1:7]<-c("184A1","184B5","21MT1","21MT2","21NT","21PT","600MPE")
single_pathway_best_icbp$phenotype<-NA
single_pathway_best_icbp[,1:8]<- scale(single_pathway_best_icbp[,1:8],scale = T,center = T)

akt<-apply(single_pathway_best_icbp[,c(1,4,5)],1,mean)
egfr<-apply(single_pathway_best_icbp[,c(2,3,6:8)],1,mean)
single_pathway_best_icbp[,1:7]<- scale(single_pathway_best_icbp[,1:7],scale = T,center = T)

akt<-apply(single_pathway_best_icbp[,c(1,4,5)],1,mean)
egfr<-apply(single_pathway_best_icbp[,c(2,3,6:8)],1,mean)

for(i in 1:nrow(single_pathway_best_icbp))
{ print(i)
  if(akt[[i]] > egfr[[i]]){
    single_pathway_best_icbp[i,9] <-"AKT Phenotype"
  }
  else{
    single_pathway_best_icbp[i,9] <- "EGFR Phenotype"
  }
}

#write.table(cbind(akt,egfr,single_pathway_best_icbp),"~/Dropbox/akt_egfr_phenotype_0220.txt",sep='\t',quote=F, col.names = NA)
exp_best<-merge_drop(single_pathway_best_icbp,exp)
colnames(exp_best)
#pdf("~/Dropbox/bild_signatures/Validations/ICBP/phenotype_gene_expression_0221.pdf")

par(mfrow = c(1,1),lwd=4)
for (i in 10:ncol(exp_best)) {
  tmp = t.test(exp_best[,i] ~ exp_best[,9])
  boxplot(exp_best[,i] ~ exp_best[,9],main = paste(colnames(exp_best)[i],"\np-value =",tmp$p.value))
  
}
#dev.off()

#######Drug response assay based analysis#######
assay_ec50<- read.table("~/Dropbox/Bild drug screen 2015/Plate_Layouts/Supplementary/BR_OV_EC50.txt",sep='\t',header=1,row.names=1)
assay_ec50_log10<- -apply((assay_ec50+0.0000000001),2,log10)##converting the EC50 values to -log10(EC50 in microM)
assay_ec50_preds<-merge_drop(assay_ec50_log10,single_pathway_best)
colnames(assay_ec50_preds)<-gsub("EC50_","",colnames(assay_ec50_preds))
colnames(assay_ec50_preds)
assay_ec50_preds$MK2206<-NULL
assay_ec50_preds$Paclitaxel<-NULL
assay_ec50_preds$NAV.4471
assay_ec50_preds$NAV.4471<-NULL
assay_ec50_preds$RA190<-NULL
assay_ec50_preds$Sorafinib<-NULL
assay_ec50_preds$SAHA<-NULL
assay_ec50_preds$Carboplatin<-NULL
colnames(assay_ec50_preds)
colnames(assay_ec50_preds)[7]="Navitaclax"
colnames(assay_ec50_preds)[10]="Sigma AKTi"
colnames(assay_ec50_preds)[9]="Palbociclib"
colnames(assay_ec50_preds)
#heatmap.2(as.matrix(cor(assay_ec50_preds[,1:14],assay_ec50_preds[,15:22],use="pairwise")),margins =c(10,10), col=bluered,dendrogram="none", trace="none",main=paste("Pathway-drug response","in breast cell lines",sep='\n '))#,cellnote = round(cors,digits = 2),notecol = 'black',density.info = 'none')
#heatmap.2(as.matrix(cor(assay_ec50_preds[,1:14],assay_ec50_preds[,15:22],use="pairwise")),margins =c(10,10), col=bluered, trace="none",main="",scale="row")#,cellnote = round(cors,digits = 2),notecol = 'black',density.info = 'none')
heatmap.2(as.matrix(cor(assay_ec50_preds[,1:11],assay_ec50_preds[,12:19],use="pairwise")),margins =c(10,10), col=bluered, trace="none",main="",scale="row")#,cellnote = round(cors,digits = 2),notecol = 'black',density.info = 'none')



library(mclust)
library(reshape2)
library(ggplot2)


par(mfrow=c(1,1),lwd=4)
#pdf("~/Dropbox/drug_assay_phenotype_boxplots.pdf")
for(i in 1:11){
    tmp<-t.test(assay_ec50_preds[,i]~assay_ec50_preds$phenotype)
    boxplot2(assay_ec50_preds[,i]~assay_ec50_preds$phenotype,main=paste(colnames(assay_ec50_preds)[i],"\np-value:",round(tmp$p.value,digits = 5),sep=' '),ylab="Sensitivity")#,notch=T,col=c("red","green"),horizontal = T,type='l',add=T)       
      
}
#dev.off()

########ICBP Phenotype based Protein analysis 
prot<-read.table("~/Dropbox/bild_signatures/bild_signatures/Datasets/proteomics.txt",sep='\t',header=1,row.names=1)

rownames(prot)[1:3]<-c("184A1","184B5","600MPE")
pred_prot<-merge_drop(single_pathway_best_icbp,prot)
dim(prot)
dim(pred_prot)
#pdf("~/Dropbox/bild_signatures/Validations/ICBP/icbp_protein_phenotype_boxplots_0221.pdf")
par(lwd=4)
  tmp<-t.test(pred_prot$Bcl2~pred_prot$phenotype)
  boxplot2(pred_prot$Bcl2~pred_prot$phenotype,main=paste("BCL2","\np-value:",round(tmp$p.value,digits = 5),sep=' '),ylab="RPPA Score")#,notch=T,col=c("red","green"),horizontal = T,type='l',add=T)       
  boxplot(pred_prot$Bcl2~pred_prot$phenotype,main=paste("BCL2","\np-value:",round(tmp$p.value,digits = 5),sep=' '),ylab="RPPA Score")#,notch=T,col=c("red","green"),horizontal = T,type='l',add=T)       
  

#dev.off()


########TCGA Analysis#######
single_pathway_best_tcga<- gatherFile("~/Desktop/20160201_TCGA_single_anchorGenes/")[,c(
  "akt_20_gene_list/adap_adap_single/pathway_activity_testset.csv/akt"
  ,"bad_250_gene_list/adap_adap_single/pathway_activity_testset.csv/bad"
  ,"egfr_25_gene_list/adap_adap_single/pathway_activity_testset.csv/egfr"
  ,"her2_10_gene_list/adap_adap_single/pathway_activity_testset.csv/her2"
  ,"igf1r_100_gene_list/adap_adap_single/pathway_activity_testset.csv/igf1r"
  ,"krasgv_175_gene_list/adap_adap_single/pathway_activity_testset.csv/krasgv"
  ,"krasqh_300_gene_list/adap_adap_single/pathway_activity_testset.csv/krasqh"
  ,"raf_350_gene_list/adap_adap_single/pathway_activity_testset.csv/raf"
)]

colnames(single_pathway_best_tcga) <-
  toupper(gsub(pattern = "_.*",replacement = "",colnames(single_pathway_best_tcga)))
colnames(single_pathway_best_tcga)


single_pathway_best_tcga$phenotype<-NA
single_pathway_best_tcga[,1:8]<- scale(single_pathway_best_tcga[,1:8],scale = T,center = T)

akt<-apply(single_pathway_best_tcga[,c(1,4,5)],1,mean)
egfr<-apply(single_pathway_best_tcga[,c(2,3,6:8)],1,mean)
for(i in 1:nrow(single_pathway_best_tcga))
{ print(i)
  if(akt[[i]] > egfr[[i]]){
    single_pathway_best_tcga[i,9] <-"AKT Phenotype"
  }
  else{
    single_pathway_best_tcga[i,9] <- "EGFR Phenotype"
  }
}

prot_brca<-read.csv("~/Dropbox/bild_signatures/Datasets/TCGA-BRCA-RBN.csv",header=1,check.names = F)
prot_brca$TCGA_patient_barcode<-short_to_long_TCGA_id(shortnames=gsub(pattern="\\.","-",as.character(prot_brca$TCGA_patient_barcode)),longnames=gsub(pattern="\\.","-",rownames(data_tcga)))


data_prot<-merge(single_pathway_best_tcga,prot_brca,by.x=0,by.y=1)
prot_of_interest<-c("BAX","BAK1","BID","BIM","BAD","BIK","NOXA","HRK","PUMA","BMF","BCL2","BCLXL","BCLW","MCL1","BFL1")
colnames(data_prot)<-toupper(colnames(data_prot))
colnames(data_prot)<-gsub(pattern = "-",replacement = "",colnames(data_prot))
#pdf("~/Dropbox/bild_signatures/Validations/TCGA/TCGA_phenotype_protein_rppa_0312.pdf")
for(i in 1:length(prot_of_interest)){
  print(paste("analyzing",prot_of_interest[i]))
  if(prot_of_interest[i]%in%colnames(data_prot)){
    print(paste(prot_of_interest[i],"data available"))
    tmp<-t.test(data_prot[,prot_of_interest[i]]~data_prot$PHENOTYPE)
    boxplot(data_prot[,prot_of_interest[i]]~data_prot$PHENOTYPE,main=paste(prot_of_interest[i],"\np-value:",tmp$p.value,sep=' '),ylab="RPPA Score")#,notch=T,col=c("red","green"),horizontal = T,type='l',add=T)       
    
  }
}
tmp<-t.test(data_prot$BAD_PS112~data_prot$PHENOTYPE)
boxplot(data_prot$BAD_PS112~data_prot$PHENOTYPE,main=paste("BAD_ps112","\np-value:",tmp$p.value,sep=' '),ylab="RPPA Score")

#dev.off()

###TCGA gene expression based phenotype analysis
test<-data.frame(fread("~/Dropbox/bild_signatures/Datasets/PANCAN24_BRCA_1119_TPMlog2.txt"), check.names=F,row.names=1)

exp_tcga<-t(test[c("BAX","BAK1","BID","BCL2L11","BAD","BIK","PMAIP1","HRK","BBC3","BMF","BCL2","BCL2L1","BCL2L2","MCL1","BCL2A1"),])

exp_best<-merge_drop(single_pathway_best_tcga,exp_tcga)
colnames(exp_best)
#pdf("~/Dropbox/bild_signatures/Validations/TCGA/phenotype_gene_expression_0221.pdf")

par(mfrow = c(1,1), lwd=4)
for (i in 10:ncol(exp_best)) {
  tmp = t.test(exp_best[,i] ~ exp_best[,9])

  boxplot2(exp_best[,i] ~ exp_best[,9],main = paste(colnames(exp_best)[i],"\np-value =",tmp$p.value))
  boxplot(exp_best[,i] ~ exp_best[,9],main = paste(colnames(exp_best)[i],"\np-value =",tmp$p.value))
  
  
}
#dev.off()

########Heatmap of selected cell lines used in Western Blot protein analysis
selected_celllines<-c("AU565","BT474","BT483","BT549","CAMA1","HCC1143","HCC1419","HCC1428","HCC1806","HCC1937","HCC1954","HCC38","HCC70","HS578T","JIMT1","MCF7","MDAMB175","MDAMB231","T47D","ZR751")
scaled_selected<-single_pathway_best_icbp[selected_celllines,]
#pdf("~/Dropbox/bild_signatures/Validations/ICBP/heatmap_pathways_20_celllines_0221.pdf")
par(lwd=4)
heatmap.2(as.matrix(scaled_selected[,1:8]),margins =c(10,14), col=my_palette, trace="none",main="",dendrogram = "column",ColSideColors = c("coral3","aquamarine4","aquamarine4","coral3","coral3","aquamarine4","aquamarine4","aquamarine4"),scale = "row")#,cellnote = round(cors,digits = 2),notecol = 'black',density.info = 'none')
par(lend = 1)           # square line ends for the color legend
legend("topright",legend = c("HER2/IGF1R/AKT phenotype", "BAD/EGFR/KRAS/RAF phenotype"), col = c("coral3","aquamarine4"),  lty= 1,lwd = 10,cex = 0.5)
#dev.off()