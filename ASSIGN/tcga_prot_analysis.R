
# Name:    tcga_prot_analysis.R
#
# Purpose: Runs Spearman correlation between ASSIGN TCGA predictions and RPPA
#          protein data.
#
# Usage:   Rscript tcga_prot_analysis.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-01-11
#
# Requirements: Key_ASSIGN_functions_balancedsig.R, TCGA-BRCA-RBN.csv
#               Modify the locations of these scripts to their locations on your
#               system.
# 
################################################################################

source("/restricted/projectnb/pathsig/20150929_bild_paper_new_ASSIGN/scripts/Key_ASSIGN_functions.R")

prot_correlation<-function(predPath=NULL,protPath=NULL,nprot=190,offset=4){#offset=# of column in the protein file before protein data starts 
  data<-gatherFile(predPath)
  colnames(data)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data))
  prot<-read.csv(protPath,header=1)
  prot$TCGA_patient_barcode<-short_to_long_TCGA_id(shortnames=gsub(pattern="\\.","-",as.character(prot$TCGA_patient_barcode)),longnames=gsub(pattern="\\.","-",rownames(data)))
  rownames(data)<-gsub(pattern="\\.","-",rownames(data))
  data_prot<-merge(data,prot,by.x=0,by.y=1)
  cor_mat=p_mat=matrix(0,ncol(data),nprot)

  rownames(cor_mat)=rownames(p_mat)=colnames(data_prot)[2:(ncol(data)+1)]
  colnames(cor_mat)=colnames(p_mat)=colnames(data_prot)[(ncol(data)+1+offset):ncol(data_prot)]
  for(i in 2:(ncol(data)+1)){
    for(j in 1:190){
      print(j)
      temp=cor.test(data_prot[,i],data_prot[,(j+offset+ncol(data))],use="pairwise",method="spearman")
      cor_mat[(i-1),j]=temp$estimate
      p_mat[(i-1),j]=temp$p.value
    }
  }
  colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
  cor_p_mat<-cbind(cor_mat,p_mat)
  order(colnames(cor_p_mat))
  cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
  cur_dir<-getwd()
  setwd(predPath)
  write.table(cor_mat, "cor.txt",col.names=NA,sep='\t',quote=F)
  write.table(p_mat, "p_value.txt",col.names=NA,sep='\t',quote=F)
  write.table(cor_p_mat,"cor_p.txt",col.names=NA,sep='\t',quote=F)
  setwd(cur_dir)
  print((date()))
  cor_mat_df <- data.frame(cor_mat)
  cor_mat_df <- cor_mat_df[,colnames(cor_mat_df) %in% c("Akt","Akt_pS473",
                           "Akt_pT308","S6_pS235_S236","PDK1","PDK1_pS241","Bad_pS112",
                           "Bcl.xL","Bcl.2","XRP44X","HER2","HER2_pY1248","IRS1",
                           "IGFBP2","MAPK_pT202_Y204","C.Raf","C.Raf_pS338","A.Raf_pS299",
                           "B.Raf","MEK1_pS217_S221","EGFR","EGFR_pY1068","EGFR_pY1173",
                           "MEK1","STAT3_pY705","STAT5.alpha","AMPK_alpha","AMPK_pT172",
                           "PKC.alpha","PKC.alpha_pS657")]
  cor_mat_names <- t(sapply(strsplit(rownames(cor_mat_df),"/"),as.character))
  colnames(cor_mat_names) <- c("run", "adapt","pathway")
  cor_mat_df <- cbind(cor_mat_names,cor_mat_df)
  rownames(cor_mat_df) <- NULL
  write.table(cor_mat_df,"cor_prot_subset.txt",row.names=F,col.names = T,quote=F,sep='\t')
  return(cor_p_mat)
} 

datasets_dir <- "/data"
brca_prot    <- paste(datasets_dir,"TCGA-BRCA-RBN.csv",sep="/")

brca_one_step<- getwd()
brca_one_step_preds<-gatherFile(brca_one_step)
colnames(brca_one_step_preds)<-paste(colnames(brca_one_step_preds),"one_step",sep="_")
brca_one_step_cor_p<-prot_correlation(predPath=brca_one_step,protPath=brca_prot,nprot=190)
