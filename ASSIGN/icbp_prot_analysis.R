
# Name:    icbp_prot_analysis.R
#
# Purpose: Runs Spearman correlation between ASSIGN ICBP predictions and RPPA
#          protein data.
#
# Usage:   Rscript icbp_prot_analysis.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-01-11
#
# Requirements: Key_ASSIGN_functions_balancedsig.R, proteomics.txt
#               Modify the locations of these scripts to their locations on your
#               system.
#
################################################################################

source("/restricted/projectnb/pathsig/20150929_bild_paper_new_ASSIGN/scripts/Key_ASSIGN_functions.R")

prot_correlation<-function(predPath=NULL,protPath=NULL,nprot=190,offset=4){#offset=# of column in the protein file before protein data starts 
  data<-gatherFile(predPath)
  for(i in grep("^[0-9]",rownames(data))){
    rownames(data)[i] <- paste("X",rownames(data)[i], sep="")
  }
  colnames(data)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data))
  prot<-read.table(protPath,header=1)
  data_prot<-merge(data,prot,by.x=0,by.y=0)
  cor_mat=p_mat=matrix(0,ncol(data),nprot)

  rownames(cor_mat)=rownames(p_mat)=colnames(data_prot)[2:(ncol(data)+1)]
  colnames(cor_mat)=colnames(p_mat)=colnames(data_prot)[(ncol(data)+1+offset):ncol(data_prot)]
  for(i in 2:(ncol(data)+1)){
    for(j in 1:nprot){
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
  write.table(cor_mat, "cor_prot.txt",col.names=NA,sep='\t',quote=F)
  write.table(p_mat, "p_value_prot.txt",col.names=NA,sep='\t',quote=F)
  write.table(cor_p_mat,"cor_p_prot.txt",col.names=NA,sep='\t',quote=F)
  setwd(cur_dir)
  print((date()))
  cor_mat_df <- data.frame(cor_mat)
  cor_mat_df <- cor_mat_df[,colnames(cor_mat_df) %in% c("Akt","Aktp308",
                           "Aktp473","Bcl2","EGFR","EGFRp1068","HER2","HER2p1248",
                           "PDK1","S6p235.236","S6p240.244","IGFR1","PDK1p241","MAPKp", "AMPKp", "MEK1","PKCalpha","PKCalphap657")]
  cor_mat_names <- t(sapply(strsplit(rownames(cor_mat_df),"/"),as.character))
  colnames(cor_mat_names) <- c("run", "adapt","pathway")
  cor_mat_df <- cbind(cor_mat_names,cor_mat_df)
  rownames(cor_mat_df) <- NULL
  write.table(cor_mat_df,"cor_prot_subset.txt",row.names=F,col.names = T,quote=F,sep='\t')
  return(cor_p_mat)
} 

icbp_prot    <- "proteomics.txt"

brca_one_step<- getwd()
brca_one_step_preds<-gatherFile(brca_one_step)
brca_one_step_cor_p<-prot_correlation(predPath=brca_one_step,protPath=icbp_prot,nprot=70,offset=1)
