
# Name:    icbp_cor.R
#
# Purpose: Runs Spearman correlation between ASSIGN ICBP predictions and drug
#          response data.
#
# Usage:   Rscript icbp_cor.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-01-11
#
# Requirements: Key_ASSIGN_functions_balancedsig.R, ICBP_drugs.txt
#               Modify the locations of these scripts to their locations on your
#               system.
#
################################################################################

filenames<-system("ls */*/pathway_activity_testset*", intern=TRUE)

for(i in 1:length(filenames))
  {
   f<-read.csv(filenames[i], header=1,row.names=1) ###reading in the filess one at a time
   colnames(f)<-paste(filenames[i],colnames(f),sep='/')
   if(i!=1){
     print(i)
     data_icbp<-cbind(data_icbp,f)
    }
   else{
     data_icbp<-f
     }
  }

colnames(data_icbp)<-gsub(pattern = "/pathway_activity_testset.csv",replacement = "",x = colnames(data_icbp))
#head(data_icbp)
drugs<-read.delim("ICBP_drugs.txt", header=1, sep='\t',row.names=1)
key_assign_file     <- "Key_ASSIGN_functions_balancedsig.R"
source(key_assign_file)
icbp_drug<-merge_drop(data_icbp,drugs)
#colnames(icbp_drug)
cor_mat=p_mat=matrix(0,ncol(data_icbp),90)
rownames(cor_mat)=rownames(p_mat)=colnames(icbp_drug)[1:ncol(data_icbp)]
colnames(cor_mat)=colnames(p_mat)=colnames(icbp_drug)[(ncol(data_icbp)+11):ncol(icbp_drug)]

for(i in 1:ncol(data_icbp)){
  for(j in 1:90){
  temp=cor.test(icbp_drug[,i],icbp_drug[,(j+ncol(data_icbp)+10)],use="pairwise",method="spearman")
  print(j)
  print(temp)
  cor_mat[i,j]=temp$estimate
  p_mat[i,j]=temp$p.value
  }
}

cor_mat_df <- data.frame(cor_mat)
cor_mat_df <- cor_mat_df[,colnames(cor_mat_df) %in% c("Sigma.AKT1.2.inhibitor",
                         "Triciribine","BEZ235","Everolimus","GSK2126458",
                         "GSK2141795","Lapatinib","PF.4691502","Rapamycin",
                         "AG1478","AZD6244","BIBW2992","Erlotinib",
                         "GSK1120212","Gefitinib","PD98059","XRP44X",
                         "GSK1838705","Tykerb.IGF1R..1.1.",
                         "ERKi.II..FR180304.","Sorafenib")]
cor_mat_names <- t(sapply(strsplit(rownames(cor_mat_df),"/"),as.character))
colnames(cor_mat_names) <- c("run", "adapt","pathway")
cor_mat_df <- cbind(cor_mat_names,cor_mat_df)
rownames(cor_mat_df) <- NULL
write.table(cor_mat,"cor_drug_mat.txt",col.names = NA,quote=F,sep='\t')
write.table(cor_mat_df,"cor_drug_mat_subset.txt",row.names=F,col.names = T,quote=F,sep='\t')
colnames(p_mat)=paste(colnames(p_mat),"p_value",sep="_")
write.table(p_mat,"p_drug_mat.txt",col.names = NA,quote=F,sep='\t')
cor_p_mat<-cbind(cor_mat,p_mat)
order(colnames(cor_p_mat))
cor_p_mat<-cor_p_mat[,order(colnames(cor_p_mat))]
write.table(cor_p_mat,"cor_p_drug_mat.txt",col.names = NA,quote=F,sep='\t')
