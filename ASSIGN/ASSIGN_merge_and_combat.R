library(devtools)
library(sva)
library(data.table)

# Name:    ASSIGN_merge_and_combat.R
#
# Purpose: Merge together the signature data, and a test dataset and perform
#          the reference version of ComBat and save a session that can be used
#          with ASSIGN_run_predictions_single.R to run ASSIGN
#
# Usage:   Rscript ASSIGN_merge_and_combat.R
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-29
#
################################################################################

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
signatures_dir      <- "/signatures"
expr_file           <- paste(signatures_dir,"GSE83083_GFP18_AKT_BAD_HER2_IGF1R_RAF_GFP30_KRAS-G12V_tpmlog.txt",sep="/")
control_egfr_l_file <- paste(signatures_dir,"18_GFP_EGFR_TPMlog2.txt",sep="/")
key_assign_file     <- "/scripts/Key_ASSIGN_functions_balancedsig.R"
testFile            <- "/data/test_data/icbp_Rsubread_tpmlog.txt"

#--------------------------------------#
#Output Files (modify these every time)#
#--------------------------------------#
working_dir         <- "/results/ASSIGN"
output_rda          <- "results.rda"

#---------#
#Load Data#
#---------#
source(key_assign_file)
setwd(working_dir)
expr<-as.matrix(read.table(expr_file,sep='\t',row.names=1,header=1))
expr_batch_1<-subset(expr, select=GFP18.1:RAF.6)
expr_batch_1_f <-expr_batch_1[apply(expr_batch_1[,1:41]==0,1,mean) < 0.85,]
control_egfr_l<-read.table(control_egfr_l_file, sep='\t', header=1, row.names=1)
gfp_egfr_multi_f <- merge_drop(control_egfr_l,expr_batch_1_f)
gfp_kras<-subset(expr,select=GFP30.1:KRAS_G12V.9)
gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f,gfp_kras)
#load in test data frame
test<-data.frame(fread(testFile), check.names=F,row.names=1)

#------#
#ComBat#
#------#
sub<-c(6,6,12,6,6,5,6,6,9,9)
pdf("batch_twostep.pdf")
pcaplot(gfp_egfr_kras_multi_f,sub)
bat<-as.data.frame(cbind(c(rep(1,12),rep(2,41),rep(3,18)),c(rep(1,6),rep(2,6),rep(1,12),rep(3,6),rep(4,6),rep(5,5),rep(6,6),rep(7,6),rep(1,9),rep(8,9))))
colnames(bat)<-c("Batch","Model")
rownames(bat)<-colnames(gfp_egfr_kras_multi_f)
mod <- model.matrix(~as.factor(bat$Model))
combat_expr<-ComBat(dat=gfp_egfr_kras_multi_f, batch=bat[,1], mod=mod, ref.batch=2)
pcaplot(combat_expr,sub)
dat<-merge_drop(combat_expr,test)
sub<-c(6,6,12,6,6,5,6,6,9,9,ncol(test))
pcaplot(dat,sub)
bat<-as.matrix(cbind(colnames(dat),c(rep(1,ncol(gfp_egfr_kras_multi_f)),rep(2,ncol(test)))))
combat_expr1<-ComBat(dat=dat,batch=bat[,2], mod=NULL, ref.batch=1)
pcaplot(combat_expr1,sub)
dev.off()
c_gfp<-subset(combat_expr1, select=GFP18.1:GFP18.12)
c_akt<-subset(combat_expr1, select=AKT.1:AKT.6)
c_bad<-subset(combat_expr1, select=BAD.1:BAD.6)
c_her2<-subset(combat_expr1, select=HER2.1:HER2.6)
c_igf1r<-subset(combat_expr1, select=IGF1R.1:IGF1R.6)
c_raf<-subset(combat_expr1, select=RAF.1:RAF.6)
train_egfr<-combat_expr1[,1:12]
c_egfr_gfp <- train_egfr[,1:6]
c_egfr <- train_egfr[,7:12]
c_kras_gfp<-subset(combat_expr1,select=GFP30.1:GFP30.9)
c_krasgv<-subset(combat_expr1,select=KRAS_G12V.1:KRAS_G12V.9)
c_test<-combat_expr1[,(ncol(gfp_egfr_kras_multi_f)+1):ncol(combat_expr1)]

save.image(file=output_rda)
