
# Name:    ASSIGN_run_predictions_single.R
#
# Purpose: Load a session created with ASSIGN_merge_and_combat.R and run a 
#          pathway prediction. This script runs one pathway at a time depending
#          on a number provided when running this script.
#
# Usage:   Rscript ASSIGN_run_predictions.R <pathway_number> <num_genes>
#
# Author:  David Jenkins (modified from ASSIGN scripts from Mumtahena Rahman)
# Date:    2015-09-30
#
################################################################################

if (length(commandArgs(trailingOnly=T)) < 1){
  print("ERROR: Submit the correct options!")
  print("Rscript ASSIGN_run_predictions.R <pathway_number> <num_genes>")
  quit(save = "no", status = 1, runLast = FALSE)
}

library(ASSIGN)

#----------------------------------------------------#
#Input Files (modify these locations for your system)#
#----------------------------------------------------#
working_dir         <- "/results"
basedir             <- "ASSIGN"
input_rda           <- "/results/ASSIGN/results.rda"
run_pathway         <- as.numeric(commandArgs(trailingOnly = T)[1]) 
num_genes           <- as.numeric(commandArgs(trailingOnly = T)[2])

#----------------------------------------------------#
#Parameters (modify these to change ASSIGN functions)#
#----------------------------------------------------#
sigma_sZero    <- 0.05
sigma_sNonZero <- 0.5
S_zeroPrior    <- FALSE

#---------#
#Load Data#
#---------#
setwd(working_dir)
load(input_rda)

#-------------------------------------------#
#running ASSIGN with helper function testSig#
#-------------------------------------------#
if(run_pathway == 1){
  trainingLabela <- list(control=list(akt=1:12),akt=13:18)
  sub_dir <- paste(basedir,paste("akt_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_akt),
                    test=c_test,
                    trainingLabel1=trainingLabela,
                    anchorGenes=list(akt=c("AKT1")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 2){
  trainingLabelb <- list(control=list(bad=1:12),bad=13:18)
  sub_dir <- paste(basedir,paste("bad_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_bad),
                    test=c_test,
                    trainingLabel1=trainingLabelb,
                    anchorGenes=list(bad=c("BAD")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 3){
  trainingLabel <- list(control=list(egfr=1:6),egfr=7:12)
  sub_dir <- paste(basedir,paste("egfr_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=train_egfr,
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(egfr=c("EGFR")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 4){
  trainingLabelh <- list(control=list(her2=1:12),her2=13:17)
  sub_dir <- paste(basedir,paste("her2_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_her2),
                    test=c_test,
                    trainingLabel1=trainingLabelh,
                    anchorGenes=list(her2=c("ERBB2")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 5){
  trainingLabeli <- list(control=list(igf1r=1:12),igf1r=13:18)
  sub_dir <- paste(basedir,paste("igf1r_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_igf1r),
                    test=c_test,
                    trainingLabel1=trainingLabeli,
                    anchorGenes=list(igf1r=c("IGF1R")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 6){
  trainingLabel <- list(control=list(krasgv=1:9),krasgv=10:18)
  sub_dir <- paste(basedir,paste("krasgv_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_kras_gfp,c_krasgv),
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(krasgv=c("KRAS")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 7){
  trainingLabel <- list(control=list(krasqh=1:9),krasqh=10:18)
  sub_dir <- paste(basedir,paste("krasqh_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_kras_gfp,c_krasqh),
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(krasqh=c("KRAS")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 8){
  trainingLabel <- list(control=list(kraswt=1:9),kraswt=10:18)
  sub_dir <- paste(basedir,paste("kraswt_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_kras_gfp,c_kraswt),
                    test=c_test,
                    trainingLabel1=trainingLabel,
                    anchorGenes=list(kraswt=c("KRAS")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}

if(run_pathway == 9){
  trainingLabelr <- list(control=list(raf=1:12),raf=13:18)
  sub_dir <- paste(basedir,paste("raf_",num_genes,"_gene_list", sep=""),sep='/')
  dir.create(sub_dir)
  assign_easy_multi_anchor_exclude(trainingData=cbind(c_gfp,c_raf),
                    test=c_test,
                    trainingLabel1=trainingLabelr,
                    anchorGenes=list(raf=c("RAF1")),
                    g=num_genes,
                    out_dir_base=sub_dir,
                    single=1,
                    sigma_sZero = sigma_sZero,
                    sigma_sNonZero = sigma_sNonZero,
                    S_zeroPrior=S_zeroPrior)
}
