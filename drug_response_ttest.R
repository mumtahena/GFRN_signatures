#Input files (Modify These to Your System's Locations)
icbp_predictions_file <- file.path("optimized_single_pathway_icbp.txt")
drug_response_file <- file.path("ICBP_drugs.txt")

#input icbp pathway data
single_pathway_best_icbp <-read.table(icbp_predictions_file, stringsAsFactors=FALSE,
                                      header=1, row.names=1,sep='\t')
rownames(single_pathway_best_icbp)[50] <- "T47D.Kbluc"

#input drug response data, subset to drugs of interest
drugs<-read.delim(drug_response_file, header=1, sep='\t',row.names=1)

rownames(drugs)[38] <- "SUM1315"
rownames(drugs)[45] <- "T47D.Kbluc"
rownames(drugs)[52] <- "MDAMB175"

drugs_all <- drugs[rownames(single_pathway_best_icbp),11:100]

drugs <- drugs[rownames(single_pathway_best_icbp),c("GSK2141795","PF.4691502","Everolimus","Sigma.AKT1.2.inhibitor","GSK1059868",
                                                    "Rapamycin","Temsirolimus","GSK2126458" ,"GSK2119563","GSK1059615","Lapatinib",
                                                    "AG1478","Gefitinib","AZD6244","Erlotinib","GSK1120212","Docetaxel","CGC.11047",
                                                    "Paclitaxel","Cisplatin")]
#Perform k-means clustering with 4 clusters
set.seed(123)
kmeans.4 <- kmeans(scale(single_pathway_best_icbp)[,c(1,2,3,4)], 4, iter.max=100000, nstart = 100)
kmeans.4.clusters <-data.frame(kmeans.4$cluster)
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "3",] <- "EGFR/BAD low"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "1",] <- "EGFR/BAD high"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "4",] <- "AKT/HER2 high"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "2",] <- "AKT/HER2 low"
colnames(kmeans.4.clusters) <- "KMeans.Clusters"

ttest_boxplots <- function(drugdata, clusters, selected_clustera, selected_clusterb, selected_drug){
  numa <- sum(!is.na(drugdata[rownames(clusters)[clusters == selected_clustera], selected_drug]))
  numb <- sum(!is.na(drugdata[rownames(clusters)[clusters == selected_clusterb], selected_drug]))
  if(numa < 2 | numb < 2){
    return(NA)
  }
  result <- t.test(drugdata[rownames(clusters)[clusters == selected_clustera], selected_drug],
                   drugdata[rownames(clusters)[clusters == selected_clusterb], selected_drug], na.rm=T)
  return(result$p.value)
}

akt.ttest.results <- unlist(lapply(colnames(drugs), function(x) ttest_boxplots(drugs, kmeans.4.clusters, "AKT/HER2 high", "AKT/HER2 low", x)))
akt.ttest.results.fdr <- p.adjust(akt.ttest.results, method="fdr")

egfr.ttest.results <- unlist(lapply(colnames(drugs), function(x) ttest_boxplots(drugs, kmeans.4.clusters, "EGFR/BAD high", "EGFR/BAD low", x)))
egfr.ttest.results.fdr <- p.adjust(egfr.ttest.results, method="fdr")

ttest_results <- data.frame(akt.ttest.results=akt.ttest.results, akt.ttest.results.fdr=akt.ttest.results.fdr, egfr.ttest.results=egfr.ttest.results, egfr.ttest.results.fdr=egfr.ttest.results.fdr, row.names=colnames(drugs))
