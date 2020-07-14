#### aPrioriGenes-Corr HPCC ###

## Load packages and data
library(tidyverse)
library(edgeR)

setwd("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript")

## Load data
annot <- read.delim("./extData/Trinotate.csv", header = TRUE, sep = "\t")
cnt.trans <- read.table("./extData/rsem.isoform.counts.matrix", header = TRUE)
cnt.gene <- read.table("./extData/rsem.gene.counts.matrix", header = TRUE)
covars <- read.csv("./extData/BWTE54_SSgroup_Raw_Pool.csv", header = TRUE)
targets <- read_delim("./extData/aPrioriTranscripts_V2.csv", delim = ",") %>%
  filter(include == "Yes" | include == "yes")

## Clean data
cnt.bursa.gene <- cnt.gene %>%
  select(-alignEstimateAbundance_BWTE_Ileum_36_S50,
         -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
         -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("gene") %>%
  separate(gene, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  column_to_rownames("gene") %>%
  select(contains("Bursa"))

cnt.ileum.gene <- cnt.gene %>%
  select(-alignEstimateAbundance_BWTE_Ileum_36_S50,
         -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
         -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("gene") %>%
  separate(gene, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  column_to_rownames("gene") %>%
  select(contains("Ileum"))

cnt.bursa.trans <- cnt.trans %>%
  select(-alignEstimateAbundance_BWTE_Ileum_36_S50,
         -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
         -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("transcript") %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", "isoform")) %>%
  unite(transcript, pt1, pt2, gene, isoform) %>%
  column_to_rownames("transcript") %>%
  select(contains("Bursa"))

cnt.ileum.trans <- cnt.trans %>%
  select(-alignEstimateAbundance_BWTE_Ileum_36_S50,
         -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
         -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("transcript") %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", "isoform")) %>%
  unite(transcript, pt1, pt2, gene, isoform) %>%
  column_to_rownames("transcript") %>%
  select(contains("Ileum"))

covars <- covars %>%
  filter(!bird %in% c("36", "19")) %>%
  arrange(bird) %>%
  mutate(group = str_remove(group, "-"))

annot <- annot %>%
  separate(transcript_id, into = c(NA, "pt1", "pt2", "gene", "isoform")) %>%
  unite(gene_id, pt1, pt2, gene, remove = FALSE) %>%
  unite(transcript_id, pt1, pt2, gene, isoform)

#### Calculate log(CPM) and assemble master DFs ####
#Convert to DGEList object
dge.bursa.trans <- DGEList(counts=cnt.bursa.trans)
dge.bursa.gene <- DGEList(counts=cnt.bursa.gene)
dge.ileum.trans <- DGEList(counts=cnt.ileum.trans)
dge.ileum.gene <- DGEList(counts=cnt.ileum.gene)

#CPM and log-CPM
lcpm.bursa.trans <- cpm(dge.bursa.trans, log = TRUE)
lcpm.bursa.gene <- cpm(dge.bursa.gene, log = TRUE)
lcpm.ileum.trans <- cpm(dge.ileum.trans, log = TRUE)
lcpm.ileum.gene <- cpm(dge.ileum.gene, log = TRUE)

## Master lcpm tibs
# Trans
lcpm.bursa.tmp <- lcpm.bursa.trans %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("transcript")
lcpm.ileum.tmp <- lcpm.ileum.trans %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("transcript")

lcpm.trans <- lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0") %>%
  filter(transcript %in% targets$transcript_id) %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  mutate(levelGT = "transcript", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, virus.sac, log.virus.sac, SSgroup.virus.sac) %>%
  summarize(meanTiter1to5 = log10(mean(titer, na.rm = TRUE)+1))


# Gene
lcpm.bursa.tmp <- lcpm.bursa.gene %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene")
lcpm.ileum.tmp <- lcpm.ileum.gene %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene")

lcpm.all <- lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0")  %>%
  filter(gene %in% targets$gene_id) %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  mutate(levelGT = "gene", identifier = gene) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, virus.sac, log.virus.sac, SSgroup.virus.sac) %>%
  summarize(meanTiter1to5 = log10(mean(titer, na.rm = TRUE)+1)) %>%
  bind_rows(lcpm.trans)

#### Gene query database ####
annot.all <- annot %>%
  select(transcript_id,
         gene_id,
         sprot_Top_BLASTX_hit,
         sprot_Top_BLASTP_hit,
         gene_ontology_BLASTX,
         gene_ontology_BLASTP,
         Kegg,
         eggnog,
         Pfam) %>%
  separate(sprot_Top_BLASTX_hit, into = c("sprot_geneName_BlastX", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
  separate(sprot2, sep = "=", into = c(NA, "sprot_geneFunction_BlastX")) %>%
  separate(sprot_geneFunction_BlastX, sep = ";", into = c("sprot_geneFunction_BlastX", NA)) %>%
  separate(sprot_Top_BLASTP_hit, into = c("sprot_geneName_BlastP", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
  separate(sprot2, sep = "=", into = c(NA, "sprot_geneFunction_BlastP")) %>%
  separate(sprot_geneFunction_BlastP, sep = ";", into = c("sprot_geneFunction_BlastP", NA)) %>%
  separate(gene_ontology_BLASTX, sep = "\\`", into = paste("GO_BlastX", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
  separate(gene_ontology_BLASTP, sep = "\\`", into = paste("GO_BlastP", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
  separate(Pfam, sep = "\\`", into = paste("Pfam", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
  as_tibble()


#### Analysis function ####
aprioriAnalysis.mean <- function(target, targetTissue, targetLevel, ...) {
  library(tidyverse)
  datSet <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue)

  correl.mean <- cor.test(datSet$meanTiter1to5, datSet$lcpm)

  as.data.frame(cbind(correl.mean$estimate, correl.mean$p.value)) %>% as_tibble() %>% rename("Est" = 1, "pval" = 2)
}

aprioriAnalysis.sac <- function(target, targetTissue, targetLevel, ...) {
  library(tidyverse)
  datSet <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue)

  correl.sac <- cor.test(datSet$log.virus.sac, datSet$lcpm)

  as.data.frame(cbind(correl.sac$estimate, correl.sac$p.value)) %>% as_tibble() %>% rename("Est" = 1, "pval" = 2)
}



#### Run analysis loop ####
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

targetTissue <- "Ileum" ## Bursa or Ileum
targetLevel <- "gene" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

finalMatrix.mean.IG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IG = aprioriAnalysis.mean(z, targetTissue, targetLevel)
  tmpMatrix.IG
}

finalMatrix.sac.IG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IG = aprioriAnalysis.sac(z, targetTissue, targetLevel)
  tmpMatrix.IG
}


#### Run analysis loop ####
targetTissue <- "Ileum" ## Bursa or Ileum
targetLevel <- "transcript" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)
results.mean <- list()
results.sac <- list()


finalMatrix.mean.IT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IT = aprioriAnalysis.mean(z, targetTissue, targetLevel)
  tmpMatrix.IT
}

finalMatrix.sac.IT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IT = aprioriAnalysis.sac(z, targetTissue, targetLevel)
  tmpMatrix.IT
}


#### Run analysis loop ####
targetTissue <- "Bursa" ## Bursa or Ileum
targetLevel <- "transcript" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)
results.mean <- list()
results.sac <- list()


finalMatrix.mean.BT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BT = aprioriAnalysis.mean(z, targetTissue, targetLevel)
  tmpMatrix.BT
}

finalMatrix.sac.BT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BT = aprioriAnalysis.sac(z, targetTissue, targetLevel)
  tmpMatrix.BT
}



#### Run analysis loop ####
targetTissue <- "Bursa" ## Bursa or Ileum
targetLevel <- "gene" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)
results.mean <- list()
results.sac <- list()

finalMatrix.mean.BG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BG = aprioriAnalysis.mean(z, targetTissue, targetLevel)
  tmpMatrix.BG
}

finalMatrix.sac.BG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BG = aprioriAnalysis.sac(z, targetTissue, targetLevel)
  tmpMatrix.BG
}


#stop cluster
stopCluster(cl)






## Recreate the "sets" for further analysis
set.BT <- lcpm.all %>%
  filter(levelGT == "transcript", tissue == "Bursa") %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

set.BG <- lcpm.all %>%
  filter(levelGT == "gene", tissue == "Bursa") %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

set.IT <- lcpm.all %>%
  filter(levelGT == "transcript", tissue == "Ileum") %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

set.IG <- lcpm.all %>%
  filter(levelGT == "gene", tissue == "Ileum") %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)


## Process DPI1-5 mean results
finalMatrix.mean.BG.sig <- finalMatrix.mean.BG %>%
  mutate(identifier = unique(set.BG$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "mean.BG")

finalMatrix.mean.BT.sig <- finalMatrix.mean.BT %>%
  mutate(identifier = unique(set.BT$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "mean.BT")

finalMatrix.mean.IG.sig <- finalMatrix.mean.IG %>%
  mutate(identifier = unique(set.IG$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "mean.IG")

finalMatrix.mean.IT.sig <- finalMatrix.mean.IT %>%
  mutate(identifier = unique(set.IT$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "mean.IT")

## Process sacrifice day results
finalMatrix.sac.BG.sig <- finalMatrix.sac.BG %>%
  mutate(identifier = unique(set.BG$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "sac.BG")

finalMatrix.sac.BT.sig <- finalMatrix.sac.BT %>%
  mutate(identifier = unique(set.BT$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "sac.BT")

finalMatrix.sac.IG.sig <- finalMatrix.sac.IG %>%
  mutate(identifier = unique(set.IG$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "sac.IG")

finalMatrix.sac.IT.sig <- finalMatrix.sac.IT %>%
  mutate(identifier = unique(set.IT$identifier)) %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1) %>%
  mutate(comparison = "sac.IT")


#### Bind everything up ####
sigResults.mean <- bind_rows(finalMatrix.mean.BG.sig,
                             finalMatrix.mean.BT.sig,
                             finalMatrix.mean.IG.sig,
                             finalMatrix.mean.IT.sig)

sigResults.sac <- bind_rows(finalMatrix.sac.BG.sig,
                            finalMatrix.sac.BT.sig,
                            finalMatrix.sac.IG.sig,
                            finalMatrix.sac.IT.sig)

### Clean up and save ###
remove(
  annot,
  targets,
  cnt.bursa.gene,
  cnt.bursa.trans,
  cnt.gene,
  cnt.ileum.gene,
  cnt.ileum.trans,
  covars,
  dge.bursa.gene,
  dge.bursa.trans,
  dge.ileum.gene,
  dge.ileum.trans,
  lcpm.bursa.gene,
  lcpm.bursa.tmp,
  lcpm.bursa.trans,
  lcpm.ileum.gene,
  lcpm.ileum.tmp,
  lcpm.ileum.trans,
  lcpm.trans,
  annot.all,
  cl,
  cnt.trans,
  finalMatrix,
  finalMatrix.IG,
  lcpm.all,
  results.mean,
  results.sac,
  set,
  set.BG,
  set.BT,
  set.IG,
  set.IT,
  tmpMatrix
)

save.image("aPrioriGenes-Corr-log.Rws")
