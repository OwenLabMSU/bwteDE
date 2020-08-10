## Load packages and data
library(tidyverse)
library(edgeR)
library(broom)
library(kableExtra)
library(ggrepel)

setwd("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/")

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
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, group, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, log.virus.sac, group) %>%
  summarize(log.meanTiter1to5 = log10(mean(titer, na.rm = TRUE)))

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
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, group, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, log.virus.sac, group) %>%
  summarize(log.meanTiter1to5 = log10(mean(titer + 1, na.rm = TRUE))) %>%
  bind_rows(lcpm.trans) %>%
  mutate(group = recode(group,
                        C1 = "Ctl",
                        C14 = "Ctl")) %>%
  mutate(group = factor(group,
                        levels = c("Ctl", "I1", "I3", "I5", "I14")))

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


##### Analysis function ######
candidateGeneAnalysis.Reg <- function(target, targetTissue, targetLevel, ...) {
  library(tidyverse)
  filtered.dat <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue)

  sum.overall <- filtered.dat %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    summary()

  res.overall <- c(sum.overall$adj.r.squared,
                   sum.overall$coefficients[2,1],
                   sum.overall$coefficients[2,4]) %>%
    as_tibble() %>%
    mutate(result.type = c("adj.r.squared", "estimate", "pval"),
           identifier = target,
           results.value = value,
           group = "all") %>%
    select(-value)

  sum.I1 <- filtered.dat %>%
    filter(group == "I1") %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    summary()

  res.I1 <- c(sum.I1$adj.r.squared,
              sum.I1$coefficients[2,1],
              sum.I1$coefficients[2,4]) %>%
    as_tibble() %>%
    mutate(result.type = c("adj.r.squared", "estimate", "pval"),
           identifier = target,
           results.value = value,
           group = "I1") %>%
    select(-value)

  sum.I3 <- filtered.dat %>%
    filter(group == "I3") %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    summary()

  res.I3 <- c(sum.I3$adj.r.squared,
              sum.I3$coefficients[2,1],
              sum.I3$coefficients[2,4]) %>%
    as_tibble() %>%
    mutate(result.type = c("adj.r.squared", "estimate", "pval"),
           identifier = target,
           results.value = value,
           group = "I3") %>%
    select(-value)

  sum.I5 <- filtered.dat %>%
    filter(group == "I5") %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    summary()

  res.I5 <- c(sum.I5$adj.r.squared,
              sum.I5$coefficients[2,1],
              sum.I5$coefficients[2,4]) %>%
    as_tibble() %>%
    mutate(result.type = c("adj.r.squared", "estimate", "pval"),
           identifier = target,
           results.value = value,
           group = "I5") %>%
    select(-value)

  sum.I14 <- filtered.dat %>%
    filter(group == "I14") %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    summary()

  res.I14 <- c(sum.I14$adj.r.squared,
               sum.I14$coefficients[2,1],
               sum.I14$coefficients[2,4]) %>%
    as_tibble() %>%
    mutate(result.type = c("adj.r.squared", "estimate", "pval"),
           identifier = target,
           results.value = value,
           group = "I14") %>%
    select(-value)

  bind_rows(res.overall, res.I1, res.I3, res.I5, res.I14)
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

finalMatrix.IG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IG = candidateGeneAnalysis.Reg(z, targetTissue, targetLevel)
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


finalMatrix.IT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.IT = candidateGeneAnalysis.Reg(z, targetTissue, targetLevel)
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


finalMatrix.BT <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BT = candidateGeneAnalysis.Reg(z, targetTissue, targetLevel)
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

finalMatrix.BG <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix.BG = candidateGeneAnalysis.Reg(z, targetTissue, targetLevel)
  tmpMatrix.BG
}


#stop cluster
stopCluster(cl)



## Identify significant patterns using the overall dataset
finalMatrix.BG.sig <- finalMatrix.BG %>%
  filter(group == "all", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "BG")

finalMatrix.BT.sig <- finalMatrix.BT %>%
  filter(group == "all", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "BT")

finalMatrix.IG.sig <- finalMatrix.IG %>%
  filter(group == "all", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "IG")

finalMatrix.IT.sig <- finalMatrix.IT %>%
  filter(group == "all", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "IT")


## Identify significant patterns using just I1
finalMatrix.BG.I1.sig <- finalMatrix.BG %>%
  filter(group == "I1", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "BG")

finalMatrix.BT.I1.sig <- finalMatrix.BT %>%
  filter(group == "I1", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "BT")

finalMatrix.IG.I1.sig <- finalMatrix.IG %>%
  filter(group == "I1", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "IG")

finalMatrix.IT.I1.sig <- finalMatrix.IT %>%
  filter(group == "I1", result.type == "pval") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "IT")


#### Bind everything up ####
sigResults <- bind_rows(finalMatrix.BG.sig,
                        finalMatrix.BT.sig,
                        finalMatrix.IG.sig,
                        finalMatrix.IT.sig)



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
  lcpm.all,
  results.mean,
  results.sac,
  set,
  tmpMatrix
)

save.image("candidateGeneAnalysis_regression.Rws")
