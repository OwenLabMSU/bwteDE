#### aPrioriGenes-Corr HPCC ###

## Load packages and data
library(tidyverse)
library(edgeR)

setwd("/mnt/research/avian/teal")

## Load data
annot <- read.delim("./SRC/Trinotate.csv", header = TRUE, sep = "\t")
cnt.trans <- read.table("./SRC/rsem.isoform.counts.matrix", header = TRUE)
cnt.gene <- read.table("./SRC/rsem.gene.counts.matrix", header = TRUE)
covars <- read.csv("./SRC/BWTE54_SSgroup_Raw_Pool.csv", header = TRUE)
targets <- read_delim("./SRC/aPrioriTranscripts_V2.csv", delim = ",") %>%
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
lcpm.bursa.trans <- cpm(dge.bursa.trans)
lcpm.bursa.gene <- cpm(dge.bursa.gene)
lcpm.ileum.trans <- cpm(dge.ileum.trans)
lcpm.ileum.gene <- cpm(dge.ileum.gene)

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
  #filter(group == "I1") %>%
  #filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "transcript", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, virus.sac, SSgroup.virus.sac) %>%
  summarize(meanTiter1to5 = mean(titer, na.rm = TRUE))


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
  #filter(group == "I1") %>%
  #filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "gene", identifier = gene) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac, virus.dpi1, virus.dpi2, virus.dpi3, virus.dpi4, virus.dpi5) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  pivot_longer(cols = virus.dpi1:virus.dpi5,
               names_to = "dpi",
               values_to = "titer") %>%
  group_by(identifier, bird, levelGT, tissue, lcpm, virus.sac, SSgroup.virus.sac) %>%
  summarize(meanTiter1to5 = mean(titer, na.rm = TRUE)) %>%
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
  datSet <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue)

  correl.mean <- cor.test(datSet$meanTiter1to5, datSet$lcpm)

  as.data.frame(cbind(correl.mean$estimate, correl.mean$p.value)) %>% as_tibble() %>% rename("Est" = 1, "pval" = 2)
}

aprioriAnalysis.sac <- function(target, targetTissue, targetLevel, ...) {
  datSet <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue)

  correl.sac <- cor.test(datSet$virus.sac, datSet$lcpm)

  as.data.frame(cbind(correl.sac$estimate, correl.sac$p.value)) %>% as_tibble() %>% rename("Est" = 1, "pval" = 2)
}



#### Run analysis loop ####
targetTissue <- "Ileum" ## Bursa or Ileum
targetLevel <- "gene" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

results.mean <- list()
results.sac <- list()


for(z in unique(set$identifier)) {
  results.mean[[length(results.mean)+1]] <- aprioriAnalysis.mean(z, targetTissue, targetLevel)
  results.sac[[length(results.mean)+1]] <- aprioriAnalysis.sac(z, targetTissue, targetLevel)
}

## Process DPI1-5 mean results
results.IG.mean.tib <- bind_rows(results.mean, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.IG.mean.tmp <- results.IG.mean.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.IG.mean.sig <- results.IG.mean.tib %>% filter(identifier %in% results.IG.mean.tmp$identifier)

## Process day of sacrifice results
results.IG.sac.tib <- bind_rows(results.sac, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.IG.sac.tmp <- results.IG.sac.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.IG.sac.sig <- results.IG.sac.tib %>% filter(identifier %in% results.IG.sac.tmp$identifier)


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


for(z in unique(set$identifier)) {
  results.mean[[length(results.mean)+1]] <- aprioriAnalysis.mean(z, targetTissue, targetLevel)
  results.sac[[length(results.mean)+1]] <- aprioriAnalysis.sac(z, targetTissue, targetLevel)
}

## Process DPI1-5 mean results
results.IT.mean.tib <- bind_rows(results.mean, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.IT.mean.tmp <- results.IT.mean.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.IT.mean.sig <- results.IT.mean.tib %>% filter(identifier %in% results.IT.mean.tmp$identifier)

## Process day of sacrifice results
results.IT.sac.tib <- bind_rows(results.sac, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.IT.sac.tmp <- results.IT.sac.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.IT.sac.sig <- results.IT.sac.tib %>% filter(identifier %in% results.IT.sac.tmp$identifier)


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


for(z in unique(set$identifier)) {
  results.mean[[length(results.mean)+1]] <- aprioriAnalysis.mean(z, targetTissue, targetLevel)
  results.sac[[length(results.mean)+1]] <- aprioriAnalysis.sac(z, targetTissue, targetLevel)
}

## Process DPI1-5 mean results
results.BT.mean.tib <- bind_rows(results.mean, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.BT.mean.tmp <- results.BT.mean.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.BT.mean.sig <- results.BT.mean.tib %>% filter(identifier %in% results.BT.mean.tmp$identifier)

## Process day of sacrifice results
results.BT.sac.tib <- bind_rows(results.sac, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.BT.sac.tmp <- results.BT.sac.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.BT.sac.sig <- results.BT.sac.tib %>% filter(identifier %in% results.BT.sac.tmp$identifier)


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


for(z in unique(set$identifier)) {
  results.mean[[length(results.mean)+1]] <- aprioriAnalysis.mean(z, targetTissue, targetLevel)
  results.sac[[length(results.mean)+1]] <- aprioriAnalysis.sac(z, targetTissue, targetLevel)
}

## Process DPI1-5 mean results
results.BG.mean.tib <- bind_rows(results.mean, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.BG.mean.tmp <- results.BG.mean.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.BG.mean.sig <- results.BG.mean.tib %>% filter(identifier %in% results.BG.mean.tmp$identifier)

## Process day of sacrifice results
results.BG.sac.tib <- bind_rows(results.sac, .id = "column_label") %>%
  select(-column_label) %>%
  rename("Est" = 1, "pval" = 2) %>%
  mutate(identifier = unique(set$identifier))

results.BG.sac.tmp <- results.BG.sac.tib %>%
  mutate(adj.p.value = p.adjust(pval, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.1)

results.BG.sac.sig <- results.BG.sac.tib %>% filter(identifier %in% results.BG.sac.tmp$identifier)


#### Bind everything up ####
BT.tib <- length(unique(results.BT.tib$identifier))
BG.tib <- length(unique(results.BG.tib$identifier))
IT.tib <- length(unique(results.IT.tib$identifier))
IG.tib <- length(unique(results.IG.tib$identifier))

BT.sig <- length(unique(results.BT.sig$identifier))
BG.sig <- length(unique(results.BG.sig$identifier))
IT.sig <- length(unique(results.IT.sig$identifier))
IG.sig <- length(unique(results.IG.sig$identifier))

datTable <- data.frame("Tissue" = c("Bursa", "Bursa", "Ileum", "Ileum"), "Level" = c("Transcript", "Gene", "Transcript", "Gene"), "Analyzed" = c(BT.tib, BG.tib, IT.tib, IG.tib), "Significant" = c(BT.sig, BG.sig, IT.sig, IG.sig))

remove(annot)
remove(targets)
remove(cnt.bursa.gene)
remove(cnt.bursa.trans)
remove(cnt.gene)
remove(cnt.ileum.gene)
remove(cnt.ileum.trans)
remove(covars)
remove(dge.bursa.gene)
remove(dge.bursa.trans)
remove(dge.ileum.gene)
remove(dge.ileum.trans)
remove(lcpm.bursa.gene)
remove(lcpm.bursa.tmp)
remove(lcpm.bursa.trans)
remove(lcpm.ileum.gene)
remove(lcpm.ileum.tmp)
remove(lcpm.ileum.trans)
remove(lcpm.trans)
save.image("aPrioriGenes-Corr-nolog.Rws")
