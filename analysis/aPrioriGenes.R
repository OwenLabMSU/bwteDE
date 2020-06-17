library(tidyverse)
library(edgeR)
library(broom)

setwd("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/")

## Load data
annot <- read.delim("./extData/Trinotate.csv", header = TRUE, sep = "\t")
cnt.trans <- read.table("./extData/rsem.isoform.counts.matrix", header = TRUE)
cnt.gene <- read.table("./extData/rsem.gene.counts.matrix", header = TRUE)
covars <- read.csv("./extData/BWTE54_SSgroup_Raw_Pool.csv", header = TRUE)

## Clean data
cnt.bursa.gene <- cnt.gene %>%
  select(
    -alignEstimateAbundance_BWTE_Ileum_36_S50,
    -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
    -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("gene") %>%
  separate(gene, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  column_to_rownames("gene") %>%
  select(contains("Bursa"))

cnt.ileum.gene <- cnt.gene %>%
  select(
    -alignEstimateAbundance_BWTE_Ileum_36_S50,
    -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
    -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("gene") %>%
  separate(gene, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  column_to_rownames("gene") %>%
  select(contains("Ileum"))

cnt.bursa.trans <- cnt.trans %>%
  select(
    -alignEstimateAbundance_BWTE_Ileum_36_S50,
    -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
    -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("transcript") %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", "isoform")) %>%
  unite(transcript, pt1, pt2, gene, isoform) %>%
  column_to_rownames("transcript") %>%
  select(contains("Bursa"))

cnt.ileum.trans <- cnt.trans %>%
  select(
    -alignEstimateAbundance_BWTE_Ileum_36_S50,
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

#annot.all <- annot %>%
#  select(transcript_id,
#         gene_id,
#         sprot_Top_BLASTX_hit,
#         sprot_Top_BLASTP_hit,
#         gene_ontology_BLASTX,
#         gene_ontology_BLASTP,
#         Kegg,
#         eggnog,
#         Pfam) %>%
#  separate(sprot_Top_BLASTX_hit, into = c("sprot_geneName_BlastX", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
#  separate(sprot2, sep = "=", into = c(NA, "sprot_geneFunction_BlastX")) %>%
#  separate(sprot_geneFunction_BlastX, sep = ";", into = c("sprot_geneFunction_BlastX", NA)) %>%
#  separate(sprot_Top_BLASTP_hit, into = c("sprot_geneName_BlastP", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
#  separate(sprot2, sep = "=", into = c(NA, "sprot_geneFunction_BlastP")) %>%
#  separate(sprot_geneFunction_BlastP, sep = ";", into = c("sprot_geneFunction_BlastP", NA)) %>%
#  separate(gene_ontology_BLASTX, sep = "\\`", into = paste("GO_BlastX", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
#  separate(gene_ontology_BLASTP, sep = "\\`", into = paste("GO_BlastP", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
#  separate(Pfam, sep = "\\`", into = paste("Pfam", 1:5, sep = "_"), extra = "drop", fill = "right") %>%
#  as_tibble() %>%
#  pivot_longer(cols = 3:23, names_to = "dataSource", values_to = "annotation") %>%
#  drop_na(annotation) %>%
#  filter(annotation != ".")

### Calculate log(CPM) and assemble master DFs
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
  replace(., is.na(.), "0")

# Gene
lcpm.bursa.tmp <- lcpm.bursa.gene %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene")
lcpm.ileum.tmp <- lcpm.ileum.gene %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene")

lcpm.gene <- lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0")

## Gene query database
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

## targets
targets <- read_delim("./extData/aPrioriTranscripts_V2.csv", delim = ",") %>%
  filter(include == "Yes" | include == "yes")


#for(gene in unique(targets$Search_term)) {
#  name <- paste0("results.", gene)
#  assign(name, annot.all %>%
#           filter_at(vars(c(contains("geneName"), contains("geneFunction"),
#                            contains("Kegg"), contains("eggnog"))),
#                          any_vars(str_detect(., fixed(gene, ignore_case = TRUE)))) %>%
#           distinct(transcript_id, .keep_all = TRUE) %>%
#           mutate(searchTerm = gene))
#}


#annot.all %>%
#         filter_all(any_vars(str_detect(., fixed("CD225", ignore_case = TRUE)))) %>%
#         distinct(transcript_id, .keep_all = TRUE)

#for(gene in unique(targets$Search_term)) {
#  name <- paste0("results.", gene)
#  assign(name, annot.all %>%
#           filter_all(any_vars(str_detect(., fixed(gene, ignore_case = TRUE)))) %>%
#           distinct(transcript_id, .keep_all = TRUE) %>%
#           mutate(searchTerm = gene))
#}

#aPrioriCandidatesv1 <- do.call(rbind, lapply(ls(pattern = "results."), get)) %>%
#  distinct(transcript_id, .keep_all = TRUE)


#write.csv(aPrioriCandidatesv1, "aPrioriTranscripts_V2.csv")

### Run regressions on queries
tmp <- lcpm.trans %>%
  filter(transcript %in% targets$transcript_id) %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  filter(group == "I1" | group == "C1") %>%
  select(transcript, tissue, bird, lcpm, virus.sac) %>%
  mutate(log.virus.sac = log(virus.sac+1))


#### Split by tissue. function(target, tissue, ... {}

transcriptAnalysis <- function(target, tissueTarget, ...) {
  result_1 <- tmp %>%
    filter(transcript == target, tissue == tissueTarget) %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    tidy(.) %>%
    mutate(transcript = target) %>%
    filter(term == "log.virus.sac") %>%
    select(-term, -std.error, -statistic) %>%
    pivot_longer(cols = estimate:p.value,
                 names_to = "result.type",
                 values_to = "results.value")

  tmp %>%
    filter(transcript == target, tissue == tissueTarget) %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    glance(.) %>%
    mutate(transcript = target) %>%
    select(r.squared, adj.r.squared, transcript) %>%
    pivot_longer(cols = r.squared:adj.r.squared,
                 names_to = "result.type",
                 values_to = "results.value") %>%
    rbind(result_1)
}

transcriptAnalysis("DN231641_c0_g1_i1")

results <- list()

for(z in unique(tmp$transcript)[4]) {
  print(z)
  results[[length(results)+1]] <- transcriptAnalysis(z, "Bursa")
}

results.tib <- bind_rows(results, .id = "column_label") %>%
  select(-column_label)

results.tmp <- results.tib %>%
  filter(result.type == "p.value") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05)

results.sig <- results.tib %>% filter(transcript %in% results.tmp$transcript)

length(unique(results.sig$transcript))


annot.sig <- annot.all %>%
  filter(transcript_id %in% results.sig$transcript) %>%
  distinct(transcript_id, .keep_all = TRUE)


##### Plotting
plotLCPM <- function(target, ...) {
  tmp %>% filter(transcript == target) %>%
    ggplot(aes(x = log.virus.sac, y = lcpm)) +
    geom_point() +
    geom_smooth(method = "lm") +
    ylab("log(Counts per million)") +
    xlab("log(Viral titer") +
    theme_bw()
}

plotLCPM("DN231641_c0_g1_i1")
