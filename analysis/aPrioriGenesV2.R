library(tidyverse)
library(edgeR)
library(broom)

setwd("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/")

#### Load data ####
annot <- read.delim("./extData/Trinotate.csv", header = TRUE, sep = "\t")
cnt.trans <- read.table("./extData/rsem.isoform.counts.matrix", header = TRUE)
cnt.gene <- read.table("./extData/rsem.gene.counts.matrix", header = TRUE)
covars <- read.csv("./extData/BWTE54_SSgroup_Raw_Pool.csv", header = TRUE)
targets <- read_delim("./extData/aPrioriTranscripts_V2.csv", delim = ",") %>%
  filter(include == "Yes" | include == "yes")


#### Clean data ####
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
  filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "transcript", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac) %>%
  mutate(log.virus.sac = log(virus.sac+1))


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
  filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "gene", identifier = gene) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac) %>%
  mutate(log.virus.sac = log(virus.sac+1)) %>%
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


#### Analysis functions ####
aprioriAnalysis <- function(target, targetTissue, targetLevel, ...) {
  result.tmp <- lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue) %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    tidy(.) %>%
    mutate(identifier = target) %>%
    filter(term == "log.virus.sac") %>%
    select(-term, -std.error, -statistic) %>%
    pivot_longer(cols = estimate:p.value,
                 names_to = "result.type",
                 values_to = "results.value")

  lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue) %>%
    lm(lcpm ~ log.virus.sac, data = .) %>%
    glance(.) %>%
    mutate(identifier = target) %>%
    select(r.squared, adj.r.squared, identifier) %>%
    pivot_longer(cols = r.squared:adj.r.squared,
                 names_to = "result.type",
                 values_to = "results.value") %>%
    rbind(result.tmp)
}


#### Run analysis loop ####
targetTissue <- "Ileum" ## Bursa or Ileum
targetLevel <- "gene" ## gene or transcript
set <- lcpm.all %>%
  filter(levelGT == targetLevel, tissue == targetTissue) %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

results <- list()

for(z in unique(set$identifier)) {
  results[[length(results)+1]] <- aprioriAnalysis(z, targetTissue, targetLevel)
}

results.tib <- bind_rows(results, .id = "column_label") %>%
  select(-column_label)

results.tmp <- results.tib %>%
  filter(result.type == "p.value") %>%
  mutate(adj.p.value = p.adjust(results.value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05)

results.sig <- results.tib %>% filter(identifier %in% results.tmp$identifier)


#### Plotting ####
results.tib %>%
  filter(result.type == "p.value") %>%
  ggplot(aes(x = results.value)) +
  #geom_histogram(bins = 30) +
  geom_density(fill = "blue") +
  ylab("Frequency") +
  xlab("Regression p-value") +
  theme_classic() +
  labs(title="P-value distribution",
       subtitle="Bursa - Gene level")


aprioriPlotting <- function(target, targetTissue, targetLevel, ...) {
  statResults <- results.tib %>%
    filter(identifier == target)

  annotation <- annot.all %>%
    filter(transcript_id == target | gene_id == target) %>%
    filter(sprot_geneName_BlastX != ".") %>%
    select(sprot_geneName_BlastX) %>%
    slice(1)

  lcpm.all %>%
    filter(levelGT == targetLevel,
           identifier == target,
           tissue == targetTissue) %>%
    ggplot(aes(x = log.virus.sac, y = lcpm)) +
    geom_point() +
    ylab("Log(Counts per million") +
    xlab("Log(Viral titer)") +
    geom_smooth(method = "lm") +
    theme_classic() +
    labs(title=paste0(annotation[1], " - ", target),
         subtitle=paste0(targetTissue, " - ", targetLevel,
                         "; p = ", round(statResults[4,3], 5),
                         "; adj R2 = ", round(statResults[2,3], 3)))
}
