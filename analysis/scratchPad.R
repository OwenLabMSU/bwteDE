#### Top section of code is for plotting expression levels of given genes/transcripts ####
tmp <- lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0") %>%
  filter(transcript == "DN3985_c0_g1_i5") %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  filter(group == "I1") %>%
  #filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "transcript", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac) %>%
  mutate(log.virus.sac = log10(virus.sac+1))

### Line
lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0") %>%
  filter(transcript == "DN3985_c0_g1_i5") %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  filter(group == "I1") %>%
  #filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "transcript", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac) %>%
  mutate(log.virus.sac = log10(virus.sac+1)) %>%
  filter(tissue == "Ileum") %>%
  unite(birdLabel, bird, SSgroup.virus.sac, sep = "-") %>%
  ggplot(aes(x = log.virus.sac, y = lcpm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = birdLabel), ) +
  ylab("Log2(Counts per million)") +
  xlab("Log10(Viral titer + 1)") +
  theme_classic()


### Box
lcpm.bursa.tmp %>%
  full_join(lcpm.ileum.tmp) %>%
  replace(., is.na(.), "0") %>%
  filter(transcript == "DN25207_c2_g1") %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c(NA, NA, "tissue", "bird", NA)) %>%
  mutate(bird = as.integer(bird)) %>%
  left_join(covars, by = "bird") %>%
  filter(group == "I1") %>%
  #filter(group == "I1" | group == "C1") %>%
  mutate(levelGT = "gene", identifier = transcript) %>%
  select(identifier, levelGT, tissue, bird, lcpm, virus.sac, SSgroup.virus.sac) %>%
  mutate(SSgroup.virus.sac = fct_relevel(SSgroup.virus.sac, "LOW", "MODERATE", "HIGH")) %>%
  ggplot(aes(x = SSgroup.virus.sac, y = lcpm)) +
  geom_boxplot() +
  ylab("Log2(Counts per million)") +
  xlab("SSgroup.virus.sac") +
  theme_classic()



#### Using topGO to do GO analysis #####

#########################################
## Create database for gene level data ##
#########################################

goIDs <- read_delim("../extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  filter(gene %in% dt.tib$gene) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  distinct(gene, GO, .keep_all = TRUE) %>%
  group_by(gene) %>%
  mutate(GO.ID = paste0("GO", 1:length(gene))) %>%
  pivot_wider(names_from = GO.ID,
              values_from = GO,
              values_fill = NA) #%>%
  #unite("GOs", GO1:GO293, na.rm = TRUE, sep = ", ")

write.table(goIDs,
            "G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/wegoGene.txt",
            row.names=FALSE,
            na = "",
            quote=FALSE,
            col.names=FALSE,
            sep = "\t")



##########################################
## Create database for trans level data ##
##########################################

goIDs <- read_delim("../extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", "trans")) %>%
  unite(transcript, pt1, pt2, gene, trans) %>%
  filter(transcript %in% dt.tib$transcript) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  group_by(transcript) %>%
  mutate(GO.ID = paste0("GO", 1:length(transcript))) %>%
  pivot_wider(names_from = GO.ID,
              values_from = GO,
              values_fill = NA) #%>%
  #unite("GOs", GO1:GO210, na.rm = TRUE, sep = ", ")

write.table(goIDs,
            "G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/wegoTrans.txt",
            row.names=FALSE,
            quote=FALSE,
            na = "",
            col.names=FALSE,
            sep = "\t")





library(topGO)
topGO.mappings <- readMappings("./TestGO.txt", sep = "\t", IDsep = ",")

DE.results <- tfit$p.value %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene") %>%
  mutate(CtlvI1.adj = p.adjust(CtlvI1, method='fdr', n = nrow(.))) %>%
  filter(CtlvI1.adj < 0.1)


all.genes <- sort(unique(as.character(rownames(tfit$p.value))))
int.genes <- DE.results$gene
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj.BP <- new("topGOdata", ontology='BP',
              allGenes = int.genes,
              annot = annFUN.gene2GO,
              nodeSize = 5,
              gene2GO = topGO.mappings)

go.obj.MF <- new("topGOdata", ontology='MF',
                 allGenes = int.genes,
                 annot = annFUN.gene2GO,
                 nodeSize = 5,
                 gene2GO = topGO.mappings)

go.obj.CC <- new("topGOdata", ontology='CC',
                 allGenes = int.genes,
                 annot = annFUN.gene2GO,
                 nodeSize = 5,
                 gene2GO = topGO.mappings)

results <- runTest(go.obj.BP, algorithm = "elim", statistic = "fisher")

results.tab <- GenTable(object = go.obj.BP,
                        elimFisher = results,
                        topNodes = 33)

##### Create GO CSVs for Amanda  ### GENES
# This is the set of GO IDs for the DE genes
goIDs <- read_delim("../extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  filter(gene %in% dge$genes$gene_id) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  dplyr::rename(GO.ID = GO) %>%
  distinct(gene, GO.ID) %>%
  filter(gene %in% dt.tib$gene)

goIDs.all <- read_delim("../extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  filter(gene %in% dge$genes$gene_id) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  dplyr::rename(GO.ID = GO)

### Annotation info
annot2 <- annot %>%
  dplyr::select(gene_id,
    sprot_Top_BLASTX_hit) %>%
  dplyr::rename(gene = gene_id) %>%
  inner_join(dt.tib) %>%
  separate(sprot_Top_BLASTX_hit, into = c("sprot1", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
  separate(sprot2, sep = "=", into = c(NA, "SwissProt_GeneFunction")) %>%
  separate(sprot1, sep = "_", into = c("SwissProt_GeneName", NA)) %>%
  distinct(gene, SwissProt_GeneName, .keep_all = TRUE) %>%
  dplyr::select(gene, SwissProt_GeneName, SwissProt_GeneFunction) %>%
  filter(SwissProt_GeneName != ".")

#results.tib <- tibble(GO.ID = c("GO:0016823"))

GOdb <- results.tib %>%
  left_join(goIDs)# %>%
  left_join(annot2)

write.csv(GOdb, "ssgroupVirusAvg-Ileum-Gene.GOdb.csv")






##### Create GO CSVs for Amanda  ### TRANSCRIPTS
# This is the set of GO IDs for the DE genes
goIDs <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", "trans")) %>%
  unite(transcript, pt1, pt2, gene, trans) %>%
  filter(transcript %in% dge$genes$transcript_id) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  dplyr::rename(GO.ID = GO) %>%
  distinct(transcript, GO.ID) %>%
  filter(transcript %in% dt.tib$transcript)

goIDs.all <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/goIDs.txt", delim = "\t", col_names = c("transcript", "GO")) %>%
  separate(transcript, into = c(NA, "pt1", "pt2", "gene", NA)) %>%
  unite(gene, pt1, pt2, gene) %>%
  filter(gene %in% dge$genes$gene_id) %>%
  separate(GO, into = paste0("GO", 1:500), sep = ",") %>%
  pivot_longer(cols = starts_with("GO"),
               names_to = "GO.ID",
               values_to = "GO",
               values_drop_na = TRUE) %>%
  dplyr::select(-GO.ID) %>%
  dplyr::rename(GO.ID = GO)

### Annotation info
annot2 <- annot %>%
  dplyr::select(transcript_id,
                sprot_Top_BLASTX_hit) %>%
  dplyr::rename(transcript = transcript_id) %>%
  inner_join(dt.tib) %>%
  separate(sprot_Top_BLASTX_hit, into = c("sprot1", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
  separate(sprot2, sep = "=", into = c(NA, "SwissProt_GeneFunction")) %>%
  separate(sprot1, sep = "_", into = c("SwissProt_GeneName", NA)) %>%
  distinct(transcript, SwissProt_GeneName, .keep_all = TRUE) %>%
  dplyr::select(transcript, SwissProt_GeneName, SwissProt_GeneFunction) %>%
  filter(SwissProt_GeneName != ".")


GOdb <- results.tib %>%
  left_join(goIDs) %>%
  left_join(annot2)

write.csv(GOdb, "ssgroupVirusAvg-Ileum-Trans.GOdb.csv")


##########################################
topGO.mappings <- readMappings("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/GOdbTrans.txt", sep = "\t", IDsep = ",")

  DE.results <- dt.tib %>% filter(!!sym("MvH.I5") != 0)

  all.genes <- sort(unique(as.character(rownames(dge$counts))))
  int.genes <- DE.results$transcript
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  names(int.genes) = all.genes

  go.obj.BP <- new("topGOdata", ontology='BP',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   nodeSize = 5,
                   gene2GO = topGO.mappings)

  go.obj.MF <- new("topGOdata", ontology='MF',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   nodeSize = 5,
                   gene2GO = topGO.mappings)

  go.obj.CC <- new("topGOdata", ontology='CC',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   nodeSize = 5,
                   gene2GO = topGO.mappings)

  results.BP <- runTest(go.obj.BP, algorithm = "elim", statistic = "fisher")

  results.tab.BP <- GenTable(object = go.obj.BP,
                             elimFisher = results.BP,
                             topNodes = 50) %>%
    as_tibble() %>%
    mutate(pVal = as.numeric(elimFisher)) %>%
    filter(pVal < 0.01) %>%
    mutate(Domain = "BP") %>%
    mutate(Comparison = targetComp)

  results.MF <- runTest(go.obj.MF, algorithm = "elim", statistic = "fisher")

  results.tab.MF <- GenTable(object = go.obj.MF,
                             elimFisher = results.MF,
                             topNodes = 50) %>%
    as_tibble() %>%
    mutate(pVal = as.numeric(elimFisher)) %>%
    filter(pVal < 0.01) %>%
    mutate(Domain = "MF") #%>%
    mutate(Comparison = targetComp)

  results.CC <- runTest(go.obj.CC, algorithm = "elim", statistic = "fisher")

  results.tab.CC <- GenTable(object = go.obj.CC,
                             elimFisher = results.CC,
                             topNodes = 50) %>%
    as_tibble() %>%
    mutate(pVal = as.numeric(elimFisher)) %>%
    filter(pVal < 0.01) %>%
    mutate(Domain = "CC") %>%
    mutate(Comparison = targetComp)

  rbind(results.tab.BP,
        results.tab.MF,
        results.tab.CC)


results <- list()

for(z in unique(names(dt.tib)[-1])) {
  if (dt.tib %>%
      filter(!!sym(z) != 0) %>%
      nrow(.) > 0)
    results[[length(results)+1]] <- goEnrich(z)
}

results.tib <- bind_rows(results, .id = "column_label") %>%
  select(-column_label) %>%
  arrange(GO.ID, Comparison)


grep("GO:0000122", topGO.mappings, ignore.case = TRUE)

showSigOfNodes(go.obj.MF, score(results.MF), firstSigNodes = 5, useInfo ='all')
printGraph(go.obj.MF, results.MF, firstSigNodes = 5, fn.prefix = "ssgroupVirusAvg-Ileum-Trans-MF", useInfo = "all", pdfSW = TRUE)
