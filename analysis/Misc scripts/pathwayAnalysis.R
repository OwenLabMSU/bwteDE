### First cut at pathway analyses ###
### Run a DE analysis script first ###

library(clusterProfiler)

cp.dat <- annot %>%
  select(
    transcript_id,
    sprot_Top_BLASTX_hit,
    Kegg) %>%
  rename(transcript = transcript_id) %>%
  separate(sprot_Top_BLASTX_hit, into = c("sprot1", NA, NA, NA, NA, "sprot2", NA), "\\^") %>%
  separate(sprot2, sep = "=", into = c(NA, "SwissProt_GeneFunction")) %>%
  separate(sprot1, sep = "_", into = c("SwissProt_GeneName", NA)) %>%
  separate(Kegg, into = c("KEGG_ID", "KO_ID"), sep = "\\`") %>%
  distinct(transcript, SwissProt_GeneName, KEGG_ID, .keep_all = TRUE) %>%
  select(transcript, SwissProt_GeneName, SwissProt_GeneFunction, KEGG_ID, KO_ID) %>%
  distinct(transcript, .keep_all = TRUE) %>%
  filter(grepl("KEGG:", KEGG_ID)) %>%
  mutate(KEGG_ID = str_remove(KEGG_ID, "KEGG:"))


id2gene <- cp.dat %>% select(KEGG_ID, transcript) #TERM2GENE
id2name <- cp.dat %>% dplyr::select(KEGG_ID, SwissProt_GeneName) #TERM2NAME

x <- enricher(dt.tib$transcript, TERM2GENE=id2gene, TERM2NAME=id2name)
x@result
keggIDs <- x@result$ID

### Calls enriched as the ones that are rare enough in the overall dataset. If there's only 1 instance of that term in a pool of 69,053,
### what is the likelihood that it is in the DE set?

library(KEGGREST)
query1 <- keggGet(keggIDs)

query1[[1]]$BRITE %>% as.data.frame()
query1[[1]]$PATHWAY

query1[[2]]$BRITE %>% as.data.frame()
query1[[2]]$PATHWAY

query1[[3]]$BRITE %>% as.data.frame
query1[[3]]$PATHWAY

query1[[4]]$BRITE %>% as.data.frame()
query1[[4]]$PATHWAY

query1[[5]]$BRITE %>% as.data.frame()
query1[[5]]$PATHWAY

query1[[6]]$BRITE %>% as.data.frame()
query1[[6]]$PATHWAY

query1[[7]]$BRITE %>% as.data.frame()
query1[[7]]$PATHWAY

query1[[8]]$BRITE %>% as.data.frame()
query1[[8]]$PATHWAY

query1[[9]]$BRITE %>% as.data.frame()
query1[[9]]$PATHWAY

query1[[10]]$BRITE %>% as.data.frame()
query1[[10]]$PATHWAY

