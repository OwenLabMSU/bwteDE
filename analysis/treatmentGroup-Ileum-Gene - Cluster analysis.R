## Load packages and data
library(limma)
library(edgeR)
library(Glimma)
library(gplots)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(kableExtra)
library(tidyverse)

setwd("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/bwteDE")

## Load data
annot <- read.delim("../extData/Trinotate.csv", header = TRUE, sep = "\t")
cnt <- read.table("../extData/rsem.gene.counts.matrix", header = TRUE)
covars <- read.csv("../extData/BWTE54_SSgroup_Raw_Pool.csv", header = TRUE)

cnt <- cnt %>%
  select(
    -alignEstimateAbundance_BWTE_Ileum_36_S50,
    -alignEstimateAbundance_BWTE_Bursa_36_S31,-alignEstimateAbundance_BWTE_Ileum_19_S35,
    -alignEstimateAbundance_BWTE_Bursa_19_S14) %>%
  rownames_to_column("gene") %>%
  separate(gene, into = c(NA, "pt1", "pt2", "gene")) %>%
  unite(gene, pt1, pt2, gene) %>%
  column_to_rownames("gene") %>%
  select(contains("Ileum"))

covars <- covars %>%
  filter(!bird %in% c("36", "19")) %>%
  arrange(bird) %>%
  separate(group, into = c("status", "dpi")) %>%
  mutate(dpi = as.numeric(dpi)) %>%
  mutate(dpi = formatC(dpi, width = 2, format = "d", flag = "0")) %>%
  unite(group, status, dpi, sep = "")

annot <- annot %>%
  separate(transcript_id, into = c(NA, "pt1", "pt2", "gene", "isoform")) %>%
  unite(transcript_id, pt1, pt2, gene, isoform) %>%
  separate(gene_id, into = c(NA, "pt1", "pt2", "gene")) %>%
  unite(gene_id, pt1, pt2, gene)




bird <- as.factor(covars$bird)
sex <- as.factor(covars$sex)
age <- as.numeric(covars$age)
weight <- covars$wt_55
group <- recode(covars$group, C01 = "Ctl", C14 = "Ctl") %>%
  as.factor()
pool <- as.factor(covars$Pool.Ileum)




#Convert to DGEList object
dge <- DGEList(counts=cnt)
dge$genes <- annot[match(rownames(dge$counts), annot$gene_id),]

#CPM and log-CPM
cpm <- cpm(dge)
lcpm <- cpm(dge, log = TRUE)

#### Retain only genes expressed at >0.5 CPM in at least 25% of individuals
table(rowSums(dge$counts==0) == length(cnt[1,]))
keep.exprs <- rowSums(cpm>0.5) >= length(cnt[1,])/4
sum(keep.exprs)

dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
table(rowSums(dge$counts==0) == length(cnt[1,]))
dim(dge)

#TMM normalization
dge <- calcNormFactors(dge, method = "TMM")




## Establish design matrix
design = model.matrix(~ 0 +
                        group +
                        age +
                        sex +
                        weight +
                        pool)

colnames(design) <- gsub("group", "", colnames(design))

contr.matrix = makeContrasts(
  CtlvI01 = Ctl - I01,
  CtlvI03 = Ctl - I03,
  CtlvI05 = Ctl - I05,
  CtlvI14 = Ctl - I14,
  I01vI03 = I01 - I03,
  I03vI05 = I03 - I05,
  I05vI14 = I05 - I14,
  levels = colnames(design))

v <- voomWithQualityWeights(dge, design, plot = TRUE)



vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
tfit <- treat(vfit, lfc = 1.0)
plotSA(tfit)




dt <- decideTests(tfit, p.value = 0.1, adjust.method = "fdr")
summary(dt) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)

dt.tib <- as_tibble(dt, rownames = NA) %>%
  rownames_to_column("gene") %>%
  mutate_at(vars(starts_with("C")), as.numeric) %>%
  mutate_at(vars(starts_with("I")), as.numeric) %>%
  filter_at(vars(CtlvI01:I05vI14), any_vars(. != 0))




# Get the gene names for DE genes
deGeneCounts <- lcpm %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% dt.tib$gene) %>%
  column_to_rownames("gene")

newNames1 <- colnames(deGeneCounts) %>%
  as_tibble() %>%
  separate(value, into = c(NA, NA, NA, "sample", NA)) %>%
  mutate(group = fct_relevel(group, "Ctl", "I01", "I03", "I05", "I14")) %>%
  #arrange(group) %>%
  unite("sample", 2:1, sep = "_")

colnames(deGeneCounts) <- newNames1$sample
deGeneCounts <- deGeneCounts %>% select(order(colnames(.))) %>%
  as.matrix()


## Set up palette
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)


hc <- hclust(as.dist(1-cor(t(deGeneCounts))))
plot(hc)

# Plot the heatmap
heatmap.2(deGeneCounts,
          Colv = FALSE,
          Rowv = as.dendrogram(hc),
          col = rev(morecols(50)),
          trace = "none",
          colsep = c(9, 21, 33, 45),
          dendrogram = "row",
          main = "Ileum - Gene level",
          margins = c(5, 10),
          #ColSideColors = col.diff,
          scale ="row"
          )


# A complete linkage hierarchical clustering was performed using
# the Pearson's correlation coefficient between samples and
# Euclidian distances between genes with the R package gplots (3.0).

# A partition clustering analysis was done using a gap statistic to
# calculate the optimal number of clusters (k) using clustgap() from
# the R package cluster (2.0.7).

# The optimal k of 9 was used in a partitioning around medoids (pam)
# clustering analysis using pam() from the R package cluster (Figure S3).
# The k clusters were then summarized in 3 major general patterns of expression,
# and the pam clustering was repeated with k = 3.

hc <- hclust(as.dist(1-cor(t(deGeneCounts))))

library(cluster)
mydist <- function(x) as.dist(1-cor(t(x)))
mycluster <- function(x, k) list(cluster=cutree(hclust(mydist(x), method = "average"),k=k))
myclusGap <- clusGap(deGeneCounts,
                     FUN = mycluster,
                     K.max = 25,
                     B = 100)

gstab = data.frame(myclusGap$Tab, k=1:nrow(myclusGap$Tab))

ggplot(gstab, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_errorbar(aes(ymax = gap+SE.sim,
                    ymin = gap-SE.sim)) +
  scale_x_continuous(breaks = round(seq(min(0), max(25), by = 1),1)) +
  theme_classic(base_size = 18)


cluster=cutree(hclust(mydist(deGeneCounts),
                      method = "average"),
               k = 14)

clusterIDs <- cluster %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("id") %>%
  filter(value == 2)

tmp <- deGeneCounts %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("id") %>%
  filter(id %in% clusterIDs$id) %>%
  pivot_longer(cols = Ctl_07:I14_52,
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c("group", "bird")) %>%
  mutate(group = fct_relevel(group, "Ctl", "I01", "I03", "I05", "I14")) %>%
  arrange(group) %>%
  mutate(xVar = paste0(group, "-", bird),
    xVar = factor(xVar, levels = unique(xVar)))


ggplot(data = tmp, aes(y = lcpm, x = xVar, group = id)) +
  geom_line(color = "red", alpha = 0.15) +
  stat_summary(aes(y = lcpm, group = 1),
               fun.y = mean, colour = "red", geom = "line", group = 1, lwd = 2) +
  ylab("Log(Counts per million)") +
  xlab("Bird") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))



deGeneCounts
