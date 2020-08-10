library(corrr)

dat <- lcpm.all %>%
  filter(levelGT == "gene",
         identifier == "DN10010_c2_g1",
         tissue == "Bursa") %>%
  select(lcpm, virus.sac)

correl <- cor.test(dat$virus.sac, dat$lcpm)

lcpm.all %>%
  filter(levelGT == "gene",
         identifier == "DN10010_c2_g1",
         tissue == "Bursa") %>%
  ggplot(aes(x = virus.sac, y = lcpm)) +
  geom_point(size = 2) +
  ylab("Counts per millions") +
  xlab("Viral titer on day of sacrifice") +
  theme_bw(base_size = 12) +
  labs(title=paste0(annotation[1], " - ", target),
       subtitle=paste0(targetTissue, " - ", targetLevel,
                       "; Pearson's Rho = ", round(correl$estimate, 3),
                       "; p-value = ", correl$p.value))

correl <- cor.test(dat$log.virus.sac, dat$lcpm)

aprioriAnalysis("DN10010_c2_g1", "Bursa", "gene")


for(z in 1:length(xloc)) {
  correl = cor.test(Bcpm[rownames(Bcpm) == xloc[z],], Bcovars$sacdayshed)
  plot(Bcovars$sacdayshed, Bcpm[rownames(Bcpm) == xloc[z],], main = paste("STAT1", z, "of 2"), xlab = "Viral Titer on Day of Sacrifice", ylab = "Counts per Million")
  legend("topright", legend = c(paste("Pearson's rho =", round(correl$estimate, 4)), paste("p-value =", correl$p.value)), bty = "n")
}





##### Making TopGO figures #####
topGO.mappings <- readMappings("G:/Shared drives/RNA Seq Supershedder Project/BWTE DE manuscript/extData/GOdbTrans.txt", sep = "\t", IDsep = ",")

