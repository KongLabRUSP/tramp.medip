# |-----------------------------------------------------------------------|
# | Project:         TRAMP MeDIP-seq Data Analysis                        |
# | Experiment:                                                           |
# | Compound Number:                                                      |
# | Script:          Analysis of DMR                                      | 
# | Scientists:      Wenji Li,  Davit Sargsyan                            |         
# | Created:         07/31/2017                                           |
# |-----------------------------------------------------------------------|
# Source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Header----
require(data.table)
require(ggplot2)

# Read data from CSV----
# Positive vs. negative controls
dt1 <- fread("data/tramp_peaks_anno_cpg.csv")

# Keep only useful variables----
dt1 <- dt1[, c("mgi_symbol",
               "gene_id",
               "insideFeature",
               "value_1",
               "value_2",
               "log2_fold_change",
               "description")]
dt1

# 1. Remove unmapped regions----
dt1 <- subset(dt1, 
              !(is.na(mgi_symbol) |
                  mgi_symbol == ""))

# 2. Remove rows with values equal to zero----
dt1 <- subset(dt1, 
              value_1 != 0 &
                value_2 != 0)

# 3. Remove rows with less than 2-fold change----
dt1 <- subset(dt1, 
              log2_fold_change < -1 |
                log2_fold_change > 1)
hist(dt1$log2_fold_change, 100)

# Keep only upstream, inside and downstream locations----
dt1 <- droplevels(subset(dt1,
                         insideFeature %in% c("upstream",
                                              "inside",
                                              "downstream")))
dt1$insideFeature <- factor(dt1$insideFeature,
                            levels = c("upstream",
                                       "inside",
                                       "downstream"),
                            labels = c("Promoter",
                                       "Inside",
                                       "Body"))

# Concatenate gene and region IDs----
dt1$id <- paste(dt1$mgi_symbol,
                dt1$gene_id,
                sep = "_")
summary(dt1)

# Order by foldchange----
dt1 <- dt1[order(log2_fold_change)]

# Top 50 upregulated regions
upreg50 <- dt1[1:25, ]
upreg50$reg <- "Upregulated"
# Order by id
upreg50 <- upreg50[order(id)]

# Top 50 downregulated regions----
downreg50 <- dt1[(nrow(dt1) - 24):nrow(dt1)]
downreg50$reg <- "Downregulated"
# Order by id
downreg50 <- downreg50[order(id)]

# Merge
dt2 <- rbindlist(list(upreg50,
                      downreg50))
dt2$id <- factor(dt2$id,
                 levels = rev(dt2$id))
dt2$reg <- factor(dt2$reg,
                  levels = unique(dt2$reg))

# Save to CSV
write.csv(dt2,
          file = "tmp/MeDIP_Wenji_Top50Up_Top50Down.csv",
          row.names = FALSE)

# Heatmap of top 50 upregulated genes
p1 <- ggplot(data = dt2) +
  facet_wrap(~ reg,
             scales = "free_y") +
  geom_tile(aes(x =  insideFeature,
                y = id,
                fill = log2_fold_change),
            color = "black") +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "black", 
                       midpoint = 0, 
                       limit = c(-6, 6), 
                       name="") +
  scale_x_discrete("Location",
                   expand = c(0, 0)) + 
  scale_y_discrete("Gene  and Region",
                   expand = c(0, 0)) +
  ggtitle("Log2(TRAMP/C57) Methylation of Top 25 Upregulated \n and Top 25 Downregulated Regions")  +
  theme(axis.text.x = element_text(angle = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
p1
tiff(filename = "tmp/MeDIP_Wenji_Top25Up_Top25Down.tiff",
     height = 7,
     width = 7,
     units = 'in',
     res = 600,
     compression = "lzw+p")
print(p1)
graphics.off()

# Are ther overlaps in genes?
unique(upreg50$mgi_symbol) %in% unique(downreg50$mgi_symbol) 
# No overlaps