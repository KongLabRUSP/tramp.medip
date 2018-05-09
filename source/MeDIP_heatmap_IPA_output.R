# |-----------------------------------------------------------------------|
# | Project:         TRAMP MeDIP-seq Data Analysis                        |
# | Experiment:                                                           |
# | Compound Number:                                                      |
# | Script:          Analysis of DMR                                      | 
# | Scientists:      Wenji Li,  Davit Sargsyan                            |         
# | Created:         08/10/2017                                           |
# |-----------------------------------------------------------------------|
# Source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Header----
require(data.table)
require(ggplot2)

setwd("MeDIP_TRAMP_Wenji")

# Read data from CSV----
# Positive vs. negative controls
dt1 <- fread("data/all_undo_aug 2017-filtered.csv")

# Remove Location and Type(s)
dt1 <- droplevels(subset(dt1,
                         select = -c(Location,
                                     `Type(s)`)))
dt1 <- unique(dt1)

# 2. Remove rows with less than 4-fold change----
dt1 <- subset(dt1, 
              `Expr Log Ratio` < -2 |
                `Expr Log Ratio` > 2)
hist(dt1$`Expr Log Ratio`, 100)

# Keep only upstream, inside and downstream locations----
table(dt1$Note)
dt1 <- droplevels(subset(dt1,
                         Note %in% c("upstream",
                                     "inside",
                                     "downstream")))
dt1$Note <- factor(dt1$Note,
                            levels = c("upstream",
                                       "inside",
                                       "downstream"),
                            labels = c("Promoter",
                                       "Body",
                                       "Downstream"))

# Concatenate gene and region IDs----
length(unique(dt1$Symbol))
length(unique(dt1$ID))

dt1$id <- paste(dt1$ID,
                dt1$Note,
                sep = "_")
length(unique(dt1$id))
summary(dt1)

# Duplicates
tmp <- dt1[id %in% dt1$id[duplicated(dt1$id)], ]

# NOTE: multiple counts for same regions; note that
setkey(dt1, id)
dt1[, n := 1:.N,
    by = id]
dt1$id <- paste(dt1$id, 
                dt1$n,
                sep = "_")
length(unique(dt1$id))

# Order by foldchange----
dt1 <- dt1[order(`Expr Log Ratio`)]

# Top X upregulated regions
top.upreg <- dt1[1:100, ]
top.upreg$reg <- "Decreased Methylation" 
# Order by id
top.upreg <- top.upreg[order(id)]

# Top X downregulated regions----
top.downreg <- dt1[(nrow(dt1) - 99):nrow(dt1)]  
top.downreg$reg <- "Increased Methylaton"
# Order by id
top.downreg <- top.downreg[order(id)]

# Merge
dt2 <- rbindlist(list(top.upreg,
                      top.downreg))
dt2$id <- factor(dt2$id,
                 levels = rev(dt2$id))
dt2$reg <- factor(dt2$reg,
                  levels = unique(dt2$reg))

# Save to CSV
write.csv(dt2,
          file = "tmp/MeDIP_Wenji_TopUp_TopDown.csv",
          row.names = FALSE)

# Heatmap of top genes
p1 <- ggplot(data = dt2) +
  facet_wrap(~ reg,
             scales = "free_y") +
  geom_tile(aes(x =  Note,
                y = id,
                fill = `Expr Log Ratio`),
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
  ggtitle("Log2(TRAMP/C57) of Top 100 Increased\n and Top 100 Desreased Methylation Regions")  +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")

tiff(filename = "tmp/MeDIP_Wenji_TopUp_TopDown.tiff",
     height = 12,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()