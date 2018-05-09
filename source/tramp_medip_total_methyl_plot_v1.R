# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/12/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_tramp_medip_total_methyl_plot_v1.txt")

require(data.table)
require(ggplot2)

# Load data----
dt1 <- fread("data/tramp_peaks_anno_cpg-10102017.csv")

# X+Y vs. X-Y
a <- dt1$value_1 + dt1$value_2
b <- dt1$value_1 - dt1$value_2
plot(a ~ b)

dt.in <- subset(dt1, 
                insideFeature == "inside")

tmp <- data.frame(gene = dt.in$mgi_symbol)

tmp$gene_lenght_x <- dt.in$end_position.x - dt.in$start_position.x
hist(tmp$gene_lenght_x, 100)

# tmp$gene_lenght_y <- dt.in$end_position.y - dt.in$start_position.y
# hist(tmp$gene_lenght_y, 100)

tmp$pct_in_x <- dt.in$distancetoFeature/tmp$gene_lenght_x
hist(tmp$pct_in_x)

# ??
# tmp$pct_in_y <- dt.in$distancetoFeature/tmp$gene_lenght_y
# hist(tmp$pct_in_y)

tmp$val1 <- dt.in$value_1
tmp$val2 <- dt.in$value_2

summary(tmp)

plot(tmp$val2 ~ tmp$pct_in_x)
points(tmp$val1 ~ tmp$pct_in_x,
       col = "red")

m1 <- loess(tmp$val1 ~ tmp$pct_in_x, span = 0.3)
m2 <- loess(tmp$val2 ~ tmp$pct_in_x, span = 0.3)
dt2 <- data.table(x = tmp$pct_in_x,
                  y1 = m1$fitted,
                  y2 = m2$fitted)
dt2 <- dt2[order(dt2$x)]
plot(dt2$y2 ~ dt2$x, type = "l", ylim = c(0, 25))
lines(dt2$y1 ~ dt2$x, col = "red")

# Upstream (-3kb?)----
unique(dt1$insideFeature)
dt.up <- subset(dt1, 
                insideFeature == "upstream")
summary(dt.up$distancetoFeature)
table(dt.up$distancetoFeature >= -3000)
length(unique(dt.up$mgi_symbol))
dt.up <- droplevels(subset(dt.up,
                           distancetoFeature >= -3000))

tmp1 <- data.frame(gene = dt.up$mgi_symbol,
                   gene_lenght_x = NA)

tmp1$pct_in_x <- dt.up$distancetoFeature/3000
hist(tmp1$pct_in_x)

tmp1$val1 <- dt.up$value_1
tmp1$val2 <- dt.up$value_2

summary(tmp1)

plot(tmp1$val2 ~ tmp1$pct_in_x)
points(tmp1$val1 ~ tmp1$pct_in_x,
       col = "red")

m1 <- loess(tmp1$val1 ~ tmp1$pct_in_x, span = 0.3)
m2 <- loess(tmp1$val2 ~ tmp1$pct_in_x, span = 0.3)
dt3 <- data.table(x = tmp1$pct_in_x,
                  y1 = m1$fitted,
                  y2 = m2$fitted)
dt3 <- dt3[order(dt3$x)]
plot(dt3$y2 ~ dt3$x, type = "l")
lines(dt3$y1 ~ dt3$x, col = "red")

# Downstream (3kb?)----
dt.dn <- subset(dt1, 
                insideFeature == "downstream")
summary(dt.dn$distancetoFeature)
table(dt.dn$distancetoFeature <= 3000)

dt.dn <- droplevels(subset(dt.dn,
                           distancetoFeature <= 3000))

tmp2 <- data.frame(gene = dt.dn$mgi_symbol,
                   gene_lenght_x = NA)

tmp2$pct_in_x <- dt.dn$distancetoFeature/3000 + 1
hist(tmp2$pct_in_x)

tmp2$val1 <- dt.dn$value_1
tmp2$val2 <- dt.dn$value_2

summary(tmp2)

plot(tmp2$val2 ~ tmp2$pct_in_x)
points(tmp2$val1 ~ tmp2$pct_in_x,
       col = "red")

m1 <- loess(tmp2$val1 ~ tmp2$pct_in_x, span = 0.3)
m2 <- loess(tmp2$val2 ~ tmp2$pct_in_x, span = 0.3)
dt4 <- data.table(x = tmp2$pct_in_x,
                  y1 = m1$fitted,
                  y2 = m2$fitted)
dt4 <- dt4[order(dt4$x)]
plot(dt4$y2 ~ dt4$x, type = "l")
lines(dt4$y1 ~ dt4$x, col = "red")

# Combine----
dt.all <- rbind(tmp, tmp1, tmp2)
plot(log2(dt.all$val2 + 1) ~ dt.all$pct_in_x,
     main = "Methylation (%) By Location",
     xaxt = "none",
     xlab = "Location",
     ylab = "log2(Readout + 1)")
points(log2(dt.all$val1 + 1) ~ dt.all$pct_in_x,
       col = "red")
axis(side = 1,
     at = -1:2,
     labels = c("3kb Upstream",
                "TSS",
                "TTS",
                "3kb Downstream"))

m1 <- loess(dt.all$val1 ~ dt.all$pct_in_x, span = 0.25)
m2 <- loess(dt.all$val2 ~ dt.all$pct_in_x, span = 0.25)
dt5 <- data.table(x = dt.all$pct_in_x,
                  y1 = m1$fitted,
                  y2 = m2$fitted)
ggData <- dt5[order(dt5$x)]


ggData <- melt.data.table(data = ggData,
                          id.vars = "x",
                          measure.vars = c("y1", "y2"),
                          variable.name = "Treatment",
                          value.name = "Methylation(%)")
ggData$Treatment <- factor(ggData$Treatment,
                           levels = c("y1", 
                                      "y2"),
                           labels = c("C57 at 24 Weeks",
                                      "TRAMP at 24 Weeks"))

ggplot(ggData,
       aes(x = x,
           y = `Methylation(%)`,
           group = Treatment,
           colour = Treatment)) +
  geom_line() +
  scale_x_continuous("Region",
                     breaks = -1:2,
                     labels = c("3kb Upstream",
                                "TSS",
                                "TTS",
                                "3kb Downstream"))

# Remove zeros----
summary(dt.all)
m1 <- loess(dt.all$val1[dt.all$val1 > 0] ~ dt.all$pct_in_x[dt.all$val1 > 0], 
            span = 0.5)
m2 <- loess(dt.all$val2[dt.all$val2 > 0] ~ dt.all$pct_in_x[dt.all$val2 > 0],
            span = 0.5)
dt5 <- rbindlist(list(data.table(x = m1$x,
                                 y = m1$fitted,
                                 grp = 1),
                      (data.table(x = m2$x,
                                  y = m2$fitted,
                                  grp = 2))))
names(dt5) <- c("x", "y", "grp")
dt5$grp <- factor(dt5$grp, 
                  levels = 1:2,
                  labels = c("C57 at 24 Weeks",
                             "TRAMP at 24 Weeks"))
dt5
ggplot(dt5,
       aes(x = x,
           y = y,
           group = grp,
           colour = grp)) +
  geom_line() +
  scale_x_continuous("Region",
                     breaks = -1:2,
                     labels = c("3kb Upstream",
                                "TSS",
                                "TTS",
                                "3kb Downstream"))