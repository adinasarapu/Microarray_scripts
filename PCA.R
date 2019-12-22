# PCA plot using plotly r package
# https://plotly-r.com/scatter-traces.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("plotly")

library(plotly)

base.dir = "/Users/dinasarapuar/nida"
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")

# Meta data from 145 subjects
meta.file <- paste(results.base.dir, "data_objects_20191001/PCA/meta_data.txt", sep="/")
samples <- read.delim(meta.file, header = TRUE, sep = "\t", na.strings = "", fill = TRUE, quote = "")
dim(samples)

# Normlized raw gene expression data from 145 subjects
in.file.raw <- paste(results.base.dir, "data_objects_20191001/PCA/Raw_RMA_Proble_Expression.txt", sep="/")
raw.normalized <- read.delim(in.file.raw, header = TRUE, check.names = FALSE, sep = "\t", na.strings = "", fill = TRUE, quote = "")
colnames(raw.normalized)[1] <- "ProbeID"
row.names(raw.normalized) <- raw.normalized$ProbeID
raw.normalized$ProbeID <- NULL

##############################
# Find highly variable genes #
##############################

# 1 indicates rows, 2 indicates columns
compute_cv <- function(x) sd(x) / mean(x)
cv <- apply(raw.normalized, 1, compute_cv)
summary(cv)

# Select 5% of genes with highest CV.
cutoff <- 0.05
summary(cv[rank(cv) / length(cv) > 1 - cutoff])
raw.normalized.selected <- raw.normalized[rank(cv) / length(cv) > 1 - cutoff, ]
dim(raw.normalized.selected)

# There are 23 subgroups based on meta data file 
groups_order <- c("Healthy.Ctrl","Healthy.ALC","Healthy.THC","Healthy.TP","Healthy.THC.TP","Healthy.ALC.THC.TP","Healthy.Other",
                  "VL.pos.Ctrl", "VL.pos.ALC", "VL.pos.THC", "VL.pos.TP","VL.pos.THC.TP","VL.pos.ALC.THC.TP","VL.pos.Other",
                  "VL.neg.Ctrl","VL.neg.ALC","VL.neg.THC", "VL.neg.ALC.THC", "VL.neg.TP", "VL.neg.ALC.TP", "VL.neg.THC.TP",
                  "VL.neg.ALC.THC.TP", "VL.neg.Other")

samples$Group <- factor(samples$Group, levels=groups_order) 
samples <- samples[order(match(samples$Group,groups_order)),]
sorted_raw.normalized <- raw.normalized.selected[,samples$pid]

################################
# PCA of highly variable genes #
################################

group_count <- dplyr::count(samples, Group)
as.data.frame(group_count)

# Dr. Sleasman suggested color panel
custom.col <- c(
  rep("#BFBFBF",group_count$n[1]),
  rep("#808080",group_count$n[2]),
  rep("#808080",group_count$n[3]),
  rep("#808080",group_count$n[4]),
  rep("#808080",group_count$n[5]),
  rep("#808080",group_count$n[6]),
  rep("#808080",group_count$n[7]),
  rep("#0432FF",group_count$n[8]),
  rep("#4472C4",group_count$n[9]),
  rep("#4472C9",group_count$n[10]),
  rep("#4472C9",group_count$n[11]),
  rep("#4472C9",group_count$n[12]),
  rep("#4472C9",group_count$n[13]),
  rep("#4472C9",group_count$n[14]),
  rep("#FF0000",group_count$n[15]),
  rep("#FF65DF",group_count$n[16]),
  rep("#FF65DF",group_count$n[17]),
  rep("#FF65DF",group_count$n[18]),
  rep("#FF65DF",group_count$n[19]),
  rep("#FF65DF",group_count$n[20]),
  rep("#FF65DF",group_count$n[21]),
  rep("#FF65DF",group_count$n[22]),
  rep("#FF65DF",group_count$n[23]))

################
# PCA analysis #
################

hiv.pca <- prcomp(t(sorted_raw.normalized), center = TRUE,scale. = TRUE)
summary(hiv.pca)

pca.df <- as.data.frame(hiv.pca$x)
pca.df$group <- samples$Group
pca.df$VL <- samples$Viral_Load
pca.df$color <- custom.col
pca.df$text <- row.names(pca.df)

xaxis <- list(showgrid = T,title = "PC1",titlefont = list(size=18,color="black"),tickfont = list(size=14,color="black"),
              dtick = 20, tickmode = "linear")
yaxis <- list(showgrid = T,title = "PC2",titlefont = list(size=18,color="black"),tickfont = list(size=14,color="black"),
              dtick = 20, tickmode = "linear")
zaxis <- list(showgrid = T,title = "PC3",titlefont = list(size=18,color="black"),tickfont = list(size=14,color="black"),
              dtick = 20, tickmode = "linear")

scene = list(
  xaxis = xaxis,
  yaxis = yaxis,
  zaxis = zaxis)

plot_ly(pca.df, type="scatter3d", 
        # mode = "text",
        x = ~PC1, y = ~PC2, z = ~PC3,
        #text = ~text, textfont = list(color = ~color, size = 9) # showlegend=T, 
        mode = "markers", marker = list(color = ~color, size=7)
) %>% add_markers(marker = list(color = ~color, size=7))  %>% 
  layout(title = "PCA plot", scene = scene)

sessionInfo()