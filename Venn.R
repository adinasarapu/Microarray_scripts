library(VennDiagram)
#library(dplyr)
#library(tidyr)

# MAIN DIRECTORIES

base.dir = "/Users/dinasarapuar/nida"
deg.dir <- paste(base.dir, "/Analysis_2019/data_objects_20191001/DEG/genes", sep="")
results_dir <- paste(base.dir, "/Analysis_2019/data_objects_20191001/DEG/VennDiagrams", sep="")

textFileData <- function(file_name){
  data <- read.delim(paste(deg.dir,file_name,sep = "/"), header = TRUE,sep = "\t",
                     na.strings = "",fill = FALSE,quote = "",stringsAsFactors = FALSE)
  return (data )
}
  
VL.pos.Ctrl <- textFileData("original/VL.pos.Ctrl_vs_Healthy.Ctrl_0.81Delta_FDR5perc_FC1.3_Gene.txt")
VL.neg.Ctrl <- textFileData("original/VL.neg.Ctrl_vs_Healthy.Ctrl_0.74Delta_FDR5perc_FC1.3_Gene.txt")
VL.neg.THC <- textFileData("original/VL.neg.THC_vs_Healthy.Ctrl_0.7Delta_FDR5perc_FC1.3_Gene.txt")
VL.neg.THC.TP <- textFileData("original/VL.neg.THC.TP_vs_Healthy.Ctrl_0.76Delta_FDR5perc_FC1.3_Gene.txt")

col <- list()
col[["I"]] <- "#8EA9DB"
col[["II"]] <- "#521B93"
col[["III"]] <- "#e35f09"
col[["IV"]] <- "#fabe07"
col[["V"]] <- "#A9D08E" 

# Venn2 

col_1 <- col[["II"]]
col_2 <- col[["IV"]]

name_1 <- names(col)[2]
name_2 <- names(col)[4]

DF.1 <- VL.pos.Ctrl
DF.2 <- VL.neg.Ctrl

list.1 <- paste(DF.1$entrez, DF.1$symbol, DF.1$reg.type,sep = "_")
list.2 <- paste(DF.2$entrez, DF.2$symbol, DF.2$reg.type,sep = "_")

length(list.1)
length(list.2)

out_venn2 <- paste(results_dir, "/Venn2_II_III_testing.tiff", sep="")
x2_df = list(A = list.1, B = list.2)
venn.plot <- venn.diagram(
  x2_df,
  filename =  out_venn2, # lwd = 4,
  # main = "HIV VL+ vs HIV VL-",
  # category.names = c("ART/VL+","ART/VL-"),
  category.names = c(name_1,name_2),
  # main.cex = 2,
  # scaled = TRUE, cex = c(8,4,3), cat.dist = c(0.045, 0.105),cat.pos = c(-40, 8),
  scaled = FALSE, cex = c(5,5,5), cat.dist = c(0.035, 0.035),cat.pos = c(-6, 6),
  # inverted = TRUE, # color
  # reverse = TRUE,
  rotation.degree=0,
  lwd = 7,
  col = c(col_1,col_2),
  alpha = 1, # color intensity
  # label.col = "white",
  fontfamily = "serif",
  # fontface = "bold",
  cat.cex = 4,
  cat.col = c(col_1,col_2) 
  # hyper.test = TRUE 
  # total.population=17050
)

# Venn3

DF.1 <- VL.neg.Ctrl
DF.2 <- VL.neg.THC
DF.3 <- VL.neg.THC.TP

list.1 <- paste(DF.1$entrez, DF.1$reg.type,sep = "_")
list.2 <- paste(DF.2$entrez, DF.2$reg.type,sep = "_")
list.3 <- paste(DF.3$entrez, DF.3$reg.type,sep = "_")

length(list.1)
length(list.2)
length(list.3)

col_1 <- col[["III"]]
col_2 <- col[["IV"]]
col_3 <- col[["V"]]

name_1 <- names(col)[3]
name_2 <- names(col)[4]
name_3 <- names(col)[5]

out_venn3 <- paste(results_dir, "/Venn3_II_III_IV_testing.tiff", sep="")

x3_df = list(A=list.1, B=list.2, C=list.3)

venn.plot <- venn.diagram(
  x3_df,
  # main = "HIV-, no substance, THC and THC+TP", 
  # main.cex = 2,
  filename = out_venn3,
  category.names = c(name_1,name_2,name_3),
  col = c(col_1,col_2,col_3),
  lwd = 7,  scaled = TRUE,
  #rotation.degree = 180,
  #col = "transparent",
  #fill = c("red", "blue", "green"),
  alpha = 1,
  #label.col = c("darkred", "white", "darkblue", "white","white", "white", "darkgreen"),
  cex = 4.0, # gene number size
  fontfamily = "serif",
  #fontface = "bold",
  #cat.default.pos = "text",
  cat.cex = 3.0, # size of III, IV and V
  cat.fontfamily = "serif",
  #euler.d = TRUE
  cat.dist = c(0.06, 0.06, -0.46),
  cat.pos = 0, margin = 0, cat.col = c(col_1,col_2,col_3)
)

# Venn4

DF.1 <- VL.neg.Ctrl
DF.2 <- VL.neg.THC
DF.3 <- VL.neg.THC.TP
DF.4 <- VL.neg.THC.TP

list.1 <- paste(DF.1$entrez, DF.1$reg.type,sep = "_")
list.2 <- paste(DF.2$entrez, DF.2$reg.type,sep = "_")
list.3 <- paste(DF.3$entrez, DF.3$reg.type,sep = "_")
list.4 <- paste(DF.4$entrez, DF.4$reg.type,sep = "_")

col_1 <- col[["I"]]
col_2 <- col[["II"]]
col_3 <- col[["III"]]
col_4 <- col[["IV"]]

name_1 <- names(col)[2]
name_2 <- names(col)[3]
name_3 <- names(col)[4]
name_4 <- names(col)[5]

out_venn4 <- paste(results_dir, "/Venn4_II__III_IV_V_testing.tiff", sep="")

x4_df = list(A=list.1, B=list.2, C=list.3, D=list.4)

venn4_out <- venn.diagram(x4_df,
          # main = "VL<=50",
          filename = out_venn4,
          alpha = 1,
          category.names = c(name_1,name_2,name_3,name_4),
          col = c(col_1,col_2,col_3,col_4),
          #reverse = TRUE,
          cex = 2, lwd = 5,
          fontfamily = "serif", 
          #fontface = "bold",
          #category.names = c("","","",""),
          #cat.col =  c("black","green","blue","red"), 
          cat.cex = 2.5, # size of I, II, III, IV and V
          cat.pos = c(-5,1,-5,1), cat.dist = c(0.24,0.24,0.13,0.12),
          cat.fontfamily = "serif", 
          rotation.degree = 0,
          margin = 0.2
)

sessionInfo()