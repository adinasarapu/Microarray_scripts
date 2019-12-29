# Title: Differential expression analysis using samr r package

# Date: 12/28/2019
# Author: Ashok Dinasarapu Ph.D
# Email: darphd@gmail.com

# samr: Significance analysis of microarrays
# http://compbio.uthsc.edu/microarray/lecture2.htm

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("affy","samr",hgu133plus2.db")

# Load affy, samr libraries
library(affy)
library(samr)

# base directory
base.dir = "/Users/dinasarapuar/nida"

# Top level directory of CEL files 
cel.dir <- paste(base.dir, "/NIDA_Data", sep="")

# sub directry
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")

# sample name with group details from 88 subjects 
meta.file <- paste(results.base.dir, "data_objects_20191001/DEG/targets.txt", sep="/")
samples <- read.delim(meta.file, header = TRUE, sep = "\t", na.strings = "", fill = TRUE, quote = "")

# Function - collapse ids for search
collapse_pids <- function(pids){
  return (paste(as.character(pids), collapse="|"))
}

# Function - Normalizes data with 'rma' function
exportRawData <- function(your_ids,target_ids,eSet){
  file_pattern <- paste(as.character(your_ids), collapse="|")
  selected.columns <- grep(file_pattern, target_ids, ignore.case = TRUE, perl = FALSE, value = FALSE)
  my.raw.data <- exprs(eSet)
  my.selected.raw.data <- my.raw.data[,selected.columns]
  df.raw.data <- data.frame(my.selected.raw.data)
  return (df.raw.data)
}

# Function - downloads probe set ids to gene ids annotation
getEntrezSymbol <- function(){
  library(hgu133plus2.db) # load database
  #ls("package:hgu133plus2.db") # Provides a summary
  hguEntrezID <- hgu133plus2ENTREZID
  mapped_entrez_keys <- mappedkeys(hguEntrezID)
  probe_entrez <- as.list(hguEntrezID[mapped_entrez_keys])
  hguSymbol <- hgu133plus2SYMBOL
  mapped_symbol_keys <- mappedkeys(hguSymbol)
  probe_symbol <- as.list(hguSymbol[mapped_symbol_keys])
  probe_entrez_symbol <- data.frame(unlist(probe_entrez),unlist(probe_symbol))
  dim(probe_entrez_symbol)
  names(probe_entrez_symbol) <- c("entrez","symbol")
  return (probe_entrez_symbol)
}

# Get all required Affy .CEL files in subdirectories of NIDA_Data directory
all.CEL.files <- dir(cel.dir, recursive=TRUE, full.names=TRUE, pattern="\\.CEL$")
selected_files <- grep(paste(samples$pid,collapse="|"), all.CEL.files, ignore.case = TRUE, perl = FALSE, value = TRUE)

avail.CEL.data <- ReadAffy(filenames = as.vector(selected_files))
eset <- rma(avail.CEL.data) 

length(selected_files) # 88

# Save Raw/RMA expression at Proble level 
raw.normalized <- exportRawData(samples$pid, samples$pid, eset)
colnames(raw.normalized) <- samples$pid

# Get Probe set to Gene ID annotations
probeEntrezSymbol <- getEntrezSymbol()

# select two groups for differential expression analysis

# 1:  "Healthy.Ctrl"        N=25
# 2:  "Healthy.ALC.THC.TP"  N=9

# 1: "VL.pos.Ctrl"          N=8

# 1: "VL.neg.Ctrl"          N=19
# 2: "VL.neg.THC"           N=8
# 3: "VL.neg.THC.TP"        N=19

ctrl <- "Healthy.Ctrl"
case <- "VL.pos.Ctrl"

dim(samples[samples$Group == ctrl,])
dim(samples[samples$Group == case,])

ctrl_case <- c(ctrl,case)

art_case_ctrl <- samples[samples$Group %in% ctrl_case,]
art_case_ctrl$Category <- NA
art_case_ctrl[art_case_ctrl$Group == ctrl,"Category"] <- "1"
art_case_ctrl[art_case_ctrl$Group == case,"Category"] <- "2"
contrasts <- as.vector(art_case_ctrl$Category, mode = "integer")

test_pids <- as.vector(art_case_ctrl$pid)
raw.data.rma <- raw.normalized[,test_pids]

# nperms: Number of permutations used to estimate false discovery rates
# testStatistic: Either "standard" (t-statistic) or ,"wilcoxon" (Two-sample wilcoxon or Mann-Whitney test).
# regression.method: "standard", (linear least squares) or "ranks" (linear least squares on ranked data). 

sams <- c(1:length(contrasts))
data.list <- list(x=as.matrix(raw.data.rma[,sams]),
                  y=contrasts[sams], 	# replace here
                  geneid=row.names(raw.data.rma),
                  genenames=row.names(raw.data.rma), 
                  logged2=TRUE)

set.seed(12345)

# samr.obj is a list with components
samr.obj <- samr(data.list, resp.type="Two class unpaired",nperms=1000) # 

# http://www.meduniwien.ac.at/msi/biometrie/MIMB/_HOMEPAGE%20FILES/complete%20R%20Code.r
# Returns a table of the FDR and upper and lower cutpoints for various values of delta, for a SAM analysis. 

samr.delta <- samr.compute.delta.table(samr.obj,nvals=50)

# choosing the FDR under 0.05
fdr  = 0.05 
delta.max <- samr.delta[samr.delta[,6] <= fdr,1][1] 

# if(!is.na(delta.max)) then select max delta for which FDR is controlled (90th percentile of FDR < FDRthreshold)
# select min delta for which FDR is controlled (median FDR > 2* FDRthreshold)

delta.min <- samr.delta[samr.delta[,5] > 2*fdr,1]               
delta.min <- as.vector(na.omit(delta.min))
delta.min <- delta.min[length(delta.min)]
length.delta.table <- sum(-1*(is.na(samr.delta[, 5]))+1)
dels <- delta.min+(delta.max-delta.min)/length.delta.table*(0:length.delta.table)

# use rough delta table to begin
delta.table.fine <- samr.compute.delta.table(samr.obj, dels=dels)   

# select delta for which FDR is controlled (median percentile of FDR < 0.05)
delta <- delta.table.fine[delta.table.fine[,5] < fdr,1][1]   
delta

samr.siggenes <- samr.compute.siggenes.table(samr.obj,delta, data.list, samr.delta, min.foldchange=0, compute.localfdr=TRUE) # delta = 0.87
# combine both up and down regulated genes

samr.siggenes.probes <- NULL
samr.siggenes.probes <- as.data.frame(rbind(samr.siggenes$genes.up[,c(3,7)],samr.siggenes$genes.lo[,c(3,7)]))
names(samr.siggenes.probes) <- c("probe","fc")

# probe set level output
# out_file_sigprobes <- paste(paste(paste(case,ctrl,sep="_vs_"),round(delta, digits=2) ,sep="_"),"Delta_FDR5perc_Probe.txt",sep="")
# out.file.probes <- paste(results.base.dir, paste0("/data_objects_20191001/DEG/",out_file_sigprobes), sep="")
# write.table(samr.siggenes.probes, file=out.file.probes, sep="\t", quote=FALSE, col.names=T, row.names=F)

row.names(samr.siggenes.probes) <- samr.siggenes.probes[,1]
samr.siggenes.probes[,1] <- NULL

sig.merged <- merge(samr.siggenes.probes,probeEntrezSymbol, by="row.names")
names(sig.merged) <- c("probe","fc","entrez","symbol")
dim(sig.merged)

df.reg <- data.frame(reg = rep(NA,dim(sig.merged)[1]))
assign.reg <- cbind(sig.merged,df.reg)
assign.reg$reg[as.vector(assign.reg$fc) >= 1.3] <- "up"
assign.reg$reg[as.vector(assign.reg$fc) <= (1/1.3)] <- "down"  
assigned.reg <- assign.reg[!is.na(assign.reg$reg),]
assigned.reg$probe <- NULL 
assigned.reg

gene.reg <- aggregate(assigned.reg$reg~assigned.reg$entrez+assigned.reg$symbol,assigned.reg,unique)
names(gene.reg) <- c("entrez","symbol","reg.type")
gene.reg

assigned.reg$fc <-as.numeric(as.character(assigned.reg$fc))
mean.fc <- data.frame(aggregate(assigned.reg$fc, by = list(assigned.reg$entrez), data = assigned.reg, FUN=mean))
names(mean.fc) <- c("entrez","meanFC")
merged.fc.reg.list <- merge(gene.reg,mean.fc, by="entrez")	
merged.fc.reg.list

ret.list <- merged.fc.reg.list[!(merged.fc.reg.list$meanFC > 1 & merged.fc.reg.list$meanFC < 1),]

# excluded genes, due to ambigious direction of regulation
ret.list.exclude <- ret.list[!ret.list[,3] %in% c("down","up"),]
ret.list.exclude

ret.list.include <- NULL
ret.list.include <- ret.list[ret.list[,3] %in% c("down","up"),]
ret.list.include$reg.type <- sapply(ret.list.include$reg.type, as.factor)
ret.list.include

print(paste0("TOTAL = ",dim(ret.list.include[ret.list.include$reg.type %in% c("up","down"),])[1]))
print(paste0("UP = ",dim(ret.list.include[ret.list.include$reg.type == "up",])[1]))
print(paste0("DOWN = ",dim(ret.list.include[ret.list.include$reg.type == "down",])[1]))

out_file <- paste(paste(paste(case,ctrl,sep="_vs_"),round(delta, digits=2) ,sep="_"),"Delta_FDR5perc_FC1.3_Gene.txt",sep="")
out.file.name <- paste(results.base.dir, paste0("/data_objects_20191001/DEG/",out_file), sep="")
write.table(ret.list.include, file=out.file.name, sep="\t", quote=FALSE, col.names=T, row.names=F)

sessionInfo()
