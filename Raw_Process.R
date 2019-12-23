# Install affy packages

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("affy")

# Load affy library
library(affy)

####################
# MAIN DIRECTORIES #
####################

base.dir = "/Users/dinasarapuar/nida"
cel.dir <- paste(base.dir, "/NIDA_Data", sep="")
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")

##################################
# META DATA - Sample/PID details #
##################################
## Self-report (Alcohol) and Toxicology assay (THC and TP)
meta.file.name <- paste(results.base.dir, "data_objects_20191001/Raw_Process/20191001_NIDA_metadata.txt", sep="/")
meta <- read.delim(meta.file.name, header = TRUE, sep = "\t", na.strings = "", fill = TRUE, quote = "")
dim(meta)

# HEALTHY
healthy <- as.vector(meta[meta$Subject == "Healthy",])
healthy_control <- as.vector(healthy[healthy$CTRL_ALC == "CTRL" & healthy$THC == "NO" & healthy$TP == "NO" & healthy$Other == "NO","SampleID"])
healthy_alcohol <- as.vector(healthy[healthy$CTRL_ALC == "ALC" & healthy$THC == "NO" & healthy$TP == "NO" & healthy$Other == "NO","SampleID"])
healthy_thc <- as.vector(healthy[healthy$CTRL_ALC %in% c("CTRL","ALC","NO") & healthy$THC == "THC" & healthy$TP == "NO" & healthy$Other == "NO","SampleID"])
healthy_tp <- as.vector(healthy[healthy$CTRL_ALC %in% c("CTRL","ALC","NO")  & healthy$THC == "NO" & healthy$TP == "TP" & healthy$Other == "NO","SampleID"])
healthy_thc_tp <- as.vector(healthy[healthy$CTRL_ALC %in% c("NO") & healthy$THC == "THC" & healthy$TP == "TP" & healthy$Other == "NO","SampleID"])
healthy_alc_thc_tp <- as.vector(healthy[healthy$CTRL_ALC %in% c("ALC")  & healthy$THC == "THC" & healthy$TP == "TP" & healthy$Other == "NO","SampleID"])
healthy_other <- as.vector(healthy[healthy$Other %in% c("Amphetamine","Oxycodone","DXM/Fluoxetine"),"SampleID"])

# VL<50
VL.neg <- as.vector(meta[meta$Subject == "VL<50",])
VL.neg_control <- as.vector(VL.neg[VL.neg$CTRL_ALC == "CTRL" & VL.neg$THC == "NO" & VL.neg$TP == "NO" & VL.neg$Other == "NO","SampleID"])
VL.neg_alcohol <- as.vector(VL.neg[VL.neg$CTRL_ALC == "ALC" & VL.neg$THC == "NO" & VL.neg$TP == "NO" & VL.neg$Other == "NO","SampleID"])
VL.neg_alc_thc <- as.vector(VL.neg[VL.neg$CTRL_ALC == "ALC" & VL.neg$THC == "THC" & VL.neg$TP == "NO" & VL.neg$Other == "NO","SampleID"])
VL.neg_thc <- as.vector(VL.neg[VL.neg$CTRL_ALC %in% c("CTRL","NO") & VL.neg$THC == "THC" & VL.neg$TP == "NO" & VL.neg$Other == "NO","SampleID"])
VL.neg_tp <- as.vector(VL.neg[VL.neg$CTRL_ALC %in% c("NO") & VL.neg$THC == "NO" & VL.neg$TP == "TP" & VL.neg$Other == "NO","SampleID"])
VL.neg_alc_tp <- as.vector(VL.neg[VL.neg$CTRL_ALC == "ALC" & VL.neg$THC == "NO" & VL.neg$TP == "TP" & VL.neg$Other == "NO","SampleID"])
VL.neg_thc_tp <- as.vector(VL.neg[VL.neg$CTRL_ALC %in% c("NO")  & VL.neg$THC == "THC" & VL.neg$TP == "TP" & VL.neg$Other == "NO","SampleID"])
VL.neg_alc_thc_tp <- as.vector(VL.neg[VL.neg$CTRL_ALC %in% c("ALC")  & VL.neg$THC == "THC" & VL.neg$TP == "TP" & VL.neg$Other == "NO","SampleID"])
VL.neg_other <- as.vector(VL.neg[VL.neg$Other %in% c("Amphetamine","Oxycodone","DXM/Fluoxetine"),"SampleID"])

# VL<50
VL.pos <- as.vector(meta[meta$Subject == "VL>50",])
VL.pos_control <- as.vector(VL.pos[VL.pos$CTRL_ALC == "CTRL" & VL.pos$THC == "NO" & VL.pos$TP == "NO" & VL.pos$Other == "NO","SampleID"])
VL.pos_alcohol <- as.vector(VL.pos[VL.pos$CTRL_ALC == "ALC" & VL.pos$THC == "NO" & VL.pos$TP == "NO" & VL.pos$Other == "NO","SampleID"])
VL.pos_thc <- as.vector(VL.pos[VL.pos$CTRL_ALC %in% c("CTRL","ALC","NO") & VL.pos$THC == "THC" & VL.pos$TP == "NO" & VL.pos$Other == "NO","SampleID"])
VL.pos_tp <- as.vector(VL.pos[VL.pos$CTRL_ALC %in% c("CTRL","ALC","NO") & VL.pos$THC == "NO" & VL.pos$TP == "TP" & VL.pos$Other == "NO","SampleID"])
VL.pos_thc_tp <- as.vector(VL.pos[VL.pos$CTRL_ALC %in% c("NO")  & VL.pos$THC == "THC" & VL.pos$TP == "TP" & VL.pos$Other == "NO","SampleID"])
VL.pos_alc_thc_tp <- as.vector(VL.pos[VL.pos$CTRL_ALC %in% c("ALC")  & VL.pos$THC == "THC" & VL.pos$TP == "TP" & VL.pos$Other == "NO","SampleID"])
VL.pos_other <- as.vector(VL.pos[VL.pos$Other %in% c("Amphetamine","Oxycodone","DXM/Fluoxetine","Methylphenidate"),"SampleID"])

# 145 subjects
total_healthy+total_vl_neg+total_vl_pos

################################################
# Collapse PIDs to get indexes for grep search #
################################################
# Collapse function
collapse_pids <- function(pids){
  return (paste(as.character(pids), collapse="|"))
}

#################
# Read CEL data #
#################
# Get all required Affy .CEL files in subdirectories of NIDA_Data directory
all.CEL.files <- dir(cel.dir, recursive=TRUE, full.names=TRUE, pattern="\\.CEL$")

# exclude HIV (no ART samples)
hiv_only <- c("050197","120152","050194","040088","170050","010141","070226","160067","010215","050203","040159","050213")
hiv_files <- grep(paste(hiv_only,collapse="|"), all.CEL.files, ignore.case = TRUE, perl = FALSE, value = TRUE)

avail.CEL.files <- all.CEL.files[!all.CEL.files %in% hiv_files]

# RMA - Normalizes data with 'rma' function and assigns them to ExpressionSet object
#Background correcting
#Normalizing
#Calculating Expression
# Features 54675; Samples 145

avail.CEL.data <- ReadAffy(filenames = as.vector(avail.CEL.files))
eset <- rma(avail.CEL.data) 
dim(eset) 


# process file names 
# targets - lists the analyzed file names
targets <- pData(eset)
sample_names <- gsub(" ", "",row.names(targets)) 
target.pids <- substr(sample_names, 
                      regexpr("now_", sample_names, ignore.case = TRUE)+4,
                      regexpr("GRF", sample_names,	ignore.case = TRUE)-1)

targets[,2] <- target.pids
row.names(targets) <- NULL
names(targets) <- c("Group","pid")
targets
dim(targets) # 88

# Get indexes using pattern search on targets data frame

hlt_idx_ctrl <- grep(collapse_pids(healthy_control), targets[,2])
hlt_idx_alc <- grep(collapse_pids(healthy_alcohol), targets[,2])
hlt_idx_thc <- grep(collapse_pids(healthy_thc), targets[,2])
hlt_idx_tp <- grep(collapse_pids(healthy_tp), targets[,2])
hlt_idx_thc_tp <- grep(collapse_pids(healthy_thc_tp), targets[,2])
hlt_idx_alc_thc_tp <- grep(collapse_pids(healthy_alc_thc_tp), targets[,2])
hlt_idx_other <- grep(collapse_pids(healthy_other), targets[,2])

VL.neg_idx_ctrl <- grep(collapse_pids(VL.neg_control), targets[,2])
VL.neg_idx_alc <- grep(collapse_pids(VL.neg_alcohol), targets[,2])
VL.neg_idx_thc <- grep(collapse_pids(VL.neg_thc), targets[,2])
VL.neg_idx_alc_thc <- grep(collapse_pids(VL.neg_alc_thc), targets[,2])
VL.neg_idx_tp <- grep(collapse_pids(VL.neg_tp), targets[,2])
VL.neg_idx_alc_tp <- grep(collapse_pids(VL.neg_alc_tp), targets[,2])
VL.neg_idx_thc_tp <- grep(collapse_pids(VL.neg_thc_tp), targets[,2])
VL.neg_idx_alc_thc_tp <- grep(collapse_pids(VL.neg_alc_thc_tp), targets[,2])
VL.neg_idx_other <- grep(collapse_pids(VL.neg_other), targets[,2])

VL.pos_idx_ctrl <- grep(collapse_pids(VL.pos_control), targets[,2])
VL.pos_idx_alc <- grep(collapse_pids(VL.pos_alcohol), targets[,2])
VL.pos_idx_thc <- grep(collapse_pids(VL.pos_thc), targets[,2])
VL.pos_idx_tp <- grep(collapse_pids(VL.pos_tp), targets[,2])
VL.pos_idx_thc_tp <- grep(collapse_pids(VL.pos_thc_tp), targets[,2])
VL.pos_idx_alc_thc_tp <- grep(collapse_pids(VL.pos_alc_thc_tp), targets[,2])
VL.pos_idx_other <- grep(collapse_pids(VL.pos_other), targets[,2])

# update targets df with group name

targets[hlt_idx_ctrl,1] <- "Healthy.Ctrl"
targets[hlt_idx_alc,1] <- "Healthy.ALC"
targets[hlt_idx_thc,1] <- "Healthy.THC"
targets[hlt_idx_tp,1] <- "Healthy.TP"
targets[hlt_idx_thc_tp,1] <- "Healthy.THC.TP"
targets[hlt_idx_alc_thc_tp,1] <- "Healthy.ALC.THC.TP"
targets[hlt_idx_other,1] <- "Healthy.Other"

targets[VL.neg_idx_ctrl,1] <- "VL.neg.Ctrl"
targets[VL.neg_idx_alc,1] <- "VL.neg.ALC"
targets[VL.neg_idx_thc,1] <- "VL.neg.THC"
targets[VL.neg_idx_alc_thc,1] <- "VL.neg.ALC.THC"
targets[VL.neg_idx_tp,1] <- "VL.neg.TP"
targets[VL.neg_idx_alc_tp,1] <- "VL.neg.ALC.TP"
targets[VL.neg_idx_thc_tp,1] <- "VL.neg.THC.TP"
targets[VL.neg_idx_alc_thc_tp,1] <- "VL.neg.ALC.THC.TP"
targets[VL.neg_idx_other,1] <- "VL.neg.Other"

targets[VL.pos_idx_ctrl,1] <- "VL.pos.Ctrl"
targets[VL.pos_idx_alc,1] <- "VL.pos.ALC"
targets[VL.pos_idx_thc,1] <- "VL.pos.THC"
targets[VL.pos_idx_tp,1] <- "VL.pos.TP"
targets[VL.pos_idx_thc_tp,1] <- "VL.pos.THC.TP"
targets[VL.pos_idx_alc_thc_tp,1] <- "VL.pos.ALC.THC.TP"
targets[VL.pos_idx_other,1] <- "VL.pos.Other"

targets <- targets[,c("pid","Group")]

merged_meta <- merge(targets, meta, by.x="pid",by.y="SampleID",all.x = TRUE)
merged_meta <- merged_meta[match(targets$pid, merged_meta$pid),]
head(merged_meta)

# function to get raw data for the selected pids

exportRawData <- function(your_ids,target_ids,eSet){
  #your_ids <- target_pids
  #target_ids <- target_pids
  #eSet <- eset
  file_pattern <- paste(as.character(your_ids), collapse="|")
  selected.columns <- grep(file_pattern, target_ids, ignore.case = TRUE, perl = FALSE, value = FALSE)
  my.raw.data <- exprs(eSet)
  my.selected.raw.data <- my.raw.data[,selected.columns]
  df.raw.data <- data.frame(my.selected.raw.data)
  # colnames(df.raw.data)
  return (df.raw.data)
}

# Save Raw/RMA expression at Proble level 
# 54675   145
raw.rma.export <- exportRawData(merged_meta$pid, targets$pid, eset)
colnames(raw.rma.export) <- merged_meta$pid
dim(raw.rma.export) 

out.file.raw <- paste(results.base.dir, "data_objects_20191001/Raw_Process/Raw_RMA_Proble_Expression.txt", sep="/")
write.table(raw.rma.export, file=out.file.raw, sep="\t", quote=FALSE, col.names=NA, row.names=T)

sessionInfo()