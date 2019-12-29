The following microarray (affy) data analysis pipeline is intended to help my collaborators at NIH/NIAID and Duke University. Session info provides working environment like the version of R used and the loaded packages. `sessioninfo()` gives a snapshot to help people reproduce what you’ve done or spot any anomalies.  

# I). Process Affy raw data - normalization, background correction, expression quantification  

## 1. Install `R` and `RStudio`
## 2. Install `affy` and `hgu133plus2.db` R packages

if (!requireNamespace("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
BiocManager::install(c("affy","hgu133plus2.db"))  

## 3. Update raw affy data directory and meta data file paths  

Create a directory `Raw_Process`  

base.dir = "/Users/dinasarapuar/nida"  
cel.dir <- paste(base.dir, "/NIDA_Data", sep="")  
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")  
meta.file.name <- paste(results.base.dir, "data_objects_20191001/Raw_Process/20191001_NIDA_metadata.txt", sep="/")  
out.file.raw <- paste(results.base.dir, "data_objects_20191001/Raw_Process/Raw_RMA_Proble_Expression.txt", sep="/")  

## 4. Run the script

> sessionInfo()  
R version 3.6.0 (2019-04-26)  
Platform: x86_64-apple-darwin15.6.0 (64-bit)  
Running under: macOS Mojave 10.14.6  

Matrix products: default  
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib  
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib  

locale:  
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8  

attached base packages:  
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base  

other attached packages:  
[1] hgu133plus2cdf_2.18.0 affy_1.62.0           Biobase_2.44.0        BiocGenerics_0.30.0  

loaded via a namespace (and not attached):  
[1] Rcpp_1.0.2            AnnotationDbi_1.46.0  rstudioapi_0.10       IRanges_2.18.1        zlibbioc_1.30.0        
[6] bit_1.1-14            rlang_0.4.0           blob_1.2.0            tools_3.6.0           DBI_1.0.0              
[11] yaml_2.2.0            bit64_0.9-7           digest_0.6.20         tibble_2.1.3          preprocessCore_1.46.0  
[16] crayon_1.3.4          affyio_1.54.0         BiocManager_1.30.4    S4Vectors_0.22.0      vctrs_0.2.0            
[21] zeallot_0.1.0         memoise_1.1.0         RSQLite_2.1.2         compiler_3.6.0        pillar_1.4.2           
[26] backports_1.1.4       stats4_3.6.0          pkgconfig_2.0.2  

# II). Differential expression analysis - from Affy raw data   

## 1. Install `R` and `RStudio`  
## 2. Install `affy` and `hgu133plus2.db` R packages   

if (!requireNamespace("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
BiocManager::install(c("affy","samr","hgu133plus2.db"))  

## 3. Update raw affy data directory, meta data file and output file paths

Create a directory `DEG`

base.dir = "/Users/dinasarapuar/nida"
cel.dir <- paste(base.dir, "/NIDA_Data", sep="")
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")
meta.file <- paste(results.base.dir, "data_objects_20191001/DEG/targets.txt", sep="/")  
out.file.name <- paste(results.base.dir, paste0("/data_objects_20191001/DEG/",out_file), sep="")

## 4. Update the following parameters  

Comparison group names based on metadata i.e `targets.txt` file   

ctrl <- "Healthy.Ctrl"  
case <- "VL.pos.Ctrl"
  
fdr  = 0.05   

## 5. Run the script  

## 6. Pathway analysis  
Upload differentially expressed genes at IPA (Ingenuity pathway analysis) for pathway enrichment analysis.  

> sessionInfo()  
R version 3.6.0 (2019-04-26)  
Platform: x86_64-apple-darwin15.6.0 (64-bit)  
Running under: macOS Mojave 10.14.6  

Matrix products: default  
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib  
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib  

locale:  
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8  

attached base packages:  
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:  
[1] hgu133plus2.db_3.2.3  org.Hs.eg.db_3.8.2    AnnotationDbi_1.46.0  IRanges_2.18.1        S4Vectors_0.22.0      hgu133plus2cdf_2.18.0 samr_3.0               
[8] affy_1.62.0           Biobase_2.44.0        BiocGenerics_0.30.0  

loaded via a namespace (and not attached):  
[1] zip_2.0.3             Rcpp_1.0.2            pillar_1.4.2          compiler_3.6.0        later_0.8.0           BiocManager_1.30.10   tools_3.6.0            
[8] zlibbioc_1.30.0       zeallot_0.1.0         bit_1.1-14            digest_0.6.20         memoise_1.1.0         RSQLite_2.1.2         jsonlite_1.6           
[15] preprocessCore_1.46.0 tibble_2.1.3          pkgconfig_2.0.2       rlang_0.4.0           openxlsx_4.1.0.1      shiny_1.3.2           DBI_1.0.0              
[22] rstudioapi_0.10       yaml_2.2.0            vctrs_0.2.0           fs_1.3.1              shinyFiles_0.7.3      bit64_0.9-7           impute_1.58.0          
[29] R6_2.4.0              GSA_1.03.1            blob_1.2.0            magrittr_1.5          backports_1.1.4       promises_1.0.1        htmltools_0.3.6        
[36] matrixStats_0.54.0    mime_0.7              xtable_1.8-4          httpuv_1.5.1          crayon_1.3.4          affyio_1.54.0          

# III). PCA plot using plot_ly R package  

## 1. Install `R` and `RStudio`  
## 2. Install `plotly` R package  

if (!requireNamespace("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
BiocManager::install("plotly") 

## 3. Update gene expression and meta data file paths  

base.dir = "/Users/dinasarapuar/nida"  
results.base.dir <- paste(base.dir, "/Analysis_2019", sep="")  

Meta data from 145 subjects  
meta.file <- paste(results.base.dir, "data_objects_20191001/PCA/meta_data.txt", sep="/")  

in.file.raw <- paste(results.base.dir, "data_objects_20191001/PCA/Raw_RMA_Proble_Expression.txt", sep="/")  

## 4. Change the number of highly variable genes, if necessary  

Select 5% of genes with highest CV.  
cutoff <- 0.05  

## 5. Change `groups_order` based on meta data column  

## 6. Change color panle  

Dr. Sleasman suggested color panel  

## 7. Run the script  

> sessionInfo()  
R version 3.6.0 (2019-04-26)  
Platform: x86_64-apple-darwin15.6.0 (64-bit)  
Running under: macOS Mojave 10.14.6  

Matrix products: default  
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib  
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib  

locale:  
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8  

attached base packages:  
[1] stats     graphics  grDevices utils     datasets  methods   base  

loaded via a namespace (and not attached):  
[1] compiler_3.6.0 tools_3.6.0    yaml_2.2.0    
