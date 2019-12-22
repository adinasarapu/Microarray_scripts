# Microarray_scripts: PCA plot using plot_ly R package  

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

## 7. Finally run the script  

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
