# args <- commandArgs(trailingOnly=TRUE)
# 
# if(length(args)!=9){
#   
#   print("Please check that all required parameters are indicated or are correct")
#   print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
#   print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
#   stop()
# } 

setwd("/home/au/code/DeconvPaper/Benchmarking/cobos_benchmarks")
rm(T,P,C,refProfiles.var,pDataC)

method = "nnls" # args[6]
method = "SCDC" # args[6]

to_remove = "none" # args[8]

deconv_type = "sc"

marker_strategy = "all"

# for sc
datafolder = "/home/au/code/DTMwork/Benchmarking/cobos_benchmarks/data/Baron/sc/MuSiC/Baron-nz.0.001-transformation.none-ncells.1000-nmix.100-normalizationC.LogNormalize-normalizationT.LogNormalize" 
addmethod = file.path(datafolder,method)
dir.create(addmethod,recursive = TRUE)
fn_T = file.path(datafolder,"T.csv")
fn_pDataC = file.path(datafolder,"pDataC.csv")
fn_P = file.path(datafolder,"P.csv")
fn_C = file.path(datafolder,"C.csv")
fn_refProfilesVar = file.path(datafolder,"refProfilesVar.csv")

pDataC = read.csv(fn_pDataC,row.names = 1)
T = read.csv(fn_T,row.names = 1)
P = read.csv(fn_P,row.names = 1)
C = read.csv(fn_C,row.names = 1)
refProfiles.var = read.csv(fn_refProfilesVar,row.names = 1)

#-------------------------------------------------------
### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
if(deconv_type == "bulk"){

  results.list = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, marker_distrib = marker_distrib, refProfiles.var = refProfiles.var) 
} else if (deconv_type == "sc"){
  
  results.list = Deconvolution(T = T, C = C, method = method, phenoDataC = pDataC, P = P, elem = to_remove, refProfiles.var = refProfiles.var) 
}

print(paste("writing results file: ",fn_res))
loss = results.list$results %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
                                                 Pearson=cor(observed_values,expected_values) %>% round(.,4))
fn_loss = file.path(addmethod,"loss.txt")
fn_results.nomelt = file.path(addmethod,"nomelt.csv")
fn_results.results = file.path(addmethod,"results.csv")
write.csv(results.list$nomelt, fn_results.nomelt, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
write.csv(results.list$results, fn_results.results, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
write.csv(loss, fn_loss, quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

print(loss)
