args <- commandArgs(trailingOnly=TRUE)


### arguments
helper_fn = args[1] # './helper_functions.R'
working_dn = args[2] # "/home/au/code/DTMwork_new/Benchmarking/deconv_benchmark"
dataset = args[3] # "Han" or "Baron"
to_remove = args[4] # 'none', 'help', or something from unique(full_phenoData$cellType)
setwd(working_dn)

### Read data and metadata
cat("\n dataset=",dataset,"\n")
data_orig = as.matrix(read.csv(list.files(path = dataset, pattern = "csv.gz", full.names = TRUE),row.names = 1))
full_phenoData = read.csv(list.files(path = dataset, pattern = "phenoData", full.names = TRUE))

if(to_remove=='help'){
  cat('available cell types are:\n')
  print(unique(full_phenoData$cellType))
  cat('\n aborting...')
}
 


data = data_orig
print(paste0("remove ",to_remove))
transformation = "none" # args[2]
deconv_type = "sc" # args[3]
methods = c("MuSiC","SCDC") # args[6]
nbatch = 5
nmix = 3
nz_thresh = 0.001 #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
number_cells = 1000 # round(as.numeric(args[7]), digits = -2) #has to be multiple of 100
num_cores = 1 # min(5,parallel::detectCores()-1)
normalization_scC = "column" # eg TMM args[4]
normalization_scT = "LogNormalize" # eg TMM args[5]
set.seed(1)

paramstring = paste0("dataset.",dataset,"-remove.",to_remove,"-nz.",nz_thresh,"-transformation.",transformation,"-ncells.",number_cells,"-nmix.",nmix,"-nbatch.",nbatch,"-normalizationC.",normalization_scC,"-normalizationT.",normalization_scT)

#-------------------------------------------------------
### external code
source(helper_fn)


#-------------------------------------------------------
### QC 
require(dplyr); require(Matrix)

# First: cells with library size, mitochondrial or ribosomal content further than three MAD away were discarded
filterCells <- function(filterParam){
  cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
  cellsToRemove
}

libSizes <- colSums(data)
gene_names <- rownames(data)

mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mtPercent <- colSums(data[mtID, ])/libSizes
rbPercent <- colSums(data[rbID, ])/libSizes

lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
  unlist() %>% 
  unique() -> cellsToRemove

if(length(cellsToRemove) != 0){
  data <- data[,-cellsToRemove]
  full_phenoData <- full_phenoData[-cellsToRemove,]
}

# Keep only "detectable" genes: at least 5% of cells (regardless of the group) have a read/UMI count different from 0
keep <- which(Matrix::rowSums(data > 0) >= round(nz_thresh * ncol(data)))
data = data[keep,]

#-------------------------------------------------------

### Data split into training/test 
require(limma); require(dplyr); require(pheatmap)

original_cell_names = colnames(data)
colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])

# Keep CTs with >= 50 cells after QC
cell_counts = table(colnames(data))
to_keep = names(cell_counts)[cell_counts >= 50]
pData <- full_phenoData[full_phenoData$cellType %in% to_keep,]
to_keep = which(colnames(data) %in% to_keep)   
data <- data[,to_keep]
original_cell_names <- original_cell_names[to_keep]

data_clean = data
losses <- data.frame(RMSE=numeric(),
                     Pearson=numeric(),
                     method=character(),
                     batchid=numeric(),
                     stringsAsFactors=FALSE) 

for (b in 1:nbatch){
  
  # Data split into train & test  
  training <- as.numeric(unlist(sapply(unique(colnames(data_clean)), function(x) {
    sample(which(colnames(data_clean) %in% x), cell_counts[x]/2) })))
  testing <- which(!1:ncol(data_clean) %in% training)
  
  # Generate phenodata for reference matrix C
  pDataC = pData[training,]
  
  train <- data_clean[,training]
  test <- data_clean[,testing]

  train_cellID = train
  colnames(train_cellID) = original_cell_names[training]

  cellType <- colnames(train)
  group = list()
  for(i in unique(cellType)){ 
    group[[i]] <- which(cellType %in% i)
  }
  C = lapply(group,function(x) Matrix::rowMeans(train[,x])) #C should be made with the mean (not sum) to agree with the way markers were selected
  C = round(do.call(cbind.data.frame, C))
  # write.table(C, file = "C",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE,)
  
  refProfiles.var = lapply(group,function(x) train[,x])
  refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
  refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
  rownames(refProfiles.var) <- rownames(train)

  #-------------------------------------------------------
  cellType <- colnames(test)
  colnames(test) <- original_cell_names[testing]
  
  generator <- Generator(sce = test, phenoData = full_phenoData, Num.mixtures = nmix, pool.size = number_cells)
  T <- generator[["T"]]
  P <- generator[["P"]]
  
  #-------------------------------------------------------
  ### Transformation, scaling/normalization:
  
  T = Transformation(T, transformation)
  C = Transformation(train_cellID, transformation)
  
  T = Scaling(T, normalization_scT)
  C = Scaling(C, normalization_scC)
  
  #If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
  if(to_remove != "none"){
    
    T <- T[,P[to_remove,] != 0]
    C <- C[,pDataC$cellType != to_remove]
    P <- P[!rownames(P) %in% to_remove, colnames(T)]
    pDataC <- pDataC[pDataC$cellType != to_remove,]
    
  }
  print("writing sc files")
  datapath = file.path("data_batch",dataset,"sc",paramstring,b)
  dir.create(datapath,recursive = TRUE)
  fn_T = file.path(datapath,"T.csv")
  fn_C = file.path(datapath,"C.csv")
  fn_P = file.path(datapath,"P.csv")
  fn_pDataC = file.path(datapath,"pDataC.csv")
  fn_refProfilesVar = file.path(datapath,"refProfilesVar.csv")
  write.csv(T, fn_T, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
  write.csv(C, fn_C, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
  write.csv(P, fn_P, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
  write.csv(pDataC, fn_pDataC, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
  write.csv(refProfiles.var, fn_refProfilesVar, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
  for (m in methods){
    savepath = file.path(datapath,m)
    dir.create(savepath,recursive=TRUE)
    results.list = Deconvolution(T = T, C = C, method = m, phenoDataC = pDataC, P = P, elem = to_remove, refProfiles.var = refProfiles.var)
    print(paste("writing results file"))
    loss = results.list$results %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4), 
                                                     Pearson=cor(observed_values,expected_values) %>% round(.,4))
    fn_loss = file.path(savepath,"loss.txt")
    fn_results.nomelt = file.path(savepath,"nomelt.csv")
    fn_results.results = file.path(savepath,"results.csv")
    write.csv(results.list$nomelt, fn_results.nomelt, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
    write.csv(results.list$results, fn_results.results, quote=FALSE,row.names=TRUE,col.names=TRUE,sep=",")
    write.csv(loss, fn_loss, quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
    loss$method=m
    loss$batchid=b
    losses = bind_rows(losses,loss)      
    print(losses)
  }
}


