
setwd("/home/au/code/DTMwork_new/Benchmarking/deconv_benchmark/")
rm(list=ls()) 

# ############################
# ## SC, All
# ############################

datasets = c("Baron","Han") # args[1]
for (dataset in datasets){
  ### Read data and metadata
  if (file.exists(paste0(dataset,"/",tolower(dataset),".rds"))){
    data_orig = readRDS(list.files(path = dataset, pattern = "rds", full.names = TRUE))  
    full_phenoData = read.table(list.files(path = dataset, pattern = "phenoData", full.names = TRUE), header=TRUE)
  }else{
    data_orig = as.matrix(read.csv(list.files(path = dataset, pattern = "csv.gz", full.names = TRUE),row.names = 1))
    full_phenoData = read.csv(list.files(path = dataset, pattern = "phenoData", full.names = TRUE))
  }
  
  celltypes = c("none")
  for (to_remove in celltypes){
    rm(data)
    data = data_orig
    print(paste0("remove ",to_remove))
    transformation = "none" # args[2]
    deconv_type = "sc" # args[3]
    methods = c("MuSiC","SCDC") # args[6]
    nbatch = 1 # number of synthetic datasets, each containing nmix mixtures
    nmix = 10 # number of pseudobulk mixtures in each dataset
    nz_thresh = 0.01 #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
    number_cells = 1000 # round(as.numeric(args[7]), digits = -2) #has to be multiple of 100
    num_cores = 1 # min(5,parallel::detectCores()-1)
    normalization_scC = "column" # eg TMM args[4]
    normalization_scT = "LogNormalize" # eg TMM args[5]
    min.percentage = 1 # minimum mixture percentage for any given cell type (this is scaled so it's actually the ratio between min max that matters)
    max.percentage = 99 # maximum mixture percentage for any given cell type (this is scaled so it's actually the ratio between min max that matters)
    CTthresh = 0.01 # only use celltypes making up at least 0.05% of all cells
    CTall = FALSE # force generation of pseudobulk mixtures containing all (post-qc) detected cell types from original data
    set.seed(1)
    
    ############################
    ## Bulk, Baron
    ############################
    
    paramstring = paste0("dataset.",dataset,"-remove.",to_remove,"-nz.",nz_thresh,"-transformation.",transformation,
                         "-ncells.",number_cells,"-nmix.",nmix,"-nbatch.",nbatch,
                         "-minP.",min.percentage,"-maxP.",max.percentage,"-CTthresh.",CTthresh,"-CTall.",CTall,
                         "-normalizationC.",normalization_scC,"-normalizationT.",normalization_scT)
    
    #-------------------------------------------------------
    ### Helper functions + CIBERSORT external code
    source('./helper_functions.R')
    # source('CIBERSORT.R')
    
    
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
    to_keep = names(cell_counts)[cell_counts >= CTthresh*sum(cell_counts)]
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
      rm(training,testing,train,test,T,P,C,refProfiles.var,marker_distrib,loss,pDataC,train_cellID,cellType,group,refProfiles.var,keep,markers,v,fit,fit2,generator,results.list)
      
      # Data split into train & test  
      training <- as.numeric(unlist(sapply(unique(colnames(data_clean)), function(x) {
        sample(which(colnames(data_clean) %in% x), cell_counts[x]/2) })))
      testing <- which(!1:ncol(data_clean) %in% training)
      
      # Generate phenodata for reference matrix C
      pDataC = pData[training,]
      
      train <- data_clean[,training]
      test <- data_clean[,testing]
      
      # "write.table" & "saveRDS" statements are optional, for users willing to avoid generation of matrix C every time:    
      # write.table(pDataC, file = paste(dataset,"phenoDataC",sep="_"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
      
      train_cellID = train
      colnames(train_cellID) = original_cell_names[training]
      # saveRDS(object = train_cellID, file = paste(dataset,"qc_filtered_train.rds",sep="_")) #It has to contain cellID as colnames, not cellType (for scRNA-seq methods)
      # saveRDS(object = test, file = paste(dataset,"qc_filtered_test.rds",sep="_"))
      
      # reference matrix (C) + refProfiles.var from TRAINING dataset
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
      # write.table(refProfiles.var, "refProfiles.var", quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
      
      #-------------------------------------------------------
      #Normalization of "train" followed by marker selection 
      
      #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
      cellType = colnames(train) 
      keep <- sapply(unique(cellType), function(x) {
        CT_hits = which(cellType %in% x)
        size = ceiling(nz_thresh*length(CT_hits))
        Matrix::rowSums(train[,CT_hits,drop=FALSE] != 0) >= size
      })
      train = train[Matrix::rowSums(keep) > 0,]
      train2 = Normalization(train)
      
      # INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
      #[compare one group with average expression of all other groups]
      annotation = factor(colnames(train2))
      design <- model.matrix(~0+annotation)
      colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
      cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
      colnames(cont.matrix) <- colnames(design)
      diag(cont.matrix) <- (ncol(design)-1)/ncol(design)
      
      v <- limma::voom(train2, design=design, plot=FALSE) 
      fit <- limma::lmFit(v, design)
      fit2 <- limma::contrasts.fit(fit, cont.matrix)
      fit2 <- limma::eBayes(fit2, trend=TRUE)
      
      markers = marker.fc(fit2, log2.threshold = log2(2))
      
      #-------------------------------------------------------
      ### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
      cellType <- colnames(test)
      colnames(test) <- original_cell_names[testing]
      
      generator <- Generator(sce = test, phenoData = pData, Num.mixtures = nmix, pool.size = number_cells,
                             min.percentage = min.percentage, max.percentage = max.percentage, CTall = CTall)
      T <- generator[["T"]]
      P <- generator[["P"]]
      #-------------------------------------------------------
      ### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
      
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
      
      # make sure all the genes are in both reference and bulk

      genes_to_keep = sort(intersect(intersect(rownames(C),rownames(T)),rownames(refProfiles.var)))

      C = C[match(genes_to_keep,rownames(C)),]
      T = T[match(genes_to_keep,rownames(T)),]
      refProfiles.var[match(genes_to_keep,rownames(refProfiles.var)),]
      
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
  }
}


