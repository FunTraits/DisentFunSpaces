################################################################################
#  ______                _   _                 
#  |  ___|              | | (_)                
#  | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
#  |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#  | | | |_| | | | | (__| |_| | (_) | | | \__ \
#  \_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#
################################################################################

# This script summarize all R functions used in this project ----

 ### 0. Statistics ----
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}
sesandpvalue = function (obs,rand,nreps,probs=c(0.025,0.975),rnd = 2){
  SES = (obs - mean(rand)) / sd(rand)
  pValsSES = rank(c(obs,rand))[1] / (length(rand) + 1)  
  results = round(c(obs,SES, mean(rand),quantile(rand,prob=probs),pValsSES,nreps), rnd)
  names(results)= c("Observed","SES","MeanRd","CI025Rd","CI975Rd","Pval","Nreps")
  return(results)
}

error_trait = function(ColumnWithNA,phenoOriginal,pheno,phenoAll){
  mean_error = numeric()
  for (na in 1:length(ColumnWithNA)){
    imputed_sp_traits = rownames(phenoOriginal)[is.na(phenoOriginal[,ColumnWithNA[na]])]
    
    error = abs(pheno[imputed_sp_traits,ColumnWithNA[na]] -
                  phenoAll[imputed_sp_traits,ColumnWithNA[na]]) / 
      (max(phenoOriginal[,ColumnWithNA[na]],na.rm=T)+min(phenoOriginal[,ColumnWithNA[na]],na.rm=T))
    
    mean_error[na] = mean(error*100)
  }
  names(mean_error) = ColumnWithNA
  return(mean_error)
}

error_space = function(PCAtraits,PCAtraitsALL){
  x = PCAtraits$PCoA$vectors
  y = PCAtraitsALL$PCoA$vectors
  PrcTest = ade4::procuste.rtest(as.data.frame(x), as.data.frame(y[rownames(x),]), nrepet = 999)
  cor_coefficient = PrcTest$obs
  p_value = PrcTest$pvalue
  significance = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", paste0("ns"))))
  title_with_pearson = paste0(round(cor_coefficient,3), " ", significance)
}

angle = function(x, y){
  norm_x = sqrt(sum(x^2))
  norm_y = sqrt(sum(y^2))
  cosXY = round(as.numeric((x %*% y) / (norm_x * norm_y)), 8)
  #rounding avoids numerical problems with the acos function
  angle = acos(cosXY) * 180 / pi
  return(angle)
}

angleCorrelation = function(all,focal){
  x = as.numeric(as.dist(focal))
  y = as.numeric(as.dist(all[rownames(focal),rownames(focal)]))
  
  names(y) = names(x) = apply(combn(rownames(focal),2),2,function(x){paste0(x[1]," ; ",x[2])})       
  results = list(table = cbind.data.frame(y,x),
                 test = cor.test(y,x))
}

### 1. Data sorting functions. ----

# get
get_taxonomy = function (x, preferred_data_sources = c(11:15), fuzzy = TRUE, 
                         verbose = TRUE) {
  
  if (length(x) > 1) {
    out <- lapply(x, FUN = function(x) {
      out <- tryCatch(get_taxonomy(as.character(x)), error = function(e) e)
      if (inherits(out, "error")) {
        out <- data.frame(scientificName = x)
      }
      return(out)
    })
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      utils::install.packages("dplyr")
    }
    out <- dplyr::bind_rows(out)
  } else {
    if (nchar(gsub(" ", "", x)) <= 1 || is.null(x)) 
      x <- NULL
    resolved <- taxize::gnr_resolve(x, preferred_data_sources = preferred_data_sources, 
                                    best_match_only = TRUE, canonical = TRUE)
    if (length(resolved$matched_name2) == 0) {
      out <- data.frame(user_supplied_name = x)
      attributes(out)$warning <- paste("No matching species name found!")
      out$warnings <- paste(out$warnings, attributes(out)$warning, 
                            sep = "; ")
    }
    else {
      temp <- taxize::get_gbifid_(resolved$matched_name2)[[1]]
      if (any(temp$matchtype == "EXACT") || !fuzzy) 
        temp <- temp[temp$matchtype == "EXACT", ]
      if (all(temp$status == "SYNONYM")) {
        out <- tryCatch(get_taxonomy(temp$species[which.max(temp$confidence)]), 
                        error = function(e) e)
        if (inherits(out, "error")) {
          out <- data.frame(scientificName = x)
        }
        out$synonym = TRUE
        out$user_supplied_name = x
        attributes(out)$warning <- paste("Synonym provided! Automatically set to accepted species Name!")
        out$warnings <- paste(out$warnings, attributes(out)$warning, 
                              sep = "; ")
      }
      else {
        if (any(temp$status %in% c("ACCEPTED", "DOUBTFUL"))) {
          temp <- temp[temp$status %in% c("ACCEPTED", 
                                          "DOUBTFUL"), ]
          out <- temp[which.max(temp$confidence), ]
          if (!requireNamespace("dplyr", quietly = TRUE)) {
            utils::install.packages("dplyr")
          }
          if (is.null(out$species)) {
            out$species <- out$genus
          }
          out <- cbind(scientificName = x, synonym = FALSE, 
                       scientificNameStd = out$species, author = sub(paste0(out$species, 
                                                                            " "), "", out$scientificname), taxonRank = out$rank, 
                       dplyr::select(out, dplyr::one_of("confidence", 
                                                        "kingdom", "phylum", "class", "order", 
                                                        "family", "genus")), taxonomy = "GBIF Backbone Taxonomy", 
                       taxonID = paste0("http://www.gbif.org/species/", 
                                        out$usagekey, ""), warnings = "")
        }
      }
      if (out$synonym[1] & verbose) 
        warning(paste("Synonym provided! Automatically set to accepted species Name!"))
    }
  }
  class(out) <- c("data.frame")
  return(out)
}
# standard_names: R function to standardize species names, by default the names
# sources is IUCN (source = 163)
standard_names = function(species,source = 163,dim=1){
  require(taxize)
  species = gsub(" ","_",species)
  nrowTraits = length(species)
  standNamesTraits = as.data.frame(matrix(NA, nrow = nrowTraits, ncol=1,
                                          dimnames=list(species, "IUCNName")))
  batchSize = 1000
  
  for(i in 1:(ceiling(nrowTraits/batchSize))){
    cat(paste0("\nBatch ", i," out of ", (ceiling(nrowTraits/batchSize)), "\n"))
    rowsSelect = ((i-1) * batchSize + 1):min((i*batchSize), nrowTraits)
    taxoAux = as.data.frame(taxize::gnr_resolve(names = rownames(standNamesTraits)[rowsSelect],
                                                preferred_data_sources = source, ## 163 is IUCN
                                                best_match_only = TRUE,
                                                canonical = TRUE))
    standNamesTraits[taxoAux$user_supplied_name,"IUCNName"] = taxoAux$matched_name2
  }
  standNamesTraits[,"IUCNName"] = gsub(x=standNamesTraits[,"IUCNName"], 
                                        pattern = " ", replacement= "_")
  naIUCN = which(is.na(standNamesTraits[,"IUCNName"]))
  standNamesTraits[naIUCN, "IUCNName"] = rownames(standNamesTraits)[naIUCN]
  
  ### Some species could not be resolved into IUCN naming. Remove them:
  yesName = which(!is.na(standNamesTraits[,"IUCNName"]))
  if(dim == 1){
    return(species[yesName])
  }else{
    return(species[yesName,])
  }
}

# checkSpToAvonet: For birds, the names are checked according to AVONET taxonomy
checkSpToAvonet =function(AVONET_Names,avonet, chechTraits){
  traits=avonet
  traitsp=traits$Species3
  AMNIOTEsp=gsub(" ","_",unique(rownames(chechTraits)))
  zz=cbind(Old=AMNIOTEsp,
           New=AMNIOTEsp)
  
  for (i in 1:nrow(zz)){
    nnam=AVONET_Names[which(AVONET_Names$Species1_BirdLife%in%zz[i,]),]$Species3_BirdTree
    if(length(nnam)>1){
      nnam=names(which.max(table(nnam[which(nnam%in%traitsp)])))
    }
    if(length(nnam)==0){
      nnam=AVONET_Names[which(AVONET_Names$Species2_eBird%in%zz[i,]),]$Species3_BirdTree
      if(length(nnam)>1){
        nnam=names(which.max(table(nnam[which(nnam%in%traitsp)])))
      }
      if(length(nnam)==0){
        nnam=AVONET_Names[which(AVONET_Names$eBird.species.group%in%zz[i,]),]$Species3_BirdTree
        if(length(nnam)>1){
          nnam=names(which.max(table(nnam[which(nnam%in%traitsp)])))
        }
      }
    }
    if(length(nnam)==1){
      zz[i,2]=nnam
    }else{
      zz[i,2]=zz[i,1]
    }
  }
  
  AMNIOTEsp = zz
  chechTraits$Species3=zz[,2]
  chechTraits = unique(chechTraits)
  toRem = chechTraits[which(chechTraits$Species3 %in% names(which(table(chechTraits$Species3) > 1))),]
  avesTraitsOK = chechTraits[!rownames(chechTraits) %in% rownames(toRem[!rownames(toRem) == toRem$Species3,]),]
  rownames(avesTraitsOK) = avesTraitsOK$Species3
  
  return(avesTraitsOK)
  
}

# sortPCAnames
checkPCA =function(PCA,tx){
  if(tx == 1){
    PCA = list(M = PCA[["M"]],L = PCA[["L"]],D = PCA[["D"]], 
               ML = PCA[["ML"]],DL = PCA[["DL"]], MD = PCA[["MD"]],
               MDL = PCA[["MDL"]])
  }
  if(tx == 2){
    PCA = list(L = PCA[["L"]],D = PCA[["D"]],DL = PCA[["DL"]])
  }
  return(PCA)
}
### 2. TPD functions, adapted from 'funspace' packages. ----

densityProfileTPD <- function(x, probs=seq(0, 1, by=0.01)){
  TPDList <- x$TPDc$TPDc
  results <- matrix(NA, nrow=length(TPDList), ncol=length(probs),
                    dimnames = list(names(TPDList), probs))
  for(comm in 1:length(TPDList)){
    TPD <- TPDList[[comm]]
    cellSize <- x$data$cell_volume
    alphaSpace_aux <- TPD[order(TPD, decreasing = T)]
    FRicFunct <- function(TPD, alpha = 1) { ### FUNCTION TO ESTIMATE FUNCTIONAL RICHNESS FROM TPD OBJECT
      FRic <- numeric(length(alpha))
      for(i in 1:length(alpha)){
        TPDAux <- TPD
        greater_prob <- max(alphaSpace_aux[which(cumsum(alphaSpace_aux) >= alpha[i])])
        TPDAux[TPDAux < greater_prob] <- 0
        FRic[i] <- sum(TPDAux>0)* cellSize
      }
      names(FRic) <- alpha
      return(FRic)
    }
    results[comm,] <- FRicFunct(TPD = TPD, alpha = probs)  
  }
  return(results)
}

# make.PCA : Here TRAITS is a table of traits sp x tr : MUST BE SCALED AND log10
make.PCA = function(TRAITS,dimensionsAux,savePCA){
  paranAux = paran(TRAITS)
  PCA_TRAITS = list()
  PCA_TRAITS$traits = TRAITS
  PCA_TRAITS$PCA = princomp(PCA_TRAITS$traits)
  PCA_TRAITS$dimensions = dimensionsAux
  PCA_TRAITS$variance = (summary(PCA_TRAITS$PCA)[1][[1]]^2)[1:dimensionsAux] /
    sum(summary(PCA_TRAITS$PCA)[1][[1]]^2)
  PCA_TRAITS$loadings = PCA_TRAITS$PCA$loadings
  PCA_TRAITS$traitsUse = data.frame(PCA_TRAITS$PCA$scores[, 1:PCA_TRAITS$dimensions]) 
  if(!is.na(savePCA)){
    saveRDS(PCA_TRAITS,savePCA)
  }
  return(PCA_TRAITS)
}

# PCoAGraph: make.PCA adapted to different kind of traits type
# PCoAGraph = function(x,datax,namesx,groups,saveFile){
#   if (!is.list(x)) stop("Input 'x' must be a list of data frames.")
#   if (!is.null(saveFile) && !is.character(saveFile)) stop("'saveFile' must be a valid file path or NA.")
#   ktabList = ktab.list.df(x)
#   disTraits = dist.ktab(ktabList, groups, scan = FALSE,option = c("scaledBYrange")) 
#   PCoATraits = pcoa(disTraits)
#   plotenvfit = vegan::envfit(PCoATraits$vectors, datax,w = 1)
#   spp.scrs = as.data.frame(scores(plotenvfit, display = "vectors"))
#   spp.scrs = cbind(spp.scrs, Species = rownames(spp.scrs))
# 
#   rownames(PCoATraits$vectors) = namesx
#   PcoaCorrel=NULL
#   for (cor in 1:length(x)){
#     PcoaCorrel=rbind(PcoaCorrel,CorrelType(PCoATraits,x[[cor]]))
#   }
#   res=list(PCoA=PCoATraits,varPlot=spp.scrs,PCoACor=PcoaCorrel)
#   
#   if(!is.na(saveFile)){
#     saveRDS(res, saveFile)
#   }
#   
#   return(res)
# }

PCoAGraph <- function(x, datax, namesx, groups, saveFile = NA) {
  # --- Checks ---
  if (!is.list(x) || !all(vapply(x, is.data.frame, logical(1L)))) {
    stop("Argument 'x' must be a list of data frames.")
  }
  if (!is.data.frame(datax)) stop("'datax' must be a data frame.")
  if (length(namesx) != nrow(datax)) {
    stop("Length of 'namesx' must match the number of rows in 'datax'.")
  }
  if (!is.null(saveFile) && !is.character(saveFile) && !is.na(saveFile)) {
    stop("'saveFile' must be a character string or NA.")
  }
  
  # --- Create ktab object ---
  ktabList <- ktab.list.df(x)
  
  # --- Compute inter-block trait distances ---
  disTraits <- dist.ktab(ktabList, groups, scan = FALSE, option = "scaledBYrange")
  
  # --- Run PCoA ---
  PCoATraits <- pcoa(disTraits)
  rownames(PCoATraits$vectors) <- namesx
  
  # --- envfit correlation with external traits ---
  plotenvfit <- vegan::envfit(PCoATraits$vectors, datax, w = 1)
  spp.scrs <- as.data.frame(scores(plotenvfit, display = "vectors"))
  spp.scrs$Species <- rownames(spp.scrs)
  
  # --- Parallel + progress bar for CorrelType loop ---
  PcoaCorrel <- NULL
  with_progress({
    p <- progressor(steps = length(x))
    
    PcoaCorrel <- future_lapply(seq_along(x), function(i) {
      res <- CorrelType(PCoATraits, x[[i]])
      p(sprintf("Trait table %d/%d done", i, length(x)))
      return(res)
    })
  })
  
  PcoaCorrel <- do.call(rbind, PcoaCorrel)
  
  # --- Output ---
  res <- list(
    PCoA = PCoATraits,
    varPlot = spp.scrs,
    PCoACor = PcoaCorrel
  )
  
  if (!is.na(saveFile)) saveRDS(res, saveFile)
  
  return(res)
}


# make.TPD.2D.high.def: Compute the TPDs
make.TPD.2D.high.def=function(traitsUSE,dimensions,alphaUse=0.95,gridSize=100,saveFile= paste0(getwd(),"/TPDs2DHighDef.rds")){
  require(TPD)
  colnames(traitsUSE)=paste0("Comp.",1:dimensions)
  sdTraits = sqrt(diag(Hpi.diag(traitsUSE)))
  if(dimensions > 2){
    TPDsAux = TPDsMean_large(species = rownames(traitsUSE), 
                       means = traitsUSE, 
                       sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                       alpha = alphaUse,
                       n_divisions = gridSize)
  }else{
    TPDsAux = TPDsMean(species = rownames(traitsUSE), 
                       means = traitsUSE, 
                       sds = matrix(rep(sdTraits, nrow(traitsUSE)), byrow=T, ncol=dimensions),
                       alpha = alphaUse,
                       n_divisions = gridSize)
  }
  

  if(!is.na(saveFile)){
    saveRDS(TPDsAux, saveFile)
  }
  return(TPDsAux)
}


###
### TPDsMean_large: function to estimate TPDs functions using a single average value per species and a given bandwidth (standard deviation). This function is equivalent to TPD::TPDsMean, but does not record the cells with zero probability, making it more suitable for higher dimensions.
### 
TPDsMean_large<- function(species, means, sds, alpha = 0.95, samples = NULL,
                          trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {
  
  # INITIAL CHECKS:
  #   1. Compute the number of dimensions (traits):
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
         number of dimensions")
  }
  #   2. sds and means must have the same dimensions:
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  #   3. species and means must have the same "dimensions":
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means'
         and 'sds'")
  }
  #	4. NA's not allowed in means, sds & species:
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  #	5. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means'
         and 'sds'")
    }
    if (any(is.na(samples))) {
      stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  }
  
  #	6. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	6. Create the grid of cells in which the density function is evaluated:
  if (is.null(n_divisions)) {
    n_divisions_choose<- c(1000, 200, 50, 25)
    n_divisions<- n_divisions_choose[dimensions]
  }
  grid_evaluate<-list()
  edge_length <- list()
  cell_volume<-1
  for (dimens in 1:dimensions){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] -
      grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions==1){
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  # Creation of lists to store results:
  results <- list()
  # DATA: To store data and common features
  results$data <- list()
  results$data$evaluation_grid <- evaluation_grid
  results$data$cell_volume <- cell_volume
  results$data$edge_length <- edge_length
  results$data$species <- species
  results$data$means <- means
  results$data$sds <- sds
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }
  
  results$data$alpha <- alpha
  results$data$pop_means <- list()
  results$data$pop_sds <- list()
  results$data$pop_sigma <- list()
  results$data$dimensions <- dimensions
  results$data$type <- type
  results$data$method <- "mean"
  
  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()
  
  
  ########Multivariate normal density calculation
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi]<-
      unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      
      multNormAux <- mvtnorm::dmvnorm(x = evaluation_grid,
                                      mean = results$data$pop_means[[spi]],
                                      sigma = results$data$pop_sigma[[spi]])
      multNormAux <- multNormAux / sum(multNormAux)
      # Now, we extract the selected fraction of volume (alpha), if necessary
      extract_alpha <- function(x){
        # 1. Order the 'cells' according to their probabilities:
        alphaSpace_aux <- x[order(x, decreasing=T)]
        # 2. Find where does the accumulated sum of ordered probabilities becomes
        #   greater than the threshold (alpha):
        greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux ) > alpha) [1]]
        # 3. Substitute smaller probabilities with 0:
        x[x < greater_prob] <- 0
        # 5. Finally, reescale, so that the total density is finally 1:
        x <- x / sum(x)
        return(x)
      }
      if (alpha < 1){
        multNormAux <- extract_alpha(multNormAux)
      }
      notZeroIndex <- which(multNormAux != 0)
      notZeroProb <- multNormAux[notZeroIndex]
      results$TPDs[[spi]] <- cbind(notZeroIndex, notZeroProb)
      
      
    }
    if (dimensions == 1) stop("This function is intended for > 1 dimension")
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}






imageTPD = function(x, thresholdPlot = 0.99){
  TPDList = x$TPDc$TPDc
  imageTPD = list()
  for(comm in 1:length(TPDList)){
    percentile = rep(NA, length(TPDList[[comm]]))
    TPDList[[comm]]= cbind(index = 1:length(TPDList[[comm]]),
                           prob = TPDList[[comm]], percentile)
    orderTPD = order(TPDList[[comm]][,"prob"], decreasing = T)
    TPDList[[comm]] = TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] = cumsum(TPDList[[comm]][,"prob"])
    TPDList[[comm]] = TPDList[[comm]][order(TPDList[[comm]][,"index"]),]
    imageTPD[[comm]] = TPDList[[comm]]
  }
  names(imageTPD) = names(TPDList)
  spacePercentiles = matrix(data=0, nrow = nrow(x$data$evaluation_grid),
                            ncol = length(TPDList),
                            dimnames = list(1:nrow(x$data$evaluation_grid),
                                            names(TPDList)))
  trait1Edges = unique(x$data$evaluation_grid[,1])
  trait2Edges = unique(x$data$evaluation_grid[,2])
  imageMat = array(NA, c(length(trait1Edges),
                         length(trait2Edges),
                         length(imageTPD)),
                   dimnames = list(trait1Edges, trait2Edges, names(TPDList)))
  for(comm in 1:length(TPDList)){
    percentileSpace = x$data$evaluation_grid
    percentileSpace$percentile = rep(NA, nrow(percentileSpace))
    percentileSpace[, "percentile"] =imageTPD[[comm]][,"percentile"]
    for(i in 1:length(trait2Edges)){
      colAux = subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])
      imageMat[, i, comm] = colAux$percentile
    }
    imageMat[imageMat > thresholdPlot] = NA
  }
  return(imageMat)
}

make.functional.space=function(species,traitsUSE,TPDsAux,PCAImpute,imageTPD,
                               ncolors,limX,limY,ColorRamp,nonat,extirp,title=NA,cex=2,alpha=0.8,colPointsI="#8B3E2F",colPointsE="#6495ED"){
  
  par(mar=c(5,5,4,4))   
  par(mgp=c(2.2,0.1,0))   
  par(cex.axis = .9, cex.lab = 1)
  dataAnalysis=traitsUSE[rownames(traitsUSE)%in%species,c(1,2)]
  var=NULL
  var=PCAImpute$variance
  if(is.null(var)){var=PCAImpute$Variance}
  commNames = c("All")
  commMatrix = matrix(0, nrow = length(commNames), ncol = nrow(dataAnalysis),
                      dimnames = list(commNames, rownames(dataAnalysis)))
  commMatrix["All", ] = 1
  commMatrix=commMatrix[,colnames(commMatrix)%in%TPDsAux$data$species,drop = F]
  
  TPDcAux = TPDc(TPDs = TPDsAux, sampUnit = commMatrix)
  xlab = paste0("PC1 (", round(100 * var[1], 2), "%)") 
  ylab = paste0("PC2 (", round(100 * var[2], 2), "%)")
  imageMat =imageTPD(TPDcAux, thresholdPlot = 1)
  trait1Edges = unique(TPDcAux$data$evaluation_grid[,1])
  trait2Edges = unique(TPDcAux$data$evaluation_grid[,2])
  
  dataPoints=cbind.data.frame(traitsUSE[nonat,],z=rep(1,nrow(traitsUSE[nonat,])))
  dataPoints$Names=unlist(lapply(strsplit(rownames(dataPoints),"_"),function(x){paste0(substr(x[1],1,1),".",substr(x[2],1,5))}))
  dataPoints$AllNames=rownames(dataPoints)
  dataPointsNonNat=dataPoints
  
  dataPoints=cbind.data.frame(traitsUSE[extirp,],z=rep(1,nrow(traitsUSE[extirp,])))
  dataPoints$Names=unlist(lapply(strsplit(rownames(dataPoints),"_"),function(x){paste0(substr(x[1],1,1),".",substr(x[2],1,5))}))
  dataPoints$AllNames=rownames(dataPoints)
  dataPointsExtirp=dataPoints
  
  
  data=reshape2::melt(imageMat[, , "All"])
  colnames(data)=c("x","y","z")
  plot=ggplot(data, aes(x, y, z=z))+ 
    geom_raster(aes(fill = z)) +  
    geom_contour(colour="grey75",breaks = seq(0.5, 0.99, by=0.4))+
    geom_contour(colour="grey50",breaks = 0.995,lty=2)+
    scale_fill_gradientn(na.value = "white",colours="white")+
    geom_point(data=dataPointsExtirp,aes(x=Comp.1,y=Comp.2),size=cex,alpha=alpha-0.2,shape=1,colour=colPointsE,fill=colPointsE)+
    geom_point(data=dataPointsNonNat,aes(x=Comp.1,y=Comp.2),size=cex,alpha=alpha-0.1,shape=21,colour=colPointsI,fill=colPointsI)+
    xlab(xlab)+
    ylab(ylab)+
    theme(axis.line = element_blank(), 
          axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(colour = "black", size = 12), 
          legend.position = "none", 
          plot.title = element_text(size = 10,hjust=0), 
          plot.background = element_blank(), 
          panel.background = element_rect(fill="white",colour="grey50"), 
          panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = ggplot2::margin(.8,.8,.5,.5, "cm"),
          axis.title.y = element_text(size=14,margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(size=14,margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)))
  plot
}

PCA_plot = function(PCoAPlot,PCoACorPlot,multAx1,multAx,title,legend){
  Plan12=ggplot(data = data.frame(PCoAPlot$vectors)) +
    geom_point(aes(x = Axis.1, y = Axis.2),col="grey89") +
    theme_classic()+
    ggtitle(title)+
    xlab(paste0("PCoA 1 (",round(PCoAPlot$values[1,2]*100,2),"%)"))+
    ylab(paste0("PCoA 2 (",round(PCoAPlot$values[2,2]*100,2),"%)"))+
    #coord_fixed() + ## need aspect ratio of 1!
    geom_segment(data = data.frame(PCoACorPlot),
                 aes(x = 0, xend = Axis.1*multAx, y = 0, yend = Axis.2*multAx),
                 arrow = arrow(length = unit(0.25, "cm")), colour = PCoACorPlot$color) +
    geom_text(data = data.frame(PCoACorPlot), 
              aes(x =  Axis.1*multAx1, y = Axis.2*multAx1, 
                  label =  names), colour = PCoACorPlot$color,
              size = 5) +
    theme(
      plot.title = element_text(
        size = 18,       # Increase size (adjust as needed)
        face = "bold"    # Make it bold
      )
    )
  
  if(!is.null(legend)){
    Plan12 = Plan12 + annotate(geom = "text",label = legend, x = Inf, y = -Inf, hjust = 1, vjust = -0.5)      # Make text bold)
  }

  
   
  Plan12
}

PCA_plot_34 = function(PCoAPlot,PCoACorPlot,multAx1,multAx,title,legend){
  Plan12=ggplot(data = data.frame(PCoAPlot$vectors)) +
    geom_point(aes(x = Axis.3, y = Axis.4),col="grey89") +
    theme_classic()+
    ggtitle(title)+
    xlab(paste0("PCoA 3 (",round(PCoAPlot$values[3,2]*100,2),"%)"))+
    ylab(paste0("PCoA 4 (",round(PCoAPlot$values[4,2]*100,2),"%)"))+
    #coord_fixed() + ## need aspect ratio of 1!
    geom_segment(data = data.frame(PCoACorPlot),
                 aes(x = 0, xend = Axis.3*multAx, y = 0, yend = Axis.4*multAx),
                 arrow = arrow(length = unit(0.25, "cm")), colour = PCoACorPlot$color) +
    geom_text(data = data.frame(PCoACorPlot), 
              aes(x =  Axis.3*multAx1, y = Axis.4*multAx1, 
                  label =  names), colour = PCoACorPlot$color,
              size = 5) +
    theme(
      plot.title = element_text(
        size = 18,       # Increase size (adjust as needed)
        face = "bold"    # Make it bold
      )
    )
  
  if(!is.null(legend)){
    Plan12 = Plan12 + annotate(geom = "text",label = legend, x = Inf, y = -Inf, hjust = 1, vjust = -0.5)      # Make text bold)
  }
  
  
  
  Plan12
}

CorrelType=function(pcoa,type){
  
  Axis=matrix(NA,nc=ncol(pcoa$vectors),nr=ncol(type),dimnames = list(colnames(type),colnames(pcoa$vectors)))
  for (i in 1:ncol(type)){
    for (tr in 1:ncol(pcoa$vectors)){
      Axis[i,tr]=cor(pcoa$vectors[,tr],type[,i])
    }
  }
  return(Axis)
}

TPDRichness = function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results = list()
  # 1. Functional Richness
  Calc_FRich = function(x) {
    results_FR = numeric()
    if (class(x) == "TPDcomm") {
      TPD = x$TPDc$TPDc
      names_aux = names(x$TPDc$TPDc)
      cell_volume = x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD = x$TPDs
      names_aux = names(x$TPDs)
      cell_volume = x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux = TPD[[i]]
      TPD_aux[TPD_aux > 0] = cell_volume
      results_FR[i] = sum(TPD_aux)
    }
    names(results_FR) = names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities = list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness = Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species = list()
      # message("Calculating FRichness of species")
      results$species$FRichness = Calc_FRich(TPDs)
    } else {
      results$populations = list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness = Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}


### quantileTPD_large: function to transform probabilities from TPDc_large object into quantiles
### x is a TPDc_large object 
quantileTPD_large <- function(x){
  TPDList <- x$TPDc$TPDc
  results <- list()
  for(comm in 1:length(TPDList)){
    percentile <- rep(NA, nrow(TPDList[[comm]]))
    TPDList[[comm]]<- cbind(TPDList[[comm]], percentile)
    orderTPD <- order(TPDList[[comm]][,"notZeroProb"], decreasing = T)
    TPDList[[comm]] <- TPDList[[comm]][orderTPD,]
    TPDList[[comm]][,"percentile"] <- 1- cumsum(TPDList[[comm]][,"notZeroProb"])
    TPDList[[comm]] <- TPDList[[comm]][order(TPDList[[comm]][,"notZeroIndex"]),]
    results[[comm]] <- TPDList[[comm]]
  }
  names(results) <- names(TPDList)
  return(results)
}

### MAE_TPD_large: function to estimate differences of quantiles between two quantileTPD_large objects
MAE_TPD_large <- function(x, y){# x and y are made with quantileTPD_large function
  mergeXY <- merge(x, y, by="notZeroIndex", all=T)
  mergeXY[is.na(mergeXY)] <- 0
  mergeXY$change <- abs(mergeXY$percentile.x - mergeXY$percentile.y)
  MAE <- mean(mergeXY$change)
  return(MAE)
}

### TPDsMean_large: function to estimate TPDs functions using a single average value per species and a given bandwidth (standard deviation). This function is equivalent to TPD::TPDsMean, but does not record the cells with zero probability, making it more suitable for higher dimensions.
TPDsMean_large<- function(species, means, sds, alpha = 0.95, samples = NULL,
                          trait_ranges = NULL, n_divisions = NULL, tolerance = 0.05) {
  
  # INITIAL CHECKS:
  #   1. Compute the number of dimensions (traits):
  means <- as.matrix(means)
  dimensions <- ncol(means)
  if (dimensions > 4) {
    stop("No more than 4 dimensions are supported at this time; reduce the
         number of dimensions")
  }
  #   2. sds and means must have the same dimensions:
  sds <- as.matrix(sds)
  if (all(dim(means) != dim(sds))) {
    stop("'means' and 'sds' must have the same dimensions")
  }
  #   3. species and means must have the same "dimensions":
  if (length(species) != nrow(means)) {
    stop("The length of 'species' does not match the number of rows of 'means'
         and 'sds'")
  }
  #	4. NA's not allowed in means, sds & species:
  if (any(is.na(means)) | any(is.na(sds)) | any(is.na(species))) {
    stop("NA values are not allowed in 'means', 'sds' or 'species'")
  }
  #	5. Compute the species or populations upon which calculations will be done:
  if (is.null(samples)) {
    species_base <- species
    if (length(unique(species_base)) == 1){
      type <- "One population_One species"
    } else{
      type <- "One population_Multiple species"
    }
  } else {
    if (length(samples) != nrow(means)) {
      stop("The length of 'samples' does not match the number of rows of 'means'
         and 'sds'")
    }
    if (any(is.na(samples))) {
      stop("NA values are not allowed in 'samples'")
    }
    species_base <- paste(species, samples, sep = ".")
    if (length(unique(species)) == 1){
      type <- "Multiple populations_One species"
    } else{
      type <- "Multiple populations_Multiple species"
    }
  }
  
  #	6. Define trait ranges:
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (5, dimensions)
  }
  if (class(trait_ranges) != "list") {
    trait_ranges_aux <- trait_ranges
    trait_ranges <- list()
    for (dimens in 1:dimensions) {
      max_aux <- max(means[, dimens] + trait_ranges_aux[dimens] * sds[, dimens])
      min_aux <- min(means[, dimens] - trait_ranges_aux[dimens] * sds[, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  #	6. Create the grid of cells in which the density function is evaluated:
  if (is.null(n_divisions)) {
    n_divisions_choose<- c(1000, 200, 50, 25)
    n_divisions<- n_divisions_choose[dimensions]
  }
  grid_evaluate<-list()
  edge_length <- list()
  cell_volume<-1
  for (dimens in 1:dimensions){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
    edge_length[[dimens]] <- grid_evaluate[[dimens]][2] -
      grid_evaluate[[dimens]][1]
    cell_volume <- cell_volume * edge_length[[dimens]]
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  if (is.null(colnames(means))){
    names(evaluation_grid) <- paste0("Trait.",1:dimensions)
  } else {
    names(evaluation_grid) <- colnames(means)
  }
  if (dimensions==1){
    evaluation_grid <- as.matrix(evaluation_grid)
  }
  # Creation of lists to store results:
  results <- list()
  # DATA: To store data and common features
  results$data <- list()
  results$data$evaluation_grid <- evaluation_grid
  results$data$cell_volume <- cell_volume
  results$data$edge_length <- edge_length
  results$data$species <- species
  results$data$means <- means
  results$data$sds <- sds
  if (is.null(samples)){
    results$data$populations <-  NA
  } else{
    results$data$populations <-  species_base
  }
  
  results$data$alpha <- alpha
  results$data$pop_means <- list()
  results$data$pop_sds <- list()
  results$data$pop_sigma <- list()
  results$data$dimensions <- dimensions
  results$data$type <- type
  results$data$method <- "mean"
  
  # TPDs: To store TPDs features of each species/population
  results$TPDs<-list()
  
  
  ########Multivariate normal density calculation
  for (spi in 1:length(unique(species_base))) {
    # Some information messages
    if (spi == 1) { message(paste0("-------Calculating densities for ", type, "-----------\n")) }
    #Data selection
    selected_rows <- which(species_base == unique(species_base)[spi])
    results$data$pop_means[[spi]] <- means[selected_rows, ]
    results$data$pop_sds[[spi]] <- sds[selected_rows, ]
    names(results$data$pop_means)[spi] <- names(results$data$pop_sds)[spi]<-
      unique(species_base)[spi]
    if (dimensions > 1) {
      results$data$pop_sigma[[spi]] <- diag(results$data$pop_sds[[spi]]^2)
      
      multNormAux <- mvtnorm::dmvnorm(x = evaluation_grid,
                                      mean = results$data$pop_means[[spi]],
                                      sigma = results$data$pop_sigma[[spi]])
      multNormAux <- multNormAux / sum(multNormAux)
      # Now, we extract the selected fraction of volume (alpha), if necessary
      extract_alpha <- function(x){
        # 1. Order the 'cells' according to their probabilities:
        alphaSpace_aux <- x[order(x, decreasing=T)]
        # 2. Find where does the accumulated sum of ordered probabilities becomes
        #   greater than the threshold (alpha):
        greater_prob <- alphaSpace_aux[which(cumsum(alphaSpace_aux ) > alpha) [1]]
        # 3. Substitute smaller probabilities with 0:
        x[x < greater_prob] <- 0
        # 5. Finally, reescale, so that the total density is finally 1:
        x <- x / sum(x)
        return(x)
      }
      if (alpha < 1){
        multNormAux <- extract_alpha(multNormAux)
      }
      notZeroIndex <- which(multNormAux != 0)
      notZeroProb <- multNormAux[notZeroIndex]
      results$TPDs[[spi]] <- cbind(notZeroIndex, notZeroProb)
      
      
    }
    if (dimensions == 1) stop("This function is intended for > 1 dimension")
  }
  names(results$TPDs) <- unique(species_base)
  class(results) <- "TPDsp"
  return(results)
}

### TPDc_large: function to estimate Trait Probability Density of Communities based on TPDs_large objects
TPDc_large <- function(TPDs, sampUnit){
  sampUnit <- as.matrix(sampUnit)
  # 1. species names:
  if (is.null(colnames(sampUnit)) | any(is.na(colnames(sampUnit)))) {
    stop("colnames(sampUnit), must contain the names of the species.
      NA values are not allowed")
  }
  # 2. communities names:
  if (is.null(rownames(sampUnit)) | any(is.na(rownames(sampUnit)))){
    stop("rownames(sampUnit), must contain the names of the sampling units.
      NA values are not allowed")
  }
  # 3. Data values will be inherithed from TPDs, which must be of class TPD
  if (class(TPDs) != "TPDsp") {
    stop("TPDs must be an object of class 'TPDsp', created with the function
      'TPDs'")
  }
  species <- samples <- abundances <- numeric()
  for (i in 1:nrow(sampUnit)){
    samples <- c(samples, rep(rownames(sampUnit)[i], ncol(sampUnit)))
    species <- c(species, colnames(sampUnit))
    abundances <- c(abundances, sampUnit[i, ])
  }
  nonZero <- which(abundances > 0)
  samples <- samples[nonZero]
  species <- species[nonZero]
  abundances <- abundances[nonZero]
  
  # Creation of lists to store results:
  results <- list()
  results$data <- TPDs$data
  results$data$sampUnit <- sampUnit
  type <- results$data$type
  # All the species or populations in 'species' must be in the species or
  #   populations of TPDs:
  if (type == "Multiple populations_One species" |
      type == "Multiple populations_Multiple species") {
    species_base <- paste(species, samples, sep = ".")
    if (!all(unique(species_base) %in% unique(results$data$populations))) {
      non_found_pops <- which(unique(species_base) %in%
                                unique(results$data$populations) == 0)
      stop("All the populations TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_pops], collapse=" / "))
    }
  }
  if (type == "One population_One species" |
      type == "One population_Multiple species") {
    species_base <- species
    if (!all(unique(species_base) %in% unique(results$data$species))) {
      non_found_sps <- which(unique(species_base) %in%
                               unique(results$data$species) == 0)
      stop("All the species TPDs must be present in 'TPDs'. Not present:\n",
           paste(species_base[non_found_sps], collapse=" / "))
    }
  }
  # END OF INITIAL CHECKS
  # TPDc computation
  results$TPDc <- list()
  results$TPDc$species <- list()
  results$TPDc$abundances <- list()
  results$TPDc$speciesPerCell <- list()
  # results$TPDc$RTPDs <- list()
  results$TPDc$TPDc <- list()
  
  for (samp in 1:length(unique(samples))) {
    selected_rows <- which(samples == unique(samples)[samp])
    species_aux <- species_base[selected_rows]
    abundances_aux <- abundances[selected_rows] / sum(abundances[selected_rows])
    RTPDsAux <- rep(0, nrow(results$data$evaluation_grid))
    TPDs_aux <- TPDs$TPDs[names(TPDs$TPDs) %in% species_aux]
    cellsOcc <- numeric()
    for (sp in 1:length(TPDs_aux)) {
      selected_name <- which(names(TPDs_aux) == species_aux[sp])
      cellsToFill <- TPDs_aux[[selected_name]][,"notZeroIndex"]
      cellsOcc <- c(cellsOcc, cellsToFill)
      probsToFill <- TPDs_aux[[selected_name]][,"notZeroProb"] * abundances_aux[sp]
      RTPDsAux[cellsToFill] <- RTPDsAux[cellsToFill] + probsToFill
    }
    TPDc_aux <- RTPDsAux
    notZeroIndex <- which(TPDc_aux != 0)
    notZeroProb <- TPDc_aux[notZeroIndex]
    results$TPDc$TPDc[[samp]] <- cbind(notZeroIndex, notZeroProb)
    results$TPDc$species[[samp]] <- species_aux
    results$TPDc$abundances[[samp]] <- abundances_aux
    results$TPDc$speciesPerCell[[samp]] <- table(cellsOcc)
    names(results$TPDc$TPDc)[samp] <-
      names(results$TPDc$species)[samp] <- names(results$TPDc$abundances)[samp] <-
      names(results$TPDc$speciesPerCell)[samp] <- unique(samples)[samp]
  }
  class(results) <- "TPDcomm"
  return(results)
}

### TPDRichness_large: Function to estimate functional richness from TPD_large object
TPDRichness_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    # message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      # message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
    } else {
      results$populations <- list()
      # message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### REND_large: Functional Evenness, Richness and Divergence from TPD_large object
REND_large <- function(TPDc = NULL, TPDs = NULL){
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp"){
    stop("The class of one object do not match the expectations,
         Please, specify if your object is a TPDc or a TPDs")
  }
  # Creation of lists to store results:
  results <- list()
  # 1. Functional Richness
  Calc_FRich <- function(x) {
    results_FR <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_aux[TPD_aux > 0] <- cell_volume
      results_FR[i] <- sum(TPD_aux)
    }
    names(results_FR) <- names_aux
    return(results_FR)
  }
  
  # 2. Functional Evenness
  Calc_FEve <- function(x) {
    results_FE <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      TPD_eve <- rep((1 / length(TPD_aux)), times = length(TPD_aux))
      results_FE[i] <- sum(pmin(TPD_aux, TPD_eve))
    }
    names(results_FE) <- names_aux
    return(results_FE)
  }
  # 3. Functional Divergence
  Calc_FDiv <- function(x) {
    results_FD <- numeric()
    if (class(x) == "TPDcomm") {
      TPD <- x$TPDc$TPDc
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDc$TPDc)
      cell_volume <- x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD <- x$TPDs
      evaluation_grid<-x$data$evaluation_grid
      names_aux <- names(x$TPDs)
      cell_volume <- x$data$cell_volume
    }
    for (i in 1:length(TPD)) { 
      notZeroCells <- TPD[[i]][,"notZeroIndex"]
      functional_volume <- evaluation_grid[notZeroCells, , drop=F]
      # Functional volume has to be standardised so that distances are
      # independent of the scale of the axes
      for (j in 1:ncol(functional_volume)){
        functional_volume[, j] <-
          (functional_volume[, j] - min(functional_volume[, j])) /
          (max(functional_volume[, j]) - min(functional_volume[, j]))
      }
      TPD_aux <- TPD[[i]][,"notZeroProb"]
      # 1. Calculate the center of gravity
      COG <- colMeans(functional_volume, na.rm=T)
      # 2. Calculate the distance of each point in the functional volume to the
      #   COG:
      dist_COG <- function(x, COG) {
        result_aux<-stats::dist(rbind(x, COG))
        return(result_aux)
      }
      COGDist <- apply(functional_volume, 1, dist_COG, COG)
      # 3. Calculate the mean of the COGDist's
      meanCOGDist <- mean(COGDist)
      # 4. Calculate the sum of the abundance-weighted deviances for distaces
      #   from the COG (AWdistDeviances) and the absolute abundance-weighted
      #   deviances:
      distDeviances <- COGDist - meanCOGDist
      AWdistDeviances <- sum(TPD_aux * distDeviances)
      absdistDeviances <- abs(COGDist - meanCOGDist)
      AWabsdistDeviances <- sum(TPD_aux * absdistDeviances)
      #Finally, calculate FDiv:
      results_FD[i] <- (AWdistDeviances + meanCOGDist) /
        ( AWabsdistDeviances +  meanCOGDist)
    }
    names(results_FD) <- names_aux
    return(results_FD)
  }
  # IMPLEMNENTATION
  if (!is.null(TPDc)) {
    results$communities <- list()
    message("Calculating FRichness of communities")
    results$communities$FRichness <- Calc_FRich(TPDc)
    message("Calculating FEvenness of communities")
    results$communities$FEvenness <- Calc_FEve(TPDc)
    message("Calculating FDivergence of communities")
    results$communities$FDivergence <- Calc_FDiv(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" |
        TPDs$data$type == "One population_Multiple species") {
      results$species <- list()
      message("Calculating FRichness of species")
      results$species$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of species")
      results$species$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of species")
      results$species$FDivergence <- Calc_FDiv(TPDs)
    } else {
      results$populations <- list()
      message("Calculating FRichness of populations")
      results$populations$FRichness <- Calc_FRich(TPDs)
      message("Calculating FEvenness of populations")
      results$populations$FEvenness <- Calc_FEve(TPDs)
      message("Calculating FDivergence of populations")
      results$populations$FDivergence <- Calc_FDiv(TPDs)
    }
    if (TPDs$data$method == "mean") {
      message("WARNING: When TPDs are calculated using the TPDsMean function, Evenness
              and Divergence are meaningless!!")
    }
  }
  return(results)
}

### redundancy_large:  Functional Redundancy of Communities from TPD_large object
redundancy_large <- function(TPDc = NULL) {
  if (class(TPDc) != "TPDcomm") {
    stop("TPDc must be an object of class TPDcomm generated with the TPDc
		    function")
  }
  x <- TPDc
  results <- list()
  results$redundancy <- results$richness <- numeric()
  for (i in 1:length(x$TPDc$TPDc)) {
    TPDc_aux <- x$TPDc$TPDc[[i]][, "notZeroProb"]
    M <- x$TPDc$speciesPerCell[[i]]
    results$redundancy[i] <- sum(M * TPDc_aux) - 1
    results$richness[i] <- sum(x$TPDc$abundances[[i]] >0)
  }
  results$redundancyRelative <- results$redundancy / (results$richness -1)
  names(results$redundancy) <- names(results$richness) <-
    names(results$redundancyRelative) <- names(x$TPDc$TPDc)
  return(results)
}

dissim_large <- function(x = NULL) {
  # INITIAL CHECKS:
  # 1. At least one of TPDc or TPDs must be supplied.
  if (class(x) == "TPDcomm"){
    TPDType<-"Communities"
    TPDc <- x
  } else{
    if (class(x) == "TPDsp"){
      TPDType<-"Populations"
      TPDs <- x
    } else{
      stop("x must be an object of class TPDcomm or TPDsp")
    }
  }
  results <- list()
  Calc_dissim <- function(x) {
    # 1. BetaO (functional dissimilarity)
    results_samp <- list()
    if (TPDType == "Communities") {
      TPD <- x$TPDc$TPDc
      names_aux <- names(x$TPDc$TPDc)
    }
    if (TPDType == "Populations") {
      TPD <- x$TPDs
      names_aux <- names(x$TPDs)
    }
    results_samp$dissimilarity <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                         dimnames = list(names_aux, names_aux))
    results_samp$P_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                    dimnames = list(names_aux, names_aux))
    results_samp$P_non_shared <- matrix(NA,ncol= length(TPD), nrow= length(TPD),
                                        dimnames = list(names_aux, names_aux))
    for (i in 1:length(TPD)) {
      TPD_i <- TPD[[i]]
      for (j in 1:length(TPD)) {
        if (i > j) {
          TPD_j <- TPD[[j]]
          commonTPD <- rbind(TPD_i, TPD_j)
          #Now, select only those cells that appear in both i and j
          duplicatedCells <- names(which(table(commonTPD[,"notZeroIndex"])==2))
          doubleTPD <- commonTPD[which(commonTPD[,"notZeroIndex"] %in% duplicatedCells), ]
          O_aux <- sum(tapply(doubleTPD[,"notZeroProb"], doubleTPD[,"notZeroIndex"], min))
          A_aux <- sum(tapply(doubleTPD[,"notZeroProb"], doubleTPD[,"notZeroIndex"], max)) - O_aux
          only_in_i_aux <- which(TPD_i[,"notZeroIndex"] %in% setdiff(TPD_i[,"notZeroIndex"], 
                                                                     TPD_j[,"notZeroIndex"]))	
          B_aux <- sum(TPD_i[only_in_i_aux,"notZeroProb"])
          only_in_j_aux <- which(TPD_j[,"notZeroIndex"] %in% setdiff(TPD_j[,"notZeroIndex"], 
                                                                     TPD_i[,"notZeroIndex"]))	
          C_aux <- sum(TPD_j[only_in_j_aux,"notZeroProb"])
          results_samp$dissimilarity[i, j] <-
            results_samp$dissimilarity[j, i] <-	1 - O_aux
          if (results_samp$dissimilarity[j, i] == 0) {
            results_samp$P_non_shared[i, j] <- NA
            results_samp$P_non_shared[j, i] <- NA
            results_samp$P_shared[i, j] <- NA
            results_samp$P_shared[j, i] <- NA
          }	else {
            results_samp$P_non_shared[i, j] <-
              results_samp$P_non_shared[j, i] <-
              (2 * min(B_aux, C_aux)) / (A_aux + 2 * min(B_aux, C_aux))
            results_samp$P_shared[i, j] <- results_samp$P_shared[j, i] <-
              1 - results_samp$P_non_shared[i, j]
          }
        }
        if (i == j) {
          results_samp$dissimilarity[i, j] <- 0
        }
      }
    }
    return(results_samp)
  }
  # IMPLEMENTATION
  if (TPDType == "Communities") {
    message("Calculating dissimilarities between ", length(TPDc$TPDc$TPDc)," communities. It might take some time")
    results$communities <- Calc_dissim(TPDc)
  }
  if (TPDType == "Populations") {
    message("Calculating dissimilarities between ", length(TPDs$TPDs)," populations. It might take some time")
    results$populations <- Calc_dissim(TPDs)
  }
  class(results) <- "OverlapDiss"
  return(results)
}


normalNull = function(traitsUSE,TPDsObs,dimensions = 2,TPDsMean,TPDc,REND,densityProfileTPD,
                      totalNull,nspsSel){
  RawResultsNormal<- RawResultsObserved  <- list()
  RawResultsNormal[["FRic100"]] <- rep(NA, totalNull)
  RawResultsNormal[["FRic99"]] <- rep(NA, totalNull)
  RawResultsNormal[["FRic50"]] <- rep(NA, totalNull)
  RawResultsNormal[["FEve"]] <- rep(NA, totalNull)
  RawResultsNormal[["FDiv"]] <- rep(NA, totalNull)
  RawResultsObserved[["FRic99"]] <- rep(NA, totalNull)
  RawResultsObserved[["FRic50"]] <- rep(NA, totalNull)
  RawResultsObserved[["FEve"]] <- rep(NA, totalNull)
  RawResultsObserved[["FDiv"]] <- rep(NA, totalNull)
  
  
  meansAux <- colMeans(traitsUSE)
  covAux <- cov(traitsUSE)
  normalTraits <- mvtnorm::rmvnorm(n = nrow(traitsUSE), mean = meansAux, sigma=covAux)
  normalSD <- sqrt(diag(Hpi.diag(normalTraits)))
  TPDsNormal <- TPDsMean(species = names(TPDsObs$TPDs), 
                         means =  normalTraits,
                         sds = matrix(rep(normalSD, nrow(normalTraits)), 
                                      byrow=T, ncol=ncol(normalTraits)),
                         alpha = TPDsObs$data$alpha,
                         n_divisions = nrow(TPDsObs$data$evaluation_grid)^(1/ncol(normalTraits)))
  
  ### Multivariate normal
  for(i in 1:totalNull){
    cat(paste("\r Rep: ", i, "out of", totalNull, "\r"))
    spSample <- sample(1:nrow(traitsUSE), nspsSel)
    commNames <- c("Sample")
    commMatrixAux <- matrix(0, nrow = length(commNames), ncol = nrow(normalTraits),
                            dimnames = list(commNames, names(TPDsObs$TPDs)))
    commMatrixAux[1, spSample] <- 1
    TPDcObsAux <- TPDc(TPDs = TPDsObs, sampUnit = commMatrixAux)
    RENDObsAux <- REND(TPDc=TPDcObsAux)
    densProfObsAux <- densityProfileTPD(TPDcObsAux, probs=c(0.5, 0.99))
    RawResultsObserved[["FRic99"]][i] <- densProfObsAux[2]
    RawResultsObserved[["FRic50"]][i] <- densProfObsAux[1]
    RawResultsObserved[["FEve"]][i] <- RENDObsAux$communities$FEvenness
    RawResultsObserved[["FDiv"]][i] <- RENDObsAux$communities$FDivergence
    TPDcObsAux <- NULL
    
    TPDcNormalAux <- TPDc(TPDs = TPDsNormal, sampUnit = commMatrixAux)
    RENDNormalAux <- REND(TPDc=TPDcNormalAux)
    densProfNormalAux <- densityProfileTPD(TPDcNormalAux, probs=c(0.5, 0.99,1))
    RawResultsNormal[["FRic100"]][i] <- densProfNormalAux[3]
    RawResultsNormal[["FRic99"]][i] <- densProfNormalAux[2]
    RawResultsNormal[["FRic50"]][i] <- densProfNormalAux[1]
    RawResultsNormal[["FEve"]][i] <- RENDNormalAux$communities$FEvenness
    RawResultsNormal[["FDiv"]][i] <- RENDNormalAux$communities$FDivergence
    TPDcNormalAux <- NULL
  }  
  return(list(norm = RawResultsNormal, obser = RawResultsObserved))
}

### 3. Calculate functional indices, adapted from 'TPD' package ----

Calc_uniq = function(x, y) {
  TPDs_aux = x$TPDs
  TPDc_aux = y$TPDc$TPDc
  names_aux_com = names(y$TPDc$TPDc)
  names_aux_sp = names(x$TPDs)
  results = matrix(NA, ncol = length(TPDs_aux), nrow = length(TPDc_aux), 
                   dimnames = list(names_aux_com, names_aux_sp))
  for (i in 1:length(TPDc_aux)) {
    TPD_ci = TPDc_aux[[i]]
    for (j in 1:length(TPDs_aux)) {
      TPD_sj = TPDs_aux[[j]]
      O_aux = sum(pmin(TPD_ci, TPD_sj))
      results[i, j] = 1 - O_aux
    }
  }
  return(results)
}
REND = function (TPDc = NULL, TPDs = NULL) {
  if (is.null(TPDc) & is.null(TPDs)) {
    stop("At least one of 'TPDc' or 'TPDs' must be supplied")
  }
  if (!is.null(TPDc) & class(TPDc) != "TPDcomm") {
    stop("The class of one object do not match the expectations,\n         Please, specify if your object is a TPDc or a TPDs")
  }
  if (!is.null(TPDs) & class(TPDs) != "TPDsp") {
    stop("The class of one object do not match the expectations,\n         Please, specify if your object is a TPDc or a TPDs")
  }
  results = list()
  Calc_FRich = function(x) {
    results_FR = numeric()
    if (class(x) == "TPDcomm") {
      TPD = x$TPDc$TPDc
      names_aux = names(x$TPDc$TPDc)
      cell_volume = x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD = x$TPDs
      names_aux = names(x$TPDs)
      cell_volume = x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux = TPD[[i]]
      TPD_aux[TPD_aux > 0] = cell_volume
      results_FR[i] = sum(TPD_aux)
    }
    names(results_FR) = names_aux
    return(results_FR)
  }
  Calc_FEve = function(x) {
    results_FE = numeric()
    if (class(x) == "TPDcomm") {
      TPD = x$TPDc$TPDc
      names_aux = names(x$TPDc$TPDc)
      cell_volume = x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD = x$TPDs
      names_aux = names(x$TPDs)
      cell_volume = x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      TPD_aux = TPD[[i]][TPD[[i]] > 0]
      TPD_eve = rep((1/length(TPD_aux)), times = length(TPD_aux))
      results_FE[i] = sum(pmin(TPD_aux, TPD_eve))
    }
    names(results_FE) = names_aux
    return(results_FE)
  }
  Calc_FDiv = function(x) {
    results_FD = numeric()
    if (class(x) == "TPDcomm") {
      TPD = x$TPDc$TPDc
      evaluation_grid = x$data$evaluation_grid
      names_aux = names(x$TPDc$TPDc)
      cell_volume = x$data$cell_volume
    }
    if (class(x) == "TPDsp") {
      TPD = x$TPDs
      evaluation_grid = x$data$evaluation_grid
      names_aux = names(x$TPDs)
      cell_volume = x$data$cell_volume
    }
    for (i in 1:length(TPD)) {
      functional_volume = evaluation_grid[TPD[[i]] > 0, 
                                          , drop = F]
      for (j in 1:ncol(functional_volume)) {
        functional_volume[, j] = (functional_volume[, 
                                                    j] - min(functional_volume[, j]))/(max(functional_volume[, 
                                                                                                             j]) - min(functional_volume[, j]))
      }
      TPD_aux = TPD[[i]][TPD[[i]] > 0]
      COG = colMeans(functional_volume, na.rm = T)
      dist_COG = function(x, COG) {
        result_aux = stats::dist(rbind(x, COG))
        return(result_aux)
      }
      COGDist = apply(functional_volume, 1, dist_COG, 
                      COG)
      meanCOGDist = mean(COGDist)
      distDeviances = COGDist - meanCOGDist
      AWdistDeviances = sum(TPD_aux * distDeviances)
      absdistDeviances = abs(COGDist - meanCOGDist)
      AWabsdistDeviances = sum(TPD_aux * absdistDeviances)
      results_FD[i] = (AWdistDeviances + meanCOGDist)/(AWabsdistDeviances + 
                                                         meanCOGDist)
    }
    names(results_FD) = names_aux
    return(results_FD)
  }
  if (!is.null(TPDc)) {
    results$communities = list()
    results$communities$FRichness = Calc_FRich(TPDc)
    results$communities$FEvenness = Calc_FEve(TPDc)
    results$communities$FDivergence = Calc_FDiv(TPDc)
  }
  if (!is.null(TPDs)) {
    if (TPDs$data$type == "One population_One species" | 
        TPDs$data$type == "One population_Multiple species") {
      results$species = list()
      results$species$FRichness = Calc_FRich(TPDs)
      results$species$FEvenness = Calc_FEve(TPDs)
      results$species$FDivergence = Calc_FDiv(TPDs)
    }
    else {
      results$populations = list()
      results$populations$FRichness = Calc_FRich(TPDs)
      results$populations$FEvenness = Calc_FEve(TPDs)
      results$populations$FDivergence = Calc_FDiv(TPDs)
    }
  }
  return(results)
}
dissim = function (x = NULL) {
  if (class(x) == "TPDcomm") {
    TPDType = "Communities"
    TPDc = x
  }
  else {
    if (class(x) == "TPDsp") {
      TPDType = "Populations"
      TPDs = x
    }
    else {
      stop("x must be an object of class TPDcomm or TPDsp")
    }
  }
  results = list()
  Calc_dissim = function(x) {
    results_samp = list()
    if (TPDType == "Communities") {
      TPD = x$TPDc$TPDc
      names_aux = names(x$TPDc$TPDc)
    }
    if (TPDType == "Populations") {
      TPD = x$TPDs
      names_aux = names(x$TPDs)
    }
    results_samp$dissimilarity = matrix(NA, ncol = length(TPD), 
                                        nrow = length(TPD), dimnames = list(names_aux, names_aux))
    results_samp$P_shared = matrix(NA, ncol = length(TPD), 
                                   nrow = length(TPD), dimnames = list(names_aux, names_aux))
    results_samp$P_non_shared = matrix(NA, ncol = length(TPD), 
                                       nrow = length(TPD), dimnames = list(names_aux, names_aux))
    for (i in 1:length(TPD)) {
      TPD_i = TPD[[i]]
      for (j in 1:length(TPD)) {
        if (i > j) {
          TPD_j = TPD[[j]]
          O_aux = sum(pmin(TPD_i, TPD_j))
          shared_aux = which(TPD_i > 0 & TPD_j > 0)
          A_aux = sum(pmax(TPD_i[shared_aux], TPD_j[shared_aux])) - 
            O_aux
          only_in_i_aux = which(TPD_i > 0 & TPD_j == 
                                  0)
          B_aux = sum(TPD_i[only_in_i_aux])
          only_in_j_aux = which(TPD_i == 0 & TPD_j > 
                                  0)
          C_aux = sum(TPD_j[only_in_j_aux])
          results_samp$dissimilarity[i, j] = results_samp$dissimilarity[j, 
                                                                        i] = 1 - O_aux
          if (results_samp$dissimilarity[j, i] == 0) {
            results_samp$P_non_shared[i, j] = NA
            results_samp$P_non_shared[j, i] = NA
            results_samp$P_shared[i, j] = NA
            results_samp$P_shared[j, i] = NA
          }
          else {
            results_samp$P_non_shared[i, j] = results_samp$P_non_shared[j, 
                                                                        i] = (2 * min(B_aux, C_aux))/(A_aux + 
                                                                                                        2 * min(B_aux, C_aux))
            results_samp$P_shared[i, j] = results_samp$P_shared[j, 
                                                                i] = 1 - results_samp$P_non_shared[i, 
                                                                                                   j]
          }
        }
        if (i == j) {
          results_samp$dissimilarity[i, j] = 0
        }
      }
    }
    return(results_samp)
  }
  if (TPDType == "Communities") {
    results$communities = Calc_dissim(TPDc)
  }
  if (TPDType == "Populations") {
    results$populations = Calc_dissim(TPDs)
  }
  class(results) = "OverlapDiss"
  return(results)
}

# FD_index_clades: FD indices (FRic, FEve, FDiv) of clade (need to be specify,'order','family'),
# and null model
# Parallelized version of FD_index_clades
# Parallelized and memory-optimized version of FD_index_clades
FD_index_clades <- function(species, traits, taxo, TPDsAux, clade,type,
                                     index = c("FRic", "All"),dims=FALSE,
                                     boot = 1000, saveFile = NULL, n_cores = 4) {
  require(TPD)
  require(future.apply)
  require(progress)
  # Match species
  species <- species[species %in% TPDsAux$data$species]
  taxo <- taxo[taxo$genusspecies %in% species, ]
  
  # Create occurrence matrix by clade
  clade_col <- match.arg(clade, choices = c("order", "family", "genus"))
  biogeo_clade <- sort(table(taxo[[clade_col]]))
  occurences <- matrix(0, ncol = length(species), nrow = length(biogeo_clade),
                       dimnames = list(names(biogeo_clade), species))
  
  for (i in seq_along(biogeo_clade)) {
    clade_species <- taxo[which(taxo[[clade_col]] == names(biogeo_clade[i])), "genusspecies"]
    occurences[i, clade_species] <- 1
  }
  
  gc()
  n_communities <- nrow(occurences)
  names_comm <- rownames(occurences)
  metrics <- if (index == "FRic") "FRichness" else c("FRichness", "FEvenness", "FDivergence")
  
  functionTPDRichness <- if (dims) TPDRichness_large else TPDRichness
  functionTPDc <- if (dims) TPDc_large else TPDc
  
  cat(paste0("\n-------------------- Clade - Observed ", clade, " START --------------------\n"))
  
  # Observed values
  TPDc_occurences <- functionTPDc(TPDsAux, occurences)
  FD_obs <- functionTPDRichness(TPDc_occurences)
  
  if(!is.null(boot)){ 
    obs_values <- lapply(metrics, function(met) {
      out <- matrix(0, nrow = boot + 1, ncol = n_communities,
                    dimnames = list(c("occ", paste0("boot", 1:boot)), names_comm))
      out[1, names(FD_obs$communities[[met]])] <- FD_obs$communities[[met]]
      out
    })
    names(obs_values) <- gsub("FRichness", "FRic", gsub("FEvenness", "FEve", gsub("FDivergence", "FDiv", metrics)))
    
  }else{
    obs_values <- lapply(metrics, function(met) {
      out <- matrix(0, nrow = 1, ncol = n_communities,
                    dimnames = list("occ", names_comm))
      out[1, names(FD_obs$communities[[met]])] <- FD_obs$communities[[met]]
      out
    })
    names(obs_values) <- gsub("FRichness", "FRic", gsub("FEvenness", "FEve", gsub("FDivergence", "FDiv", metrics)))
    
  }

  rm(TPDc_occurences)
  gc()
  cat(paste0("\n-------------------- Clade - Observed ", clade, " DONE --------------------\n"))
  
if(!is.null(boot)){  # Set up parallel plan
  options(future.globals.maxSize = 12 * 1024^3)  # 4 GiB
  old_plan <- plan(multisession, workers = n_cores)
  on.exit(plan(old_plan), add = TRUE)
  
  # Parallel bootstraps
  boot_res <- local({
    TPDsAux_local <- TPDsAux  # Bind locally
    occurences_local <- occurences
    names_comm_local <- names_comm
    metrics_local <- metrics
    
    future_lapply(1:boot, function(b) {
      Rand <- occurences_local
      for (i in seq_len(nrow(Rand))) {
        sp_present <- which(Rand[i, ] == 1)
        Rand[i, ] <- 0
        Rand[i, sample(ncol(Rand), length(sp_present))] <- 1
      }
      
      TPDc_rand <- if (dims) TPDc_large(TPDsAux_local, Rand) else TPDc(TPDsAux_local, Rand)
      FD_rand <- if (dims) TPDRichness_large(TPDc_rand) else REND(TPDc_rand)
      
      met_values <- lapply(metrics_local, function(met) {
        out <- rep(NA, length(names_comm_local))
        names(out) <- names_comm_local
        if (!is.null(FD_rand$communities[[met]])) {
          out[names(FD_rand$communities[[met]])] <- FD_rand$communities[[met]]
        }
        out
      })
      names(met_values) <- metrics_local
      return(met_values)
    }, future.seed = TRUE)
  })
  
  
  # Fill in bootstrap values
  for (m in seq_along(metrics)) {
    met_name <- names(obs_values)[m]
    for (b in seq_len(boot)) {
      obs_values[[met_name]][b + 1, ] <- boot_res[[b]][[m]]
    }
  }
}  
  if (!is.null(saveFile)) {
    saveRDS(obs_values, saveFile)
  }
  
  return(obs_values)
}


Dissim_index_clades = function(species, traits, taxo, TPDsAux, clade,type, dim=F , saveFile) {
  # Filtrer les espèces présentes dans les TPDs
  species = species[species %in% TPDsAux$data$species]
  taxo = taxo[taxo$genusspecies %in% species, ]
  
  # Initialiser la matrice d'occurrence
  if (clade == "order") {
    clade_tab = sort(table(taxo$order))
    clade_names = names(clade_tab)
    occurences = matrix(0, ncol = length(species), nrow = length(clade_tab),
                        dimnames = list(clade_names, species))
    for (i in seq_along(clade_names)) {
      all_sp = taxo[which(taxo$order == clade_names[i]), ]$genusspecies
      occurences[i, all_sp] = 1
    }
  } else if (clade == "family") {
    clade_tab = sort(table(taxo$family))
    clade_names = names(clade_tab)
    occurences = matrix(0, ncol = length(species), nrow = length(clade_tab),
                        dimnames = list(clade_names, species))
    for (i in seq_along(clade_names)) {
      all_sp = taxo[which(taxo$family == clade_names[i]), ]$genusspecies
      occurences[i, all_sp] = 1
    }
  } else if (clade == "genus") {
    clade_tab = sort(table(taxo$genus))
    clade_names = names(clade_tab)
    occurences = matrix(0, ncol = length(species), nrow = length(clade_tab),
                        dimnames = list(clade_names, species))
    for (i in seq_along(clade_names)) {
      all_sp = taxo[which(taxo$genus == clade_names[i]), ]$genusspecies
      occurences[i, all_sp] = 1
    }
  } else {
    stop("Clade level not recognized. Use 'order', 'family', or 'genus'.")
  }
  
  gc()
  cat(paste0("\n -------------------- Clade - Dissimilarity ", clade, " START -------------------- \n"))

  functionTPDDissim <- if (dim) dissim_large else TPD::dissim
  functionTPDc <- if (dim) TPDc_large else TPD::TPDc
  
  # Calcul de la dissimilarité fonctionnelle
  TPDc_occurences = functionTPDc(TPDsAux, occurences)
  FDiss_occurences = functionTPDDissim(TPDc_occurences)
  
  if (!is.null(saveFile)) {
    saveRDS(FDiss_occurences, saveFile)
  }
  
  return(FDiss_occurences)
}


biPlot_Correlation_orig <- function(d1, d2, realm = NULL, title=NA, xlab=NA, ylab=NA,
                                    xlim=c(0,1), ylim=c(0,1), angle=NA, spear = TRUE, place="topleft",
                                    colpoints = NULL, alpha = 1, ratio = FALSE, smi = FALSE, size=1, vert = FALSE, abl = FALSE){
  
  if (!is.null(realm)) {
    realm <- factor(realm)
    palette_realm <- paletteer::paletteer_d('khroma::sunset',11)[c(1,3,4,5,7,9,11)]
    colpoints <- palette_realm[as.numeric(realm)]
    your_data <- data.frame(
      data1 = d1,
      data2 = d2,
      realm = realm,
      colpoints = colpoints,
      alpha = rep(alpha, length(d1)),
      size = rep(size, length(d1))
    )
  } else {
    if (is.null(colpoints)) colpoints <- "grey9"
    your_data <- data.frame(
      data1 = d1,
      data2 = d2,
      colpoints = colpoints,
      alpha = rep(alpha, length(d1)),
      size = rep(size, length(d1))
    )
  }
  
  # Définir noms corrects
  if (!is.null(realm)) {
    names(your_data) <- c("data1", "data2", "realm", "colpoints", "alpha", "size")
  } else {
    names(your_data) <- c("data1", "data2", "colpoints", "alpha", "size")
  }
  
  your_data <- na.omit(your_data)
  
  pearson_test <- cor.test(na.omit(your_data$data1),
                           na.omit(your_data$data2), method = "spearman")
  
  cor_coefficient <- round(pearson_test$estimate, 2)
  p_value <- pearson_test$p.value
  
  significance <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
  
  # Construction du plot
  if (!is.null(realm)) {
    plot_res <- ggplot(your_data, aes(x = data1, y = data2, color = realm)) +
      geom_point(alpha = your_data$alpha, size = your_data$size) +
      scale_color_manual(values = palette_realm, name = "")+
      theme_classic(base_size = 11) +
      xlim(xlim) + ylim(ylim) +
      ggtitle(label = title) +
      theme(
        plot.margin = unit(c(0,0.3,0.4,0.4), "cm"),
        plot.title = element_text(size = 13, face = "bold", hjust = -0.14, vjust = -4),
        text = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 8, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 8, r = 8, b = 0, l = 0)),
        axis.ticks.length = unit(-0.15, "cm"),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA),
        legend.title = element_blank()
      ) +
      xlab(xlab) + ylab(ylab)
  } else {
    plot_res <- ggplot(your_data, aes(x = data1, y = data2)) +
      geom_point(color = your_data$colpoints, alpha = your_data$alpha, size = your_data$size) +
      scale_color_identity()  +
      theme_classic(base_size = 11) +
      xlim(xlim) + ylim(ylim) +
      ggtitle(label = title) +
      theme(
        plot.margin = unit(c(0,0.3,0.4,0.4), "cm"),
        plot.title = element_text(size = 13, face = "bold", hjust = -0.14, vjust = -4),
        text = element_text(size = 12),
        axis.title.y = element_text(margin = margin(t = 8, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 8, r = 8, b = 0, l = 0)),
        axis.ticks.length = unit(-0.15, "cm")
      ) +
      xlab(xlab) + ylab(ylab)
  }
  
  if(ratio){
    plot_res <- plot_res + geom_segment(color="grey19",x = xlim[1], y = 1, xend = xlim[2], yend = 1,linetype = 3)
  }
  
  if(smi){
    plot_res <- plot_res + geom_segment(color="grey19",x = xlim[1], y = 0, xend = xlim[2], yend = 0,linetype = 3)
  }
  
  if(spear){
    title_with_spearman <- paste0("Rho = ", cor_coefficient, " ", significance)
    if(place == "topleft"){
      plot_res <- plot_res + annotate("text", x = xlim[1], y = ylim[2], label = title_with_spearman, hjust = 0.02, vjust = 0.3)
    }
    if(place == "topright"){
      plot_res <- plot_res + annotate("text", x = xlim[2], y = ylim[2], label = title_with_spearman, hjust = 1, vjust = 0.4)
    }
    if(place == "bottomright"){
      plot_res <- plot_res + annotate("text", x = xlim[2], y = ylim[1], label = title_with_spearman, hjust = 1, vjust = -0.3)
    }
  }
  
  if(vert){
    plot_res <- plot_res + geom_segment(color="grey19",x = 0, y = ylim[1], xend = 0, yend = ylim[2],linetype = 3)
  }
  
  if(abl){
    plot_res <- plot_res + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey38", linewidth = 0.6) + annotate("text", 
                                                                                                                                    x = min(xlim), 
                                                                                                                                    y = min(xlim), 
                                                                                                                                    label = "1:1", 
                                                                                                                                    hjust = -0.1, vjust = -1, 
                                                                                                                                    size = 4, 
                                                                                                                                    color = "black",
                                                                                                                                    angle = 45)
  }
  
  if(!is.na(angle)){
    plot_res <- plot_res +
      scale_y_continuous(expand = c(0.03, 0.03), limits=ylim, breaks = seq(0, 180, by = 30)) +
      scale_x_continuous(expand = c(0.03, 0.03), limits=xlim, breaks = seq(0, 180, by = 30)) +
      geom_segment(x = 90, y = 0, xend = 90, yend = 180, color = "grey19", linetype = 3) +
      geom_segment(x = 0, y = 90, xend = 180, yend = 90, color = "grey19", linetype = 3)
  }
  
  plot_res
}

### 4. Graphics ----
biPlot_Correlation <- function(d1, d2, xlab = "", ylab = "", title = "",
                               xlim = NULL, ylim = NULL,
                               abl = TRUE,
                               orders = NULL,
                               order_colors = NULL,
                               orders_to_label = NULL) {
  library(ggplot2)
  library(ggrepel)
  
  df <- data.frame(d1 = d1, d2 = d2)
  if (!is.null(orders)) df$order <- orders
  if (!is.null(orders_to_label)) {
    df$label <- ifelse(rownames(df) %in% orders_to_label, rownames(df), NA)
  }
  
  p <- ggplot(df, aes(x = d1, y = d2)) +
    geom_point(aes(color = order), size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey38", linewidth = 0.6) +
    annotate("text", 
             x = min(xlim), 
             y = min(xlim), 
             label = "1:1", 
             hjust = -0.1, vjust = -1, 
             size = 4, 
             color = "black",
             angle = 45)+
    coord_fixed() +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 1) # réduit marge gauche
    )
  
  # Calcul de la régression
  reg <- lm(d2 ~ d1, data = df)
  r2 <- summary(reg)$r.squared
  
  if(summary(reg)$coef[2,4]<0.05){
   p <- p + geom_smooth(method = "lm", se = F, color = "black", linewidth = 1.2) 
  }else{
    p <- p + geom_smooth(method = "lm", se = F, color = "grey56", linewidth = 1.2, linetype = "dashed")
  }
  
  # Ajouter texte R²
  annot_label <- bquote(italic(R)^2 == .(format(r2, digits = 2)))
  p <- p + annotate("text", 
                    x = min(xlim), 
                    y = max(xlim), 
                    label = as.expression(annot_label), 
                    hjust = 0, vjust = 0, 
                    size = 4.5, fontface = "italic")
  
  if (!is.null(order_colors)) {
    p <- p + scale_color_manual(values = order_colors)
  }
  
  # Ajout des lignes de liaison vers les labels
  if (!is.null(orders_to_label)) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = label),
        box.padding = 0.4,
        point.padding = 0.2,
        segment.color = "gray30",
        segment.size = 0.4,
        size = 3.5,
        na.rm = TRUE,
        max.overlaps = Inf
      )
  }
  
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(ylim)) p <- p + ylim(ylim)
  
  return(p)
}





# 
boxPlotCdi_all=function(cdi_all,title="NA",values){
  
  your_data = data.frame(
    Group = rep(colnames(cdi_all),each=nrow(cdi_all)),
    Value =unlist(cdi_all)  # Replace with your actual values
  )
  
  your_data$Group = factor(your_data$Group,
                           levels = unique(your_data$Group))
  
  
  # Create the boxplot with jitter using ggplot2
  boxplot_res = ggplot(your_data, aes(x = Group, y = Value,fill=Group)) +
    coord_flip()+
    geom_boxplot(position = position_dodge(width = 0.15), width = 0.7,na.rm = T,outlier.color="grey69",outlier.size = 0.5,
                 linewidth=0.1,color="grey7") +
    stat_summary(fun = mean, geom="point", colour="black", size=2,
                 position = position_dodge2(width = 0.75))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "#D8511DFF", size = 0.8)+
    scale_fill_manual(values = values)+
    theme_classic(  base_size = 11,
                    base_family = "",
                    base_line_size = 0.1,
                    base_rect_size = 0.1)+
    ggtitle(label=title)+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),plot.title = element_text(size = 13, face = "bold"),legend.position="none")+
    xlab("")+ylab("")
  boxplot_res
}

boxplot_g = function(group,values,levels=c("NORMAL","OBSERVED"),xlab="",ylab="",
                     title = ""){
  your_data <- data.frame(
    Group = group,
    Value = values
    
  )
  
  your_data$Group = factor(your_data$Group,
                           levels = levels)
  # Create the boxplot with jitter using ggplot2
  boxplot_res = ggplot(your_data, aes(x = Group, y = Value,fill=Group)) +
    geom_boxplot(position = position_dodge(width = 0.15), width = 0.7,na.rm = T,outlier.color="grey69",outlier.size = 0.5,
                 linewidth=0.1,color="grey7")+
    scale_fill_paletteer_d("MetBrewer::Austria")+
    theme_classic(  base_size = 11,
                    base_family = "",
                    base_line_size = 0.1,
                    base_rect_size = 0.1)+
    ggtitle(label=title)+
    theme(plot.margin = unit(c(0,0.3,0.4,0.4), "cm"),
          plot.title = element_text(size = 13, face = "bold",hjust = -0.1,vjust = 0),
          text=element_text(size=12),
          axis.title.y = element_text(margin = margin(t = 8, r = 7, b = 10, l = 0),vjust = .5),
          axis.title.x = element_text(margin = margin(t = 8, r = 8, b = 0, l = 0)),
          axis.ticks.length  = unit(-0.15, "cm"),legend.position="none")+
    xlab(xlab)+ylab(ylab)
  boxplot_res
}

#
Pairs_Correlation=function(data_pairs,type,title=NA,l){
  col=rev(paletteer_c("ggthemes::Red-Blue Diverging", 100))
  maxXY = length(type)
  plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(0,length(type)),
       ylim=c(0,length(type)))
  # par(mai=c(0,1,1,0),mar=c(0,1,1,1))
  # mtext(2,at=c(0.5:6.5),text=rev(type), line = -1,cex=1,las=1)
  # mtext(3,at=c(0.5:6.5),text=type, line =-1,cex=1,las=1)
  mtext(3,adj=0,text=paste0("(",l,") ",title), line = -0.5,cex=1.25,las=1)
  
  for(j in 1:(length(type))){
    x=data_pairs[[j]]
    for(i in 1:(length(type))){
      y=data_pairs[[i]]
      if(i>j){
        
        your_data = data.frame(
          data1 = x,
          data2 = y
        )
        names(your_data) = c("data1","data2")
        pearson_test = cor.test(your_data$data1, your_data$data2, method = "spearman")
        
        # Extract correlation coefficient and p-value
        cor_coefficient = round(pearson_test$estimate, 2)
        p_value = pearson_test$p.value
        
        # Determine significance level
        significance = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
        
        # Create title string
        title_with_spearman = paste0(cor_coefficient, " ", significance)
        
        rect(j-1,((maxXY+1)-i)-1,j,((maxXY+1) - i),col=col[round(cor_coefficient*100)])
        text(j-1,((maxXY+0.5)-i),title_with_spearman,pos=4,cex=1.1)
      }
      if(i==j){
        rect(j-1,((maxXY+1)-i)-1,j,((maxXY+1)-i),col="white")
        text(j-0.8,((maxXY+0.5)-i),type[i],pos=4,cex=1.1)
      }
    }
  } 
}

#
Pairs_Correlation_Diff <- function(data1, data2, trait_names, 
                                   title1, title2, title = NA, label = NA,
                                   col = rev(paletteer_c("ggthemes::Red-Blue Diverging", 100)),
                                   show_table = FALSE, diag_col = "white") {
  
  n_traits <- length(trait_names)
  
  # Set up empty plot
  plot(0, type = "n", axes = FALSE, ann = FALSE,
       xlim = c(0, n_traits), ylim = c(0, n_traits))
  
  # Titles
  if (!is.na(title)) {
    mtext(3, adj = 0, text = paste0("(", label, ") ", title), line = 1, cex = 1.25)
  }
  mtext(2, adj = 0.5, text = title1, line = -0.5, cex = 1.25)
  mtext(3, adj = 0.5, text = title2, line = -0.5, cex = 1.25)
  
  # Table to store results (optional)
  cor_table <- matrix(NA, ncol = 2,
                      nrow = choose(n_traits, 2),
                      dimnames = list(
                        apply(combn(trait_names, 2), 2, paste, collapse = " - "),
                        c(title1, title2)
                      ))
  
  k <- 1
  for (j in seq_len(n_traits)) {
    for (i in seq_len(n_traits)) {
      if (i > j) {
        # Correlation for dataset 1
        test1 <- cor.test(data1[[j]], data1[[i]], method = "spearman")
        cor1 <- round(test1$estimate, 2)
        sig1 <- ifelse(test1$p.value < 0.001, "***",
                       ifelse(test1$p.value < 0.01, "**",
                              ifelse(test1$p.value < 0.05, "*", "ns")))
        label1 <- paste0(cor1, " ", sig1)
        col_index1 <- round((cor1 + 1) * 50)
        rect(j - 1, n_traits - i, j, n_traits - i + 1, col = col[col_index1])
        text(j - 0.9, n_traits - i + 0.5, label1, pos = 4, cex = 1.1)
        cor_table[paste(trait_names[j], trait_names[i], sep = " - "), 1] <- label1
        
        # Correlation for dataset 2
        test2 <- cor.test(data2[[j]], data2[[i]], method = "spearman")
        cor2 <- round(test2$estimate, 2)
        sig2 <- ifelse(test2$p.value < 0.001, "***",
                       ifelse(test2$p.value < 0.01, "**",
                              ifelse(test2$p.value < 0.05, "*", "ns")))
        label2 <- paste0(cor2, " ", sig2)
        col_index2 <- round((cor2 + 1) * 50)
        rect(i - 1, n_traits - j, i, n_traits - j + 1, col = col[col_index2])
        text(i - 0.9, n_traits - j + 0.5, label2, pos = 4, cex = 1.1)
        cor_table[paste(trait_names[j], trait_names[i], sep = " - "), 2] <- label2
      }
      
      if (i == j) {
        rect(j - 1, n_traits - i, j, n_traits - i + 1, col = diag_col)
        text(j - 0.8, n_traits - i + 0.5, trait_names[i], pos = 4, cex = 1.1)
      }
    }
  }
  
  if (show_table) return(cor_table)
}


#
Pairs_Correlation_Lat=function(data_pairs1,data_pairs2,Lat,Lat1,type,title1,title2,table = F){
  
  cor_coefficient_1 = matrix(NA, ncol = 2,
                             nrow = (length(data_pairs1)*(length(data_pairs1)-1))/2,
                             dimnames = list(apply(combn(type,2),2,function(x){paste0(x[1]," - ",x[2])}),
                                             c(title1,title2)))
  col=rev(paletteer_c("ggthemes::Red-Blue Diverging", 100))
  maxXY = length(type)
  plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(0,length(type)),
       ylim=c(0,length(type)))
  # par(mai=c(0,1,1,0),mar=c(0,1,1,1))
  # mtext(2,at=c(0.5:6.5),text=rev(type), line = -1,cex=1,las=1)
  # mtext(3,at=c(0.5:6.5),text=type, line =-1,cex=1,las=1)
  # mtext(3,adj=0,text=paste0("(",l,") ",title), line = -0.5,cex=1.25,las=1)
  mtext(2,adj=0.5,text=title1, line =-0.5,cex=1.25,las=0)
  mtext(3,adj=0.5,text=title2, line =-0.5,cex=1.25,las=0)
  
  for(j in 1:(length(type))){
    x=data_pairs1[[j]]
    x1=data_pairs2[[j]]
    for(i in 1:(length(type))){
      xy=data_pairs1[[i]]
      xy1=data_pairs2[[i]]
      y=Lat
      y1=Lat1
      if(i>j){
        
        your_data = data.frame(
          data1 = x-xy,
          data2 = y
        )
        names(your_data) = c("data1","data2")
        pearson_test = cor.test(your_data$data1, your_data$data2, method = "spearman")
        
        # Extract correlation coefficient and p-value
        cor_coefficient = round(pearson_test$estimate, 2)
        p_value = pearson_test$p.value
        
        # Determine significance level
        significance = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
        
        # Create title string
        if( significance == "ns"){
          title_with_spearman = paste0(significance)
        }else{
          title_with_spearman = paste0(cor_coefficient, " ", significance)
        }
        
        rect(j-1,((maxXY+1)-i)-1,j,((maxXY+1) - i),col=col[round(cor_coefficient*100)])
        text(j-1,((maxXY+0.5)-i),title_with_spearman,pos=4,cex=1.1)
        
        cor_coefficient_1[paste0(type[j]," - ",type[i]),1] = title_with_spearman
        
        ##
        your_data = data.frame(
          data1 = x1-xy1,
          data2 = y1
        )
        names(your_data) = c("data1","data2")
        pearson_test = cor.test(your_data$data1, your_data$data2, method = "spearman")
        
        # Extract correlation coefficient and p-value
        cor_coefficient = round(pearson_test$estimate, 2)
        p_value = pearson_test$p.value
        
        # Determine significance level
        significance = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
        
        # Create title string
        if( significance == "ns"){
          title_with_spearman = paste0(significance)
        }else{
          title_with_spearman = paste0(cor_coefficient, " ", significance)
        }
        
        
        rect(i-1,((maxXY+1) - j)-1,i,((maxXY+1) - j),
             col=col[round(cor_coefficient*100)])
        text(i-1,((maxXY+0.5) - j),title_with_spearman,pos=4,cex=1.1)
        cor_coefficient_1[paste0(type[j]," - ",type[i]),2] = title_with_spearman
      }
      if(i==j){
        rect(j-1,((maxXY+1)-i)-1,j,((maxXY+1)-i),col="white")
        text(j-0.8,((maxXY+0.5)-i),type[i],pos=4,cex=1.1)
      }
      
    }
  } 
  if(table){
    return(cor_coefficient_1)
  }
}


## 5. Biogeo ----
# dggs7 <- dgconstruct(res=7, metric=FALSE, resround='down',topology="HEXAGON")
# geog7 = readRDS("data/geographic borders/geogInfo_dggs7.RDS")
# R functions to extract and estimate coordinates of each species
# The functions below are used to extract the raster of each species and estimate the coordinates of each occurrence for each species

### 1. Create and save raster of each species
extract_shp <- function(fileShp,pathSave){
  rasterSp::rasterizeRange(dsn=fileShp,id ="new",
                           resolution=0.5, origin=1, presence=c(1, 2, 3),
                           save=TRUE, path=pathSave)
}

#### 2. Estimate coordinates of each occurrence for each species
est_coord_sp = function(filedir,fileList,pathSave){
  coordsSpecies = list()
  for(i in 1:length(fileList)){
    cat(paste("\rSpecies ", i, "of ", length(fileList), "        \r"))
    sfilename = paste0(filedir, "/SpeciesData/", fileList[i])
    imported_raster = raster(sfilename)
    # Get rows and columns of cells with presences:
    coordsAux = which(!is.na(as.matrix(imported_raster)), arr.ind=T)
    latAux = yFromRow(imported_raster, row = coordsAux[,1])
    longAux = xFromCol(imported_raster, col = coordsAux[,2])
    latLongMat = cbind(longAux, latAux)
    colnames(latLongMat) = c("long", "lat")
    coordsSpecies[[i]] = latLongMat
    spNameAux = unlist(strsplit(fileList[i], split="_"))
    spNameAux = paste(spNameAux[1], spNameAux[2], sep="_")
    names(coordsSpecies)[i] = spNameAux
  }
  saveRDS(coordsSpecies, file = pathSave)
  return(coordsSpecies)
}
#### 3.  ASSIGN HEXAGONS FROM DGGS to coordinates:
assign_hexagons = function(coordsSpecies, dggs7, sitesInRealms7,savePath){
  speciesRemove = numeric()
  for(i in 1:length(coordsSpecies)){
    cat(paste("\rSpecies ", i, "of ", length(coordsSpecies), "        \r"))
    if(nrow(coordsSpecies[[i]]) > 0){
      spdggs7 = dgGEO_to_SEQNUM(dggs7, coordsSpecies[[i]][, "long"], coordsSpecies[[i]][, "lat"])$seqnum
      keepPoints = which(spdggs7 %in% sitesInRealms7)
      spdggs7 = spdggs7[keepPoints]
      if(length(spdggs7) > 0){
        coordsSpecies[[i]] = coordsSpecies[[i]][keepPoints, , drop=FALSE]
        coordsSpecies[[i]] = cbind(coordsSpecies[[i]], spdggs7)
        colnames(coordsSpecies[[i]]) = c("long", "lat", "dggs7")
      }else{
        speciesRemove = c(speciesRemove, names(coordsSpecies)[i])
      }
    }
  }
  for(sp in 1:length(speciesRemove)){
    coordsSpecies[[speciesRemove[sp]]] = NULL 
  }
  saveRDS(coordsSpecies, file = savePath)
  return(coordsSpecies)
}
#### 4. Lets reorder the data. We want a list for each resolution, with an element for each 
### hexagon, listing the species present there.
reorder_hexagon = function(coordsSpecies,savePath){
  uniquedggs7 = numeric()
  for(i in 1:length(coordsSpecies)){
    cat(paste("\rSpecies ", i, "of ", length(coordsSpecies), "        \r"))
    if(nrow(coordsSpecies[[i]]) > 0){
      uniquedggs7 = c(uniquedggs7, coordsSpecies[[i]][, "dggs7"])
    }
  }
  uniquedggs7 = unique(uniquedggs7)
  saveRDS(uniquedggs7, file = savePath)
  return(uniquedggs7)
}
### 5. Create a matrix with the number of species in each hexagon
siteswithspecies = function(coordsSpecies, uniquedggs7,savePath){
  sitesdggs7 = list()
  for(i in 1:length(uniquedggs7)){
    sitesdggs7[[i]] = character()
  }
  names(sitesdggs7) = uniquedggs7
  for(i in 1:length(coordsSpecies)){# Fill species in list of hexagons
    cat(paste("\rSpecies ", i, "of ", length(coordsSpecies), "        \r"))
    if(nrow(coordsSpecies[[i]]) > 0){
      spAux = names(coordsSpecies)[i]
      dggs7Aux = as.character(unique(coordsSpecies[[i]][, "dggs7"]))
      for(site in 1:length(dggs7Aux)){
        sitesdggs7[[dggs7Aux[site]]] = c(sitesdggs7[[dggs7Aux[site]]], spAux)
      }
    }
  }
  saveRDS(sitesdggs7, file = savePath)
  return(sitesdggs7)
}
### 6. Compute the TPDs
TPDs_compute = function(TraitsPCA, sitesdggs7,savePath,sampleComms = 800,  alphaUse = 0.95,
                        gridSize = 100){
  require(TPD)
  sdMatAux = numeric() 
  commsSampled = sample(1:length(sitesdggs7), sampleComms) 
  for(i in 1:sampleComms){
    cat(paste0("\r Community ", i," out of ", sampleComms, "\r"))
    sitesdggs7Aux = sitesdggs7[[commsSampled[i]]]
    if(length(sitesdggs7Aux) > 2){
      traitsAux = na.omit(TraitsPCA[which(rownames(TraitsPCA)%in%sitesdggs7Aux),,drop = F])
      if(nrow(traitsAux) > 2){
        sdMatAux = rbind(sdMatAux, sqrt(diag(Hpi.diag(traitsAux))))
      }
    }
  }
  
  sdMatUse = matrix(rep(colMeans(sdMatAux, na.rm=T), nrow(TraitsPCA)),
                     nrow=nrow(TraitsPCA), ncol=ncol(TraitsPCA), byrow=T)
  TPDssdggs7 = TPDsMean(species = rownames(TraitsPCA),
                           means = TraitsPCA,
                           sds = sdMatUse,
                           alpha = alphaUse,
                           n_divisions = gridSize)
  saveRDS(TPDssdggs7, savePath)
}

TPDs_compute_large = function(TraitsPCA, sitesdggs7,savePath,sampleComms = 800,  alphaUse = 0.95,
                        gridSize = 100){
  require(TPD)
  sdMatAux = numeric() 
  commsSampled = sample(1:length(sitesdggs7), sampleComms) 
  for(i in 1:sampleComms){
    cat(paste0("\r Community ", i," out of ", sampleComms, "\r"))
    sitesdggs7Aux = sitesdggs7[[commsSampled[i]]]
    if(length(sitesdggs7Aux) > 10){
      traitsAux = na.omit(TraitsPCA[which(rownames(TraitsPCA)%in%sitesdggs7Aux),,drop = F])
      if(nrow(traitsAux) > 10){
        sdMatAux = rbind(sdMatAux, sqrt(diag(Hpi.diag(traitsAux))))
      }
    }
  }
  
  sdMatUse = matrix(rep(colMeans(sdMatAux, na.rm=T), nrow(TraitsPCA)),
                    nrow=nrow(TraitsPCA), ncol=ncol(TraitsPCA), byrow=T)
  TPDssdggs7 = TPDsMean_large(species = rownames(TraitsPCA),
                        means = TraitsPCA,
                        sds = sdMatUse,
                        alpha = alphaUse,
                        n_divisions = gridSize)
  saveRDS(TPDssdggs7, savePath)
}

compute_cell_indices <- function(
    TraitsPCA, IUCN, threat, sitesdggs7, TPDs_sdggs7,
    functionREND, functionTPDc,
    indx = c("TD_bef","TD_aft","FRic_bef","FRic_aft",
             "FEve_bef","FEve_aft","FDiv_bef","FDiv_aft",
             "FRed_bef","FRed_aft","FUniq_bef","FUniq_aft",
             "Dissim","P_Shared"),
    progress = TRUE,
    chunk_size = 300,   # <= tune this; smaller = lower memory spike
    ...
) {
  # ---- helpers ----
  std_names <- function(x) gsub(" ", "_", as.character(x), fixed = TRUE)
  tick <- function(msg) if (progress) cat(sprintf("%s\n", msg))
  
  # 1) Threatened set (bugfix)
  IUCN$scientificNameStd <- std_names(IUCN$scientificNameStd)
  IUCNSp <- IUCN %>%
    dplyr::filter(.data$redlistCategory %in% threat) %>%
    dplyr::pull(scientificNameStd) %>%
    unique()
  
  site_names <- names(sitesdggs7)
  sp_names   <- std_names(rownames(TraitsPCA))
  
  # Species universe from TPDs
  all_species <- if (!is.null(TPDs_sdggs7$data$species)) {
    std_names(TPDs_sdggs7$data$species)
  } else if (!is.null(names(TPDs_sdggs7$TPDs))) {
    std_names(names(TPDs_sdggs7$TPDs))
  } else stop("Could not find species names in TPDs_sdggs7.")
  
  # 2) Build communities (cells × species in TraitsPCA)
  tick("Step 1/6: building community matrix …")
  build_row <- function(i) {
    spIn <- std_names(sitesdggs7[[i]])
    spInTraits <- spIn[spIn %in% sp_names]
    row <- integer(length(sp_names))
    if (length(spInTraits)) row[match(spInTraits, sp_names)] <- 1L
    row
  }
  idx <- seq_along(sitesdggs7)
  if (progress && requireNamespace("pbapply", quietly = TRUE)) {
    Comm_list <- pbapply::pblapply(idx, build_row)
  } else {
    Comm_list <- lapply(idx, build_row)
  }
  Comm <- do.call(rbind, Comm_list)
  rownames(Comm) <- site_names
  colnames(Comm) <- sp_names
  
  # 3) Reduce to TPD species + drop never-present species and empty cells
  tick("Step 2/6: reducing matrix (species/cells) …")
  keep_cols <- colnames(Comm) %in% all_species
  CommsTPD  <- Comm[, keep_cols, drop = FALSE]
  
  # Drop species never present in any cell (colSums == 0)
  cs <- colSums(CommsTPD)
  CommsTPD <- CommsTPD[, cs > 0, drop = FALSE]
  
  # Create "after" scenario and drop threatened spp
  CommsTPDExtinctions <- CommsTPD
  thr_mask <- colnames(CommsTPDExtinctions) %in% IUCNSp
  if (any(thr_mask)) CommsTPDExtinctions[, thr_mask] <- 0L
  
  # Identify non-empty cells in each scenario
  rows_bef <- which(rowSums(CommsTPD) > 0)
  rows_aft <- which(rowSums(CommsTPDExtinctions) > 0)
  
  # For speed, also drop totally empty rows from both matrices
  any_rows <- sort(unique(c(rows_bef, rows_aft)))
  CommsTPD_sub  <- CommsTPD[any_rows, , drop = FALSE]
  CommsTPDE_sub <- CommsTPDExtinctions[any_rows, , drop = FALSE]
  
  tick(sprintf("   kept %d/%d cells; %d/%d species",
               nrow(CommsTPD_sub), nrow(Comm),
               ncol(CommsTPD_sub), ncol(Comm)))
  
  # 4) Chunked TPDc + REND
  run_chunked <- function(comm_mat, label) {
    n <- nrow(comm_mat)
    out_FRic <- out_FEve <- out_FDiv <- out_FRed <- out_FUni <- out_Diss <- out_Psha <- rep(NA_real_, n)
    
    if (n == 0) return(list(FRic=out_FRic))
    
    chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    if (progress) pb <- utils::txtProgressBar(min = 0, max = length(chunks), style = 3)
    for (k in seq_along(chunks)) {
      rows <- chunks[[k]]
      if (progress) utils::setTxtProgressBar(pb, k)
      # guard: if whole chunk is empty (shouldn't happen after filtering), skip
      if (all(rowSums(comm_mat[rows, , drop = FALSE]) == 0)) next
      
      TPDc_k <- functionTPDc(TPDs = TPDs_sdggs7, sampUnit = comm_mat[rows, , drop = FALSE])
      REND_k <- functionREND(TPDc = TPDc_k)
      
      out_FRic[rows] <- REND_k$communities$FRichness
    }
    if (progress) close(pb)
    list(FRic=out_FRic)
  }
  
  tick("Step 3/6: TPDc/REND for BEFORE (chunked) …")
  bef <- run_chunked(CommsTPD_sub, "before")
  
  tick("Step 4/6: TPDc/REND for AFTER (chunked) …")
  aft <- run_chunked(CommsTPDE_sub, "after")
  
  # 5) Assemble back to full list of cells (fill NA for cells dropped as empty)
  tick("Step 5/6: assembling outputs …")
  TD_bef_full <- rowSums(CommsTPD)
  TD_aft_full <- rowSums(CommsTPDExtinctions)
  
  # allocate full-length vectors
  FR_bef_full <- FR_aft_full <- FE_bef_full <- FE_aft_full <- FDiv_bef_full <- FDiv_aft_full <-
    FRed_bef_full <- FRed_aft_full <- FUni_bef_full <- FUni_aft_full <- Dissim_full <- Pshare_full <-
    rep(NA_real_, nrow(CommsTPD))
  
  # insert chunked results back into the positions used
  FR_bef_full[any_rows]   <- bef$FRic
  FE_bef_full[any_rows]   <- bef$FEve
  FDiv_bef_full[any_rows] <- bef$FDiv
  FRed_bef_full[any_rows] <- bef$FRed
  FUni_bef_full[any_rows] <- bef$FUni
  Dissim_full[any_rows]   <- bef$Dissim
  Pshare_full[any_rows]   <- bef$P_Shared
  
  FR_aft_full[any_rows]   <- aft$FRic
  FE_aft_full[any_rows]   <- aft$FEve
  FDiv_aft_full[any_rows] <- aft$FDiv
  FRed_aft_full[any_rows] <- aft$FRed
  FUni_aft_full[any_rows] <- aft$FUni
  # Dissim / Pshare typically “before”; keep as NA or compute separately for "after" if needed
  
  out <- data.frame(
    TD_bef   = TD_bef_full,
    TD_aft   = TD_aft_full,
    FRic_bef = FR_bef_full,
    FRic_aft = FR_aft_full,
    FEve_bef = FE_bef_full,
    FEve_aft = FE_aft_full,
    FDiv_bef = FDiv_bef_full,
    FDiv_aft = FDiv_aft_full,
    FRed_bef = FRed_bef_full,
    FRed_aft = FRed_aft_full,
    FUniq_bef = FUni_bef_full,
    FUniq_aft = FUni_aft_full,
    Dissim   = Dissim_full,
    P_Shared = Pshare_full,
    check.names = FALSE
  )
  rownames(out) <- rownames(CommsTPD)
  
  # keep requested columns in required order
  missing_cols <- setdiff(indx, colnames(out))
  if (length(missing_cols)) out[missing_cols] <- NA_real_
  out <- out[, indx, drop = FALSE]
  
  tick("Step 6/6: done.")
  out
}






compute_cell_FRic_Null <- function(TraitsPCA, IUCN, threat, sitesdggs7, TPDs_sdggs7,
                                         functionREND = REND, functionTPDc = TPDc,
                                         boot = 100, 
                                         indx = c("TD_bef", "TD_aft", "FRic_bef", "FRic_aft"),
                                         n_cores = parallel::detectCores() - 1,
                                         save_prefix = "res_cell_FRic_null_",
                                         save_path = ".",  
                                         save_every = 25,
                                         batch_size = 500) {
  
  require(future.apply)
  require(progressr)
  future::plan(multicore, workers = n_cores)
  handlers("cli")
  
  IUCNSp <- names(IUCN[as.character(IUCN) %in% threat])
  sp_names <- rownames(TraitsPCA)
  site_names <- names(sitesdggs7)
  
  # === 1. Create full community matrix
  Comm <- matrix(0, nrow = length(site_names), ncol = length(sp_names),
                 dimnames = list(site_names, sp_names))
  for (i in seq_along(sitesdggs7)) {
    sp_present <- intersect(sitesdggs7[[i]], sp_names)
    Comm[i, sp_present] <- 1
  }
  
  # === 2. Filter valid communities
  valid_comm <- rowSums(Comm[, colnames(Comm) %in% TPDs_sdggs7$data$species]) > 2
  CommsTPD <- Comm[valid_comm, , drop = FALSE]
  CommsTPD <- CommsTPD[, colnames(CommsTPD) %in% TPDs_sdggs7$data$species, drop = FALSE]
  
  n_com <- nrow(CommsTPD)
  site_ids <- rownames(CommsTPD)
  
  # === 3. Preallocation
  resultIndex <- matrix(NA, ncol = length(indx), nrow = n_com,
                        dimnames = list(site_ids, indx))
  resultIndexRand <- matrix(NA, ncol = boot, nrow = n_com,
                            dimnames = list(site_ids, paste0("Boot", 1:boot)))
  
  # === 4. Batch indices
  batch_indices <- split(seq_len(n_com), ceiling(seq_len(n_com) / batch_size))
  
  # === 6. Bootstrap nulls par batch avec barre de progression
  cat("\nLancement du bootstrap (null model) par lots...\n")
  with_progress({
    p <- progressor(steps = boot * length(batch_indices))
    for (b in seq_len(boot)) {
      for (batch_id in seq_along(batch_indices)) {
        batch_rows <- batch_indices[[batch_id]]
        Rand_batch <- CommsTPD[batch_rows, , drop = FALSE]
        
        for (i in seq_len(nrow(Rand_batch))) {
          spp <- which(Rand_batch[i, ] == 1)
          Rand_batch[i, ] <- 0
          Rand_batch[i, sample(colnames(Rand_batch), length(spp))] <- 1
        }
        
        TPDc_rand <- functionTPDc(TPDs = TPDs_sdggs7, sampUnit = Rand_batch)
        FRic_rand <- functionREND(TPDc_rand)
        resultIndexRand[rownames(Rand_batch), b] <- FRic_rand$communities$FRichness
        
        p(message = sprintf("Bootstrap %d, batch %d", b, batch_id))
      }
      
      # === Sauvegarde temporaire
      if (b %% save_every == 0 || b == boot) {
        save_file <- file.path(save_path, paste0(save_prefix, "partial_boot_", b, ".RDS"))
        saveRDS(list(rand = resultIndexRand[, 1:b, drop = FALSE],
                     current_boot = b),
                save_file)
        cat(sprintf("\nSauvegarde à l’itération %d : %s\n", b, save_file))
      }
    }
  })
  
  return(rand = resultIndexRand)
}


compute_cell_FRic_Obs <- function(TraitsPCA, IUCN, threat, sitesdggs7, TPDs_sdggs7,
                                   functionREND = REND, functionTPDc = TPDc,
                                   boot = 100, 
                                   indx = c("TD_bef", "TD_aft", "FRic_bef", "FRic_aft"),
                                   n_cores = parallel::detectCores() - 1,
                                   save_prefix = "res_cell_FRic_null_",
                                  save_sufix = "_observed.RDS",
                                   save_path = ".",  
                                   save_every = 25,
                                   batch_size = 500) {
  
  require(future.apply)
  require(progressr)
  future::plan(multicore, workers = n_cores)
  handlers("cli")
  
  IUCNSp <- names(IUCN[as.character(IUCN) %in% threat])
  sp_names <- rownames(TraitsPCA)
  site_names <- names(sitesdggs7)
  
  # === 1. Create full community matrix
  Comm <- matrix(0, nrow = length(site_names), ncol = length(sp_names),
                 dimnames = list(site_names, sp_names))
  for (i in seq_along(sitesdggs7)) {
    sp_present <- intersect(sitesdggs7[[i]], sp_names)
    Comm[i, sp_present] <- 1
  }
  
  # === 2. Filter valid communities
  valid_comm <- rowSums(Comm[, colnames(Comm) %in% TPDs_sdggs7$data$species]) > 2
  CommsTPD <- Comm[valid_comm, , drop = FALSE]
  CommsTPD <- CommsTPD[, colnames(CommsTPD) %in% TPDs_sdggs7$data$species, drop = FALSE]
  
  n_com <- nrow(CommsTPD)
  site_ids <- rownames(CommsTPD)
  
  # === 3. Preallocation
  resultIndex <- matrix(NA, ncol = length(indx), nrow = n_com,
                        dimnames = list(site_ids, indx))
  resultIndexRand <- matrix(NA, ncol = boot, nrow = n_com,
                            dimnames = list(site_ids, paste0("Boot", 1:boot)))
  
  # === 4. Batch indices
  batch_indices <- split(seq_len(n_com), ceiling(seq_len(n_com) / batch_size))
  
  # === 5. Observed values par batch avec barre de progression
  cat("\nCalcul des indices fonctionnels observés par lots...\n")
  with_progress({
    p <- progressor(steps = length(batch_indices))
    for (batch_id in seq_along(batch_indices)) {
      batch_rows <- batch_indices[[batch_id]]
      Com_batch <- CommsTPD[batch_rows, , drop = FALSE]
      
      TPDc_obs_batch <- functionTPDc(TPDs = TPDs_sdggs7, sampUnit = Com_batch)
      FRic_obs_batch <- functionREND(TPDc_obs_batch)
      
      resultIndex[rownames(Com_batch), ] <- cbind(rowSums(Com_batch), 
                                                  FRic_obs_batch$communities$FRichness)
      
      p(message = sprintf("Batch observé %d/%d", batch_id, length(batch_indices)))
    }
  })
  
  save_file <- file.path(save_path, paste0(save_prefix, save_sufix))
  saveRDS(resultIndex,
          save_file)
  
  return(obs = resultIndex)
}






#### 8. Mapping
cellToMap = function(res_cell7,dggs7,sitesdggs7){
  cellcenters   = dgSEQNUM_to_GEO(dggs7, as.numeric(names(sitesdggs7)))
  
  resultsdggs7 = cbind.data.frame(cell=as.numeric(names(sitesdggs7)),
                                  Long=cellcenters$lon_deg,
                                  Lat=cellcenters$lat_deg,
                                  res_cell7[names(sitesdggs7),])
  
  grid  = dgcellstogrid(dggs7, resultsdggs7[, "cell"])
  grid = merge(grid, resultsdggs7, by.x="seqnum",by.y="cell")
    return(grid)
}

# 
# saveRDS(grid, paste0(pathGr,"Clade_",targetGroupName[gr],"_",clade,"_Results_dggs7.rds"))
# 
# coT = cor.test(grid$Density[which(grid$propExt!=0)],grid$propExt[which(grid$propExt!=0)])
# corTh1 = round(coT$estimate,3)
# pCor1 = ifelse(coT$p.value>0.05,"ns",ifelse(coT$p.value>0.01,"*",ifelse(coT$p.value>0.001,"**","***")))
# 
# coT = cor.test(grid$ChgeFRic[which(grid$propExt!=0)],grid$propExt[which(grid$propExt!=0)])
# corTh2 = round(coT$estimate,3)
# pCor2 = ifelse(coT$p.value>0.05,"ns",ifelse(coT$p.value>0.01,"*",ifelse(coT$p.value>0.001,"**","***")))
# 
# grid[which(grid$propExt==0),"Density"] = NA
# grid[which(grid$propExt==0),"ChgeFRic"] = NA
# grid[which(grid$propExt==0),"ChgeFEve"] = NA
# grid[which(grid$propExt==0),"ChgeFDiv"] = NA
# grid[which(grid$propExt==0),"propExt"] = NA
# 
# 
# png(paste0("Figures/",targetGroupName[gr],"/",targetGroupName[gr],"_",clade,"_BiPloT.png"),width = 9000, height = 7000, units = "px", pointsize = 6,
#     bg = "white", res = 600, family = "")
# plot_erosion(grid,titre=clade,titre2= paste0("Rho : ",corTh1," ",pCor1),titre3 = paste0("Rho : ",corTh2," ",pCor2))
# dev.off()   

#plot_biclass(grid,titre=paste0(clade," (Rho : ",corTh," ",pCor,")"))

# coT = cor.test(grid$SES[which(grid$propExt!=0)],grid$propExt[which(grid$propExt!=0)])
# corTh = round(coT$estimate,3)
# pCor = ifelse(coT$p.value>0.05,"ns",ifelse(coT$p.value>0.01,"*",ifelse(coT$p.value>0.001,"**","***")))
# 
# png(paste0("Figures/Clade/",clade,"_SES.png"),width = 5000, height = 5000, units = "px", pointsize = 8,
#     bg = "white", res = 600, family = "")
# plot_erosion_SES(grid[which(grid$propExt!=0),],titre=paste0(clade," (Rho : ",corTh," ",pCor,")"))
# dev.off()   
#### 10. MATRIX SUMMARIZING RESULTS:
summary_var = function(variables = c("spRic", "spTraits", "spTraitsExt", "Density", "cell", "Long", "Lat"), dggs7 = dggs7, sitesdggs7 = sitesdggs7,
                        chge_mean_aux,CommsTPDExtinctions,CommsTPD,spClade,sesMeanClade){
  
  resultsdggs7 = matrix(NA, nrow=length(sitesdggs7), ncol=length(variables),
                        dimnames=list(names(sitesdggs7), variables))
  cellcenters   = dgSEQNUM_to_GEO(dggs7, as.numeric(names(sitesdggs7)))
  resultsdggs7[, "spRic"] = unlist(lapply(sitesdggs7,length))
  resultsdggs7[, "spTraits"] = rowSums(CommsTPD[,spClade])[names(sitesdggs7)]
  resultsdggs7[, "spTraitsExt"] = (rowSums(CommsTPD[,spClade])-rowSums(CommsTPDExtinctions[,spClade]))[names(sitesdggs7)]
  resultsdggs7[, "propExt"] = resultsdggs7[, "spTraitsExt"]/resultsdggs7[, "spTraits"]
  resultsdggs7[, "Density"] = chge_mean_aux
  #resultsdggs7[, "SES"] = sesMeanClade
  resultsdggs7[, "cell"] = as.numeric(names(sitesdggs7))
  resultsdggs7[, "Long"] = cellcenters$lon_deg
  resultsdggs7[, "Lat"] = cellcenters$lat_deg
  
  return(resultsdggs7)
}
### 11. Get the grid cell boundaries for cells which had speciess


