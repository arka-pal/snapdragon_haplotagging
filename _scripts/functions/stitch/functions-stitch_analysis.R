stitch_folder <- function(snpID, chrom, K, cov, niter, ngen, r){
  tmpFolder <- Sys.glob(file.path(paste0("???",chrom), 
                                  paste0(snpID, "-*_K", K, "_cov", cov, "_*_niter", niter,"_ngen", ngen, "_r", r)))
  return(tmpFolder)
}

#----#

stitch_Rdata <- function(snpID, chrom, K, cov, niter, ngen, r){
  tmpRData <- Sys.glob(file.path(stitch_folder(snpID, chrom, K, cov, niter, ngen, r), 
                                 "RData", 
                                 "EM.all*[^withBuffer].RData"))
  return(tmpRData)
}

#----#

extract_kasp_GT <- function(snpID, chrom, GT_kasp, coords_for_stitch){
  # Extract information on which SNP to use. 
  snp_rowID <- which(coords_for_stitch$CHROM == paste0("GWHBJVT0000000",chrom) & coords_for_stitch$snpID == snpID)
  LocusName <- coords_for_stitch[snp_rowID, "LocusName"]
  POS <- coords_for_stitch[snp_rowID, "POS"]
  REF <- coords_for_stitch[snp_rowID, "REF"]
  ALT <- coords_for_stitch[snp_rowID, "ALT"]
  
  # Choose which columns to keep from the whole KASP GT data table
  cols <- c("PlantID", coords_for_stitch$LocusName[snp_rowID])
  tmpGT <- GT_kasp[, ..cols]
  colnames(tmpGT) <- c("PlantID", "GT_kasp")
  tmpGT[, PlantID := str_to_upper(PlantID)]
  #tmpGT$GT_kasp <- as.factor(tmpGT$GT_kasp)
  
  # Replace nucloetides with 0,1,2
  tmpGT[GT_kasp == paste(REF, REF, sep = ":"), GT_kasp := "0"]
  tmpGT[GT_kasp == paste(REF, ALT, sep = ":"), GT_kasp := "1"]
  tmpGT[GT_kasp == paste(ALT, ALT, sep = ":"), GT_kasp := "2"]
  # everything else is set to -1
  tmpGT[GT_kasp %in% setdiff(levels(as.factor(tmpGT$GT_kasp)), as.character(0:2)), GT_kasp := -1]
  tmpGT[, GT_kasp := as.numeric(GT_kasp)]
  
  return(tmpGT)
}

#----#

extract_stitch_GT_highConf <- function(snpID, chrom, K, cov, niter, ngen, r){
  tmpGT <- fread(Sys.glob(file.path(stitch_folder(snpID, chrom, K, cov, niter, ngen, r), 
                                    paste0(snpID,"*",chrom,"*formatted.012.txt"))), header=F)
  colnames(tmpGT) <- c("PlantID", "GT_stitch_highConf")
  tmpGT[, PlantID := str_to_upper(PlantID)]
  tmpGT[, GT_stitch_highConf := as.numeric(GT_stitch_highConf)]
  return(tmpGT)
}

#----#

extract_stitch_GT_all <- function(snpID, chrom, K, cov, niter, ngen, r){
  tmpGT <- fread(Sys.glob(file.path(stitch_folder(snpID, chrom, K, cov, niter, ngen, r), "*formatted.PL.012.txt")), header=F)
  colnames(tmpGT) <- c("PlantID", "GT_stitch_all")
  tmpGT[, PlantID := str_to_upper(PlantID)]
  tmpGT[, GT_stitch_all := as.numeric(GT_stitch_all)]
  return(tmpGT)
}

#----#

merge_GT <- function(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch){
  tmpGT_kasp <- extract_kasp_GT(snpID, chrom, GT_kasp, coords_for_stitch)
  tmpGT_stitch_highConf <- extract_stitch_GT_highConf(snpID, chrom, K, cov, niter, ngen, r)
  tmpGT_stitch_all <- extract_stitch_GT_all(snpID, chrom, K, cov, niter, ngen, r)
  
  allGT <- reduce(list(tmpGT_kasp, tmpGT_stitch_highConf, tmpGT_stitch_all), 
                  merge, by = "PlantID", all = FALSE)
  return(allGT)
}

#----#

merge_GT_acrossK <- function(snpID, chrom, cov, K, niter, ngen, r, GT_kasp, coords_for_stitch){
  ## K: vector of all K values
  
  # Extract all KASP GTs
  allGT <- extract_kasp_GT(snpID, chrom, GT_kasp, coords_for_stitch)
  
  # Merge stitch GTs by each K
  for(tmpK in K){
    tmpGT <- merge_GT(snpID, chrom, tmpK, cov, niter, ngen, r, GT_kasp, coords_for_stitch)
    cols <- c("PlantID", "GT_stitch_highConf", "GT_stitch_all")
    tmpGT <- tmpGT[,..cols]
    colnames(tmpGT) <- c("PlantID", paste0("K",tmpK,"_stitchHC"), paste0("K",tmpK,"_stitch"))
    
    allGT <- merge(allGT, tmpGT, by = "PlantID")
  }
  
  return(allGT)
}

#----#

extract_crosstabs <- function(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch){
  allGT <- merge_GT(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch)
  
  crosstabs <- list()
  
  crosstabs[["highConf"]] <- table(allGT$GT_kasp, allGT$GT_stitch_highConf)
  crosstabs[["all"]] <- table(allGT$GT_kasp, allGT$GT_stitch_all)
  
  return(crosstabs)
}

#----#

extract_crosstab_counts <- function(crosstab){
  
  # check for missing genotypes in the crosstab and fill with 0. 
  for(gt in c("-1", "0", "1", "2")) {
    if (!any(rownames(crosstab) == gt)) {
      rowsOLD <- rownames(crosstab) 
      crosstab <- rbind(crosstab, rep(0, ncol(crosstab)))
      rownames(crosstab) <- c(rowsOLD, gt)
    }
    if (!any(colnames(crosstab) == gt)) {
      colsOLD <- colnames(crosstab)
      crosstab = cbind(crosstab, rep(0, nrow(crosstab)))
      colnames(crosstab) = c(colsOLD, gt)
    }
  }
  
  # Order the rows and coloumns in crosstab
  crosstab <- crosstab[c("-1", "0", "1", "2"), c("-1", "0", "1", "2")]
  
  # Extract genotype counts from the crosstab
  totalGTs <- sum(crosstab)  
  missingGTs_kasp <- sum(crosstab["-1",])
  missingGTs_stitch <- sum(crosstab[,"-1"])
  missingGTs <- missingGTs_kasp + missingGTs_stitch - crosstab["-1", "-1"]
  GT_match <- sum(diag(crosstab)) - crosstab["-1", "-1"]
  GT_mismatch_hom <- crosstab["0", "2"] + crosstab["2", "0"]
  GT_mismatch_het <- totalGTs - missingGTs - GT_match - GT_mismatch_hom
  
  return(data.frame(totalGTs, missingGTs_kasp, missingGTs_stitch, missingGTs, GT_match, GT_mismatch_hom, GT_mismatch_het))
}

#----# 

extract_stitch_metrics <- function(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch, highCov_samples){
  
  # Merge all genotypes
  allGT <- merge_GT(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch)
  # Load STITCH Rdata
  load(stitch_Rdata(snpID, chrom, K, cov, niter, ngen, r))
  
  # Extract SNP position
  snp_rowID <- which(coords_for_stitch$CHROM == paste0("GWHBJVT0000000",chrom) & coords_for_stitch$snpID == snpID)
  POS <- coords_for_stitch[snp_rowID, "POS"]
  MAF <- coords_for_stitch[snp_rowID, "MAF"]
  
  # Extract basic info
  markerID <- paste(chrom, snpID, sep = "_")
  info_pos <- info[which(L == as.numeric(POS))]
  info_mean <- mean(info)
  info_meanADJ <- mean(info[which(info!=1)])
  
  R2_alleleFreq <- cor(estimatedAlleleFrequency, alleleCount[,3])^2
  R2_dosageHC <- summary(lm(GT_kasp~GT_stitch_highConf, 
                            allGT[GT_stitch_highConf != -1 & GT_kasp != -1, ]))$r.squared
  R2_dosage <- summary(lm(GT_kasp~GT_stitch_all, 
                          allGT[GT_stitch_all != -1 & GT_kasp != -1, ]))$r.squared
  
  ## ALL samples
  
  # Extract crosstabs
  tabs <- extract_crosstabs(snpID, chrom, K, cov, niter, ngen, r, GT_kasp, coords_for_stitch)
  tab_counts <- extract_crosstab_counts(tabs$all)
  tab_countsHC <- extract_crosstab_counts(tabs$highConf)
  
  # Extract crosstab metrics
  noGTs_kasp <- tab_counts$missingGTs_kasp
  noGTs_stitch <- tab_counts$missingGTs_stitch
  noGTs <- tab_counts$missingGTs
  noGTsHC_kasp <- tab_countsHC$missingGTs_kasp
  noGTsHC_stitch <- tab_countsHC$missingGTs_stitch
  noGTsHC <- tab_countsHC$missingGTs
  
  sensHC <- (tab_countsHC$totalGTs - tab_countsHC$missingGTs_stitch)/tab_counts$totalGTs
  
  spec <- tab_counts$GT_match/(tab_counts$totalGTs - tab_counts$missingGTs)
  specHC <- tab_countsHC$GT_match/(tab_countsHC$totalGTs - tab_countsHC$missingGTs)
  HOM_mismatch <- tab_counts$GT_mismatch_hom/(tab_counts$totalGTs - tab_counts$missingGTs)
  HOM_mismatchHC <- tab_countsHC$GT_mismatch_hom/(tab_countsHC$totalGTs - tab_countsHC$missingGTs)
  HET_mismatch <- tab_counts$GT_mismatch_het/(tab_counts$totalGTs - tab_counts$missingGTs)
  HET_mismatchHC <- tab_countsHC$GT_mismatch_het/(tab_countsHC$totalGTs - tab_countsHC$missingGTs)
  
  
  ## High coverage samples
  highCovGT <- merge(allGT, highCov_samples, by = "PlantID")
  
  # Extract crosstabs
  tabs_highCov <- list()
  tabs_highCov$all <- table(highCovGT$GT_kasp, highCovGT$GT_stitch_all)
  tabs_highCov$highConf <- table(highCovGT$GT_kasp, highCovGT$GT_stitch_highConf)
  tab_highCov_counts <- extract_crosstab_counts(tabs_highCov$all)
  tab_highCov_countsHC <- extract_crosstab_counts(tabs_highCov$highConf)
  
  # Extract crosstab metrics
  highCov_noGTs_kasp <- tab_highCov_counts$missingGTs_kasp
  highCov_noGTs_stitch <- tab_highCov_counts$missingGTs_stitch
  highCov_noGTs <- tab_highCov_counts$missingGTs
  highCov_noGTsHC_kasp <- tab_highCov_countsHC$missingGTs_kasp
  highCov_noGTsHC_stitch <- tab_highCov_countsHC$missingGTs_stitch
  highCov_noGTsHC <- tab_highCov_countsHC$missingGTs
  
  highCov_sensHC <- (tab_highCov_countsHC$totalGTs - tab_highCov_countsHC$missingGTs_stitch)/tab_highCov_counts$totalGTs
  
  highCov_spec <- tab_highCov_counts$GT_match/(tab_highCov_counts$totalGTs - tab_highCov_counts$missingGTs)
  highCov_specHC <- tab_highCov_countsHC$GT_match/(tab_highCov_countsHC$totalGTs - tab_highCov_countsHC$missingGTs)
  highCov_HOM_mismatch <- tab_highCov_counts$GT_mismatch_hom/(tab_highCov_counts$totalGTs - tab_highCov_counts$missingGTs)
  highCov_HOM_mismatchHC <- tab_highCov_countsHC$GT_mismatch_hom/(tab_highCov_countsHC$totalGTs - tab_highCov_countsHC$missingGTs)
  highCov_HET_mismatch <- tab_highCov_counts$GT_mismatch_het/(tab_highCov_counts$totalGTs - tab_highCov_counts$missingGTs)
  highCov_HET_mismatchHC <- tab_highCov_countsHC$GT_mismatch_het/(tab_highCov_countsHC$totalGTs - tab_highCov_countsHC$missingGTs)
  
  
  stitch_metrics <- data.frame(chrom, #chromosome
                               snpID, #snpID
                               markerID, #chromosome_snpID
                               K, 
                               cov, 
                               niter,
                               ngen,
                               r,
                               
                               MAF,
                               
                               info_pos, #INFO score at SNP position
                               info_mean, #mean INFO score within 200kb region around SNP
                               info_meanADJ, #mean INFO adjusted after removing INFO==1 
                               
                               R2_alleleFreq, #mpileup vs STITCH allele frequency
                               R2_dosage, #KASP vs STITCH 012 genotype 
                               R2_dosageHC, #KASP vs high confidence STITCH 012 genotype 
                               # HC: high confidence
                               
                               noGTs_kasp,
                               noGTs_stitch,
                               noGTs, #no. of missing genotypes 
                               
                               noGTsHC_kasp,
                               noGTsHC_stitch,
                               noGTsHC,
                               
                               sensHC, #sensitivity=(non-missing STITCH genotypes/all samples) 
                               specHC, 
                               spec, #specificity=(genotype matches/non-missing genotypes)
                               HET_mismatchHC,
                               HET_mismatch,#heterozygote mismatch/non-missing genotypes
                               HOM_mismatchHC,
                               HOM_mismatch,#homozygote mismatch/non-missing genotypes
                               
                               highCov_noGTs_kasp,
                               highCov_noGTs_stitch,
                               highCov_noGTs,
                               
                               highCov_noGTsHC_kasp,
                               highCov_noGTsHC_stitch,
                               highCov_noGTsHC,
                               
                               highCov_sensHC,
                               highCov_specHC,
                               highCov_spec,
                               highCov_HET_mismatchHC,
                               highCov_HET_mismatch,
                               highCov_HOM_mismatchHC,
                               highCov_HOM_mismatch
  )
  
  return(stitch_metrics)
}

#----#

plot_theme <- theme_bw() + theme(panel.grid = element_blank())