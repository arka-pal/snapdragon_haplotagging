#### Functions to polarise alleles into ancestral/derived based on outgroup
#### written by Arka Pal
#### Date: 28.03.2023
#### last update: 28.03.2023

library(data.table)
library(dplyr)

#### Count alleles from poolSeq
count_poolAlleles <- function(alleleSequence, poolCounts, allele){
    return(as.integer(unlist(strsplit(poolCounts, split = ':'))[which(alleleSequence == allele)]))    
}

#### Check if more than ref/alt alleles are present
check_RefAlt <- function(alleleSequence, poolCounts, ref, alt){
    if(sum(as.integer(unlist(strsplit(poolCounts, split = ':'))[!alleleSequence %in% c(ref,alt)])) > 0){
        check = 1}
    else{check = 0}
    return(check)
}

#### Check PoolSeq sync file and POS file
merge_seqData <- function(syncFile, posFile, start, end){
    posFile_segment = posFile[V2 >= start & V2 <= end]
    merged <- merge(posFile_segment, syncFile[,c(2,4:6)], by='V2', all.x = TRUE)
    colnames(merged) <- c("pos", "chrom", "ref", "alt", "AN", "AC", "pop1", "pop2", "pop3")
    merged = merged %>% 
            mutate(Amajus_refFreq = 1-AC/AN) %>%
            mutate(Amajus_altFreq = AC/AN) %>%
            mutate(Amajus_majFreq = ifelse(Amajus_refFreq >= Amajus_altFreq, Amajus_refFreq, Amajus_altFreq)) %>%
            mutate(Amajus_majAllele = ifelse(Amajus_refFreq >= Amajus_altFreq, ref, alt)) %>%
            mutate(Amajus_minFreq = ifelse(Amajus_refFreq < Amajus_altFreq, Amajus_refFreq, Amajus_altFreq)) %>%
            mutate(Amajus_minAllele = ifelse(Amajus_refFreq < Amajus_altFreq, ref, alt))

    return(merged)
}


#### Compare allele data between PoolSeq and posFile 
compile_AlleleDat <- function(alleleSequence, alleleDat){

	# alleleSequence <- c('A','T','C','G','N','del')
	# alleleDat <- merged

	pool = c()
	pool_numAlleles = c()
	pool_totalCounts = c()
	pool_refCounts = c()
	pool_altCounts = c()
	pool_allele3Counts = c()
	pool_check = c()
	pool_allele3 = c()
	pool_alleleShared = c()
	pool_missing = c()


	for (site in 1:nrow(alleleDat)){
	    # cat(site, alleleDat$pos[site],'\n')
	    missing = 0
	    allele3 = 'NA'
	    alleleShared = 'NA'
	    
	    ref = alleleDat$ref[site]
	    alt = alleleDat$alt[site]
	    pop1 = alleleDat$pop1[site]
	    pop2 = alleleDat$pop2[site]
	    pop3 = alleleDat$pop3[site]

	    if(is.na(pop1) & is.na(pop2) & is.na(pop3)){missing = 1}
	    pop1 = ifelse(is.na(pop1), '0:0:0:0:0:0', pop1)
	    pop2 = ifelse(is.na(pop2), '0:0:0:0:0:0', pop2)
	    pop3 = ifelse(is.na(pop3), '0:0:0:0:0:0', pop3)
	    
	    popSum = as.integer(unlist(strsplit(pop1, split = ':'))) + 
	                as.integer(unlist(strsplit(pop2, split = ':'))) +
	                as.integer(unlist(strsplit(pop3, split = ':')))

	    tmp = paste(popSum, collapse = ':')
	    numAlleles = sum(popSum != 0)
	    totalCounts = sum(popSum)
	    refCounts = ifelse(!is.na(tmp), count_poolAlleles(alleleSequence, tmp, ref), 'NA')
	    altCounts = ifelse(!is.na(tmp), count_poolAlleles(alleleSequence, tmp, alt), 'NA')
	    allele3Counts = totalCounts - (refCounts + altCounts)
	    check = ifelse(!is.na(tmp), check_RefAlt(alleleSequence, tmp, ref, alt), 'NA')

	    if (check == 1){
	        allele3 = setdiff(alleleSequence[which(unlist(strsplit(tmp, split = ':')) > 0)], c(ref,alt))
	        alleleShared = intersect(alleleSequence[which(unlist(strsplit(tmp, split = ':')) > 0)], c(ref,alt))
	        if (length(allele3) != 1) {allele3 = paste(allele3, collapse = ':')}
	        if (length(alleleShared) == 0) {alleleShared = 'none'} else {alleleShared = paste(alleleShared, collapse = ':')}
	        #else if(length(alleleShared) > 1) {alleleShared = 'poly'}
	    }
	    
	    pool = c(pool, tmp)
	    pool_numAlleles = c(pool_numAlleles, numAlleles)
	    pool_totalCounts = c(pool_totalCounts, totalCounts)
	    pool_refCounts = c(pool_refCounts, refCounts) 
	    pool_altCounts = c(pool_altCounts, altCounts)
	    pool_allele3Counts = c(pool_allele3Counts, allele3Counts)
	    pool_check = c(pool_check, check)
	    pool_allele3 = c(pool_allele3, allele3)
	    pool_alleleShared = c(pool_alleleShared, alleleShared)
	    pool_missing = c(pool_missing, missing)
	}

	alleleDat$pool = pool
	alleleDat$pool_numAlleles = pool_numAlleles
	alleleDat$pool_totalCounts = pool_totalCounts
	alleleDat$pool_refCounts = pool_refCounts
	alleleDat$pool_altCounts = pool_altCounts
	alleleDat$pool_allele3Counts = pool_allele3Counts
	alleleDat$pool_check = pool_check
	alleleDat$pool_allele3 = pool_allele3
	alleleDat$pool_alleleShared = pool_alleleShared
	alleleDat$pool_missing = pool_missing

	alleleDat = alleleDat %>% rowwise() %>% 
				mutate(pool_refFreq = pool_refCounts/pool_totalCounts) %>%
                mutate(pool_altFreq = pool_altCounts/pool_totalCounts) %>%
                mutate(pool_allele3Freq = pool_allele3Counts/pool_totalCounts) %>%
                mutate(pool_minFreq = ifelse(pool_check == 1,
                                             min(pool_refFreq, pool_altFreq, pool_allele3Freq),
                                             min(pool_refFreq, pool_altFreq))
                      ) %>%
                mutate(pool_minAllele = ifelse(pool_check == 1, 
                                               case_when(
                                                   pool_minFreq == pool_refFreq ~ ref,
                                                   pool_minFreq == pool_altFreq ~ alt,
                                                   pool_minFreq == pool_allele3Freq ~ allele3),
                                               case_when(
                                                   pool_minFreq == pool_refFreq ~ ref,
                                                   pool_minFreq == pool_altFreq ~ alt))
                      ) %>% 
                mutate(check_minAllele = ifelse(pool_minAllele == Amajus_minAllele, 0, 1))

    return(alleleDat)
}


#### Process comparison for each chromosome
compile_AlleleDat_chromSegments <- function(alleleSequence, baseDIR, chrom, stitchRun){
  # alleleSequence =  c('A','T','C','G','N','del')
  # baseDIR = '~/snap_hap'
  # chrom='Chr6'
  # stitchRun = 'stitchRun1'
  
  filepath_chromSegment = file.path(baseDIR, 'ref_genome', 'chromSegments', paste0(chrom,'_segments.txt'))
  filepath_syncFile = file.path(baseDIR, 'Amolle_syncFiles', paste0(chrom,'.sync'))
  filepath_posFile = file.path(baseDIR,'Amajus_alleles-ancestral','pos', paste0(stitchRun, '_', chrom, '.final.pos'))

  # chromSegments = data.frame(V1 = seq(11001,20001,1000), V2 = seq(12000,21000,1000))  
  chromSegments = fread(filepath_chromSegment)
  syncFile = fread(filepath_syncFile)
  posFile = fread(filepath_posFile)
  
  alleleDat = foreach(segment = c(1:nrow(chromSegments)), .combine = rbind
                      ) %do% {
                        start = chromSegments$V1[segment]
                        end = chromSegments$V2[segment]
                        
                        tic()
                        alleleDat_raw = merge_seqData(syncFile, posFile, start, end)
                        alleleDat_processed = compile_AlleleDat(alleleSequence, alleleDat_raw)
                        cat(chrom, start, end, paste0('\t', nrow(alleleDat_processed), ' sites'), '\t')
                        toc()
                        
                        as.data.table(alleleDat_processed)
                      }
  
  return(alleleDat)
}

#### Assign alleles to ancestral/derived
assign_alleles_AncDer = function(alleleDat){
    # alleleDat is the dataset compiled by the function: compile_AlleleDat_chromSegments
    alleleDat = alleleDat %>% 

        mutate(

            ## Type of ancestral allele: major, minor, shared or unresolved 
            anc_alleleType = case_when(
                # Pattern A: 1 allele fixed in Amolle
                (pool_numAlleles == 1 & pool_check == 0 & check_minAllele == 0) ~ 'major',
                (pool_numAlleles == 1 & pool_check == 0 & check_minAllele == 1) ~ 'minor',

                # Pattern B: Both share both alleles, major alleles match
                (pool_numAlleles == 2 & pool_check == 0 & pool_minFreq == 0.5) ~ 'major',
                (pool_numAlleles == 2 & pool_check == 0 & pool_minFreq != 0.5 & check_minAllele == 0) ~ 'major',
                (pool_numAlleles > 2 & pool_check == 1 & (!pool_alleleShared %in% c('A','T','C','G','none',NA)) & pool_minFreq == 0.5) ~ 'major',
                (pool_numAlleles > 2 & pool_check == 1 & (!pool_alleleShared %in% c('A','T','C','G','none',NA)) & pool_minFreq != 0.5 & check_minAllele == 0) ~ 'major',
            
                # Pattern C: Both share at least 1 allele
                (pool_numAlleles == 2 & pool_check == 1 & pool_alleleShared != 'none') ~ 'shared',
                (pool_numAlleles > 2 & pool_check == 1 & pool_alleleShared %in% c('A','T','C','G')) ~ 'shared',

                # Pattern D-F
                TRUE ~ 'unresolved'
                ),

            ## Ancestral allele (for sites that are resolved)
            anc_allele_resolvedOnly = case_when(
                (anc_alleleType == 'major') ~ Amajus_majAllele,
                (anc_alleleType == 'minor') ~ Amajus_minAllele,
                (anc_alleleType == 'shared') ~ pool_alleleShared,
                (anc_alleleType == 'unresolved') ~ 'NA'
                ),
            
            ## Derived allele (for sites that are resolved)
            der_allele_resolvedOnly = case_when(
                (anc_allele_resolvedOnly == 'NA') ~ 'NA',
                TRUE ~ ifelse(anc_allele_resolvedOnly == ref, alt, ref)
                ),
            
            ## Ancestral allele (for all sites)
            anc_allele = case_when(
                (anc_alleleType == 'unresolved') ~ Amajus_majAllele,
                TRUE ~ anc_allele_resolvedOnly
                ),
            
            ## Derived allele (for all sites)
            der_allele = ifelse(anc_allele == ref, alt, ref)

            )

    return(alleleDat)
}