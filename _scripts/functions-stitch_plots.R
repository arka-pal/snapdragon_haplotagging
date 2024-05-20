
args = commandArgs(trailingOnly=TRUE)

##usage: R_plots_stitch.R FOLDER_WITH_RESULTS chromosome

file = file.path(args[1], "RData", paste0("EM.all.", 'Chr6', ".RData"))
outputdir <- args[1]

load(file)

regionName <- args[2] ##!!!

## plot HWE, info, MAF, coloured by discrepency
plotMetricsForPostImputationQC <- function(
  iSample,
  highCovList,
  gen,
  gen_imp,
  alleleCount,
  chr,
  L,
  estimatedAlleleFrequency,
  info,
  outputdir,
  colour = TRUE,
  hwe,
  regionName
) {
  if (is.null(iSample)==FALSE) {
    m1A <- cbind(gen[,iSample] / 2, gen_imp[, iSample])
    m1A[alleleCount[,3]>0.5,] <- 1 - m1A[alleleCount[,3]>0.5,]
    dist <- abs(m1A[,1]-m1A[,2])
    dist[dist > 0.1] <- 0.1
  } else {
    dist <- rep(0, T)
  }
  colfunc <- colorRampPalette(c("blue", "red"))
  col=colfunc(11)[round(dist,2)*100+1]
  ## plot
  jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQC.",regionName,".sample",iSample,".jpg")),height=1000,width=3000,qual=100)
  par(mfrow=c(1,3))
  plot(info,log10(hwe) ,cex.lab=3 ,col=col)
  plot(info,estimatedAlleleFrequency,cex.lab=3,col=col)
  plot(estimatedAlleleFrequency,log10(hwe) ,cex.lab=3,col=col )
  dev.off()
  ## now plot along the chromosome as well
  jpeg(file.path(outputdir, "plots", paste0("metricsForPostImputationQCChromosomeWide.",regionName,".sample",iSample,".jpg")),height=3000,width=3000,qual=100)
  par(mfrow=c(3,1))
  plot(L,log10(hwe) ,cex.lab=3 ,col=col)
  plot(L,info,cex.lab=3,col=col)
  plot(L,estimatedAlleleFrequency,cex.lab=3,col=col )
  dev.off()
}

### plot estimated against actual allele frequency
plotEstimatedAgainstReal <- function(outputdir,alleleCount,estimatedAlleleFrequency,which,chr,regionName) {
  m1 <- cbind(alleleCount[, 3], estimatedAlleleFrequency)
  jpeg(file.path(outputdir, "plots", paste0("r2.",regionName,".goodonly.jpg")),height=800,width=800,qual=100)
  m1 <- m1[which, ]
  corr <- suppressWarnings(
    round(cor(m1[,1], m1[,2], use="complete.obs") ** 2, 3)
  )
  main <- paste("nSnps: ",dim(m1)[1],", r2:",corr,sep="")
  plot(x=m1[,1],y=m1[,2],xlim=c(0,1),ylim=c(0,1),main=main,xlab="Real Allele Frequency",ylab="Estimated allele frequency",pch=20,cex.main=1.5,cex.lab=1.5)
  abline(0,1)
  dev.off()
}

## make interim plots of things that might be useful to understand performance
## hapSumCurrent_t, regular and log10
## alphaMatCurrent_t
interim_plotter <- function(
  outputdir,
  regionName,
  iteration,
  L_grid,
  hapSumCurrent_tc,
  alphaMatCurrent_tc,
  sigmaCurrent_m,
  N,
  final_iteration = FALSE,
  is_reference = FALSE
) {
  ## easiest - just break up by S
  nGrids <- ncol(alphaMatCurrent_tc) + 1
  K <- nrow(alphaMatCurrent_tc)
  S <- dim(hapSumCurrent_tc)[3]
  for(s in 1:S) {
    alphaMatCurrent_t <- array(0, c(K, nGrids - 1))
    alphaMatCurrent_t[] <- alphaMatCurrent_tc[, , s]
    hapSumCurrent_t <- array(0, c(K, nGrids))
    hapSumCurrent_t[] <- hapSumCurrent_tc[, , s]
    plotHapSumCurrent_t(
      L_grid = L_grid,
      K = K,
      hapSumCurrent_t = hapSumCurrent_t,
      nGrids = nGrids,
      N = N,
      outputdir = outputdir,
      iteration = iteration,
      regionName = regionName,
      s = s,
      S = S,
      final_iteration = final_iteration,
      is_reference = is_reference
    )
    plotHapSumCurrent_t_log(
      L_grid = L_grid,
      K = K,
      hapSumCurrent_t = hapSumCurrent_t,
      nGrids = nGrids,
      N = N,
      outputdir = outputdir,
      regionName = regionName,
      iteration = iteration,
      s = s,
      S = S,
      final_iteration = final_iteration,
      is_reference = is_reference            
    )
    plotAlphaMatCurrent_t(
      L_grid = L_grid,
      alphaMatCurrent_t = alphaMatCurrent_t,
      sigmaCurrent = sigmaCurrent_m[, s],
      outputdir = outputdir,
      iteration = iteration,
      regionName = regionName,
      s = s,
      S = S,
      final_iteration = final_iteration,
      is_reference = is_reference            
    )
  }
  return(NULL)
}

## plot hapSumCurrent along genome
plotHapSumCurrent_t <- function(
  outname,
  L_grid,
  K,
  hapSumCurrent_t,
  nGrids,
  N,
  outputdir,
  regionName,
  iteration,
  s,
  S,
  final_iteration,
  is_reference
) {
  ##
  outname <- interim_plot_name("hapSum", regionName, outputdir, s, S, iteration, final_iteration, is_reference)
  ## 
  width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)    
  png(outname, height = 10, width = width, res = 100, units = "in")    
  colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  nCols <- length(colStore)
  sum <- array(0, nGrids)
  xlim <- range(L_grid)
  ylim <- c(0, 1)
  ## OK so if there are grids, use the grid points
  plot(x = L_grid[1], y = 0, xlim = xlim, ylim = ylim, axes = FALSE)
  x <- c(L_grid[1], L_grid, L_grid[length(L_grid):1])
  m <- array(0, c(nGrids, K + 1))
  for(i in 1:K) {
    m[, i + 1] <- m[, i] + hapSumCurrent_t[i, , drop = FALSE] / N
  }
  for(i in K:1) {
    polygon(
      x = x, y = c(m[1, i], m[, i + 1], m[nGrids:1, i]),
      xlim = xlim, ylim = ylim, col = colStore[(i %% nCols) + 1]
    )
  }
  dev.off()
}

## plot hapSumCurrent along genome
plotHapSumCurrent_t_log <- function(
  L_grid,
  K,
  hapSumCurrent_t,
  nGrids,
  N,
  outputdir,
  regionName,
  iteration,
  s,
  S,
  final_iteration,
  is_reference
) {
  outname <- interim_plot_name("hapSum_log", regionName, outputdir, s, S, iteration, final_iteration, is_reference)
  ## 
  width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)
  png(outname, height = 10, width = width, res = 100, units = "in")
  main <- "log10 of average haplotype usage vs physical position"
  colStore <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  nCols <- length(colStore)
  sum <- array(0, nGrids)
  xlim <- range(L_grid) / 1e6    
  ylim <- c(log10(max(1, min(hapSumCurrent_t))), log10(max(hapSumCurrent_t)))
  if (sum(hapSumCurrent_t) == 0) {
    stop("Something has done wrong and an ampty hapSumCurrent_t has been passed to plotHapSumCurrent_t_log")
  }
  plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, main = main, xlab = "Physical position (Mbp)", ylab = "log10 Haplotype usage")
  axis(1)
  axis(2)
  x <- L_grid / 1e6
  for(k in 1:K) {
    points(x = x, y = log10(hapSumCurrent_t[k, ]), col = colStore[(k %% nCols) + 1], type = "l")
  }
  for(i in floor(ylim[1]):ceiling(ylim[2])) {
    abline(h = i, col = "grey")
  }
  abline(h = log10(N / K), col = "red")
  dev.off()
  ##    system(paste0("rsync -av '", outname, "' florence:~/"))
}

plotAlphaMatCurrent_t <- function(L_grid, alphaMatCurrent_t, sigmaCurrent, outputdir, iteration, regionName, s, S, final_iteration, is_reference) {
  if (S > 1) {
    suffix <- ".png"
  } else {
    suffix <- paste0(".s.", s, ".png")
  }
  nGrids <- ncol(alphaMatCurrent_t)
  K <- nrow(alphaMatCurrent_t)
  ## two views - proportional, total?
  xleft <- (L_grid)[-length(L_grid)]
  xright <- (L_grid)[-1]
  xlim <- range(L_grid)
  ##
  for(i_what in 1:2) {
    ## make barplot type thing
    if (i_what == 2) {
      ## normalize by sigmaCurrent
      main <- "alphaMatCurrent_t P(q_t = k, I_t = 1)"
      x <- alphaMatCurrent_t
      for(k in 1:K) {
        x[k, ] <- alphaMatCurrent_t[k, ] * (1 - sigmaCurrent)
      }
      fbd_store <- list(list(gammaK_t = x))
      what <- "normalized"
    } else {
      main <- "alphaMatCurrent_t P(q_t | I_t = 1)"
      fbd_store <- list(list(gammaK_t = alphaMatCurrent_t))
      what <- "all"
    }
    outname <- interim_plot_name("alphaMat", regionName, outputdir, s, S, iteration, final_iteration, is_reference, what)
    width <- min(max(20, (L_grid[length(L_grid)] - L_grid[1]) / 1e6 * 12), 200)
    png(outname, height = 10, width = width, res = 100, units = "in")
    plot_fbd_store(fbd_store, xleft, xright, xlim, main = main)
    dev.off()
  }
}

interim_plot_name <- function(name, regionName, outputdir, s, S, iteration, final_iteration, is_reference, what = "") {
  if (final_iteration) {
    it_name <- ""
  } else {
    it_name <- paste0(".iteration.", iteration)
  }
  if (S > 1) {
    suffix <- ".png"
  } else {
    suffix <- paste0(".s.", s, ".png")
  }
  if (is_reference) {
    name <- paste0("ref.", name)
  }
  if (nchar(what) > 0) {
    what <- paste0(".", what)
  }
  outname <- file.path(outputdir, "plots", paste0(name, ".", regionName, it_name, what, suffix))
  return(outname)
}

plot_fbd_store <- function(fbd_store, xleft, xright, xlim, main = "Haplotype usage per-sample", mbp_xlab = FALSE, xleft2 = NULL, xright2 = NULL) {
  if (is.null(xleft2) == FALSE) {
    xleft <- xleft2
    xright <- xright2
  }
  if (mbp_xlab) {
    xlab <- "Physical position (Mbp)"
    xleft <- xleft / 1e6
    xright <- xright / 1e6
    xlim <- xlim / 1e6
    d <- 1e1 / 1e6
  } else {
    xlab <- "Physical position"
    d <- 0
  }
  NN <- length(fbd_store)
  ylim <- c(1, NN + 1)
  plot(x = 0, y = 0, ylim = ylim, xlim = xlim, xlab = xlab, ylab = "Sample", main = main, col = "white")
  cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)
  nr <- length(xleft)
  for(iSample in 1:NN) {
    ## plot each of them, as rectangle?
    R <- fbd_store[[iSample]]$gammaK_t
    ## should sum to 1, plot rectangles of each
    ybottom <- iSample + array(0, nr)
    ytop <- iSample + array(0, nr)
    for(k in 1:nrow(R)) {
      if (ncol(R) != nr) {
        ytop <- ytop + R[k, -ncol(R)]
      } else {
        ytop <- ytop + R[k, ]
      }
      rect(
        xleft = xleft - d,
        xright = xright + d, ## should not be necessary, R artefact I think
        ybottom = ybottom,
        ytop = ytop,
        col = cbPalette[k],
        border = NA
      )
      ybottom <- ytop
    }
  }
}

alpha_col <- function(col, alpha) {
  x <- col2rgb(col) / 255
  return(rgb(x["red", 1], x["green", 1], x["blue", 1], alpha = alpha)    )
}

plotMetricsForPostImputationQC(
  iSample = NULL, highCovList = NULL, gen = gen, gen_imp = gen_imp,
  alleleCount = alleleCount, chr = chr, L = L,
  estimatedAlleleFrequency = estimatedAlleleFrequency, info = info,
  outputdir = outputdir, hwe = hwe, regionName = regionName
)
##
## plot estimated AF against real 0.1X pileups (get r2 as well)
##
if(sum(passQC) > 1) {
  #print_message("Make estimated against real")
  plotEstimatedAgainstReal(
    outputdir = outputdir,alleleCount=alleleCount,
    estimatedAlleleFrequency=estimatedAlleleFrequency,
    which=passQC,chr=chr,regionName=regionName
  )
}
##
## plot hapProbs, alphaMat, etc
##
#print_message("Make other plots")
interim_plotter(
  outputdir = outputdir,
  regionName = regionName,
  iteration = NA,
  L_grid = L_grid,
  hapSumCurrent_tc = hapSumCurrent_tc,
  alphaMatCurrent_tc = alphaMatCurrent_tc,
  sigmaCurrent_m = sigmaCurrent_m,
  N = N,
  final_iteration = TRUE
)