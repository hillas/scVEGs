scVEGs <- function(data, pVal, pFlag, species, tmpName) {
  require(locfit)
  require(MASS)
  require(hydroGOF)
  require(calibrate)
  require(locfit)
  geneName <- rownames(data)
  infoFlag <- 1
  if (species == 'hsa') {
    geneLength <- read.delim('hg19_genes_length.tsv', header = TRUE, stringsAsFactors = TRUE)
  } else if (species == 'mmu') {
    geneLength <- read.delim('mm9_genes_length.tsv', header = TRUE, stringsAsFactors = TRUE)
  } else {
    print('no gene length information')
    infoFlag <- 0
  }
  # Calculate TPM
  if(infoFlag == 1) {
    rpk <- data
    gLength <- vector(mode="numeric", length=dim(data)[1])
    for (i in 1:dim(data)[1]) {
      dx <- which(geneLength$geneName == geneName[i])
      if(length(dx) > 0) {
        gLength[i] <- max(geneLength$geneLength[dx])
      } else {
        gLength[i] <- -1
      }
    }
    dx <- which(gLength == -1)
    rpk[dx, ] <- 0
    rpk <- rpk * 1000/kronecker(matrix(1, 1, dim(data)[2]), gLength)
    sscale <- apply(rpk, 2, sum)/10^6
    tpm <- rpk / kronecker(matrix(1, dim(data)[1], 1), sscale)
    flag <- matrix(0, dim(data)[1], dim(data)[2])
    flag[tpm >= 1] <- 1
    num <- apply(flag, 1, sum)
    ix <- which(num > 0.01*dim(data)[2])
    data <- data[ix, ]
    geneName <- geneName[ix]
  }

  m <- dim(data)[1]
  std <- apply(data, 1, sd)
  avg <- apply(data, 1, mean)
  cv <- std / avg
  # over dispersion sigma  (var = u(1 + u * sigma^2))
  rm(tpm)
  rm(flag)
  xdata <- (avg)
  ydata <- log10(cv)
  xdata <- xdata[is.na(ydata) != "TRUE"]
  ydata <- ydata[is.na(ydata) != "TRUE"]
  
  fitLoc <- locfit.robust(ydata ~ lp(log10(xdata), nn = .2))
  xSeq <- seq(min(log10(xdata)), max(log10(xdata)), 0.005)
  gapNum <- matrix(0, length(xSeq), 1)
  for(i in 1:length(xSeq)) {
    cdx <- which((log10(xdata) >= xSeq[i] - 0.05) & (log10(xdata) < xSeq[i] + 0.05))
    gapNum[i,1] <- length(cdx)
  }
  cdx <- which(gapNum > m*0.005)
  xSeq <- 10 ^ xSeq
  ySeq <- predict(fitLoc,log10(xSeq))
  yDiff <- diff(ySeq)
  ix <- which(yDiff > 0 & log10(xSeq[-1]) > 0)
  if(length(ix) == 0)
  ix <- length(ySeq) - 1
  xSeq_all <- 10^seq(min(log10(xdata)), max(log10(xdata)), 0.001)
  xSeq <- xSeq[cdx[1]:ix[1] + 1]
  ySeq <- ySeq[cdx[1]:ix[1] + 1]

  b <- 1
  a <- 0
  df <- data.frame(x=xSeq, y = ySeq)
  fit = nls(y ~ 0.5 * log10(b / x + a), data = df, 
            start=list(b = b,a = a), nls.control(maxiter = 500), na.action =  'na.exclude')
  newdf <- data.frame(x = xSeq_all)
  ydataFit <- predict(fit,newdata = newdf)
  
  # Calculate CV difference
  logX <- log10(xdata)

  logXseq <- log10(xSeq_all)
  cvDist <- matrix(0,length(xdata),1)
  for (i in 1:length(logX)) {
    cx <- which(logXseq >= logX[i] - 0.2 & logXseq < logX[i] + 0.2)
    tmp <- sqrt((logXseq[cx] - logX[i])^2 + (ydataFit[cx] - ydata[i])^2)
    tx <- which.min(tmp)
    
    if(logXseq[cx[tx]] > logX[i]) {
      if(ydataFit[cx[tx]] > ydata[i]) {
        cvDist[i] <- -1*tmp[tx]
      } else {
        cvDist[i] <- tmp[tx]
      }
      cvDist[i] <- -1*tmp[tx]
    } else if (logXseq[cx[tx]] <= logX[i]) {
      if(ydataFit[cx[tx]] < ydata[i]) {
        cvDist[i] <- tmp[tx]
      } else {
        cvDist[i] <- -1*tmp[tx]
      }
    } 
  }
  cvDist <- log2(10^cvDist)
  
  # use kernel density estimate to find the peak
  dor <- density(cvDist, kernel = "gaussian")
  distMid <-dor$x[which.max(dor$y)]
  dist2 <- cvDist - distMid
  tmpDist <- c(dist2[dist2 <= 0], abs(dist2[dist2 < 0])) + distMid
  distFit <- fitdistr(tmpDist, "normal")
  pRaw <- pnorm(cvDist, mean = distFit$estimate[1], sd = distFit$estimate[2], lower.tail = FALSE)
  pAdj <- p.adjust(pRaw, 'hochberg')
  if (pFlag == 1) {
    dx <- which(pAdj < pVal)
  } else {
    dx <- which(pRaw < pVal)
  }
  
  # Plot CV-mean figure with 2 fitted lines
  filename <- paste(tmpName,"cv",sep="_")
  pdf_filename<-paste(filename,".pdf",sep="")
  pdf(pdf_filename)
  # eps_filename <- paste(filename,".eps",sep="")
  # setEPS()
  # postscript(eps_filename)
  plot(log10(avg), log10(cv), type = 'p', pch = 20, cex = 0.5)
  abline(0,-0.5)
  lines(log10(xSeq), ySeq, col = 'cyan3')
  lines(log10(xSeq_all), ydataFit, col = 'firebrick')
  points(log10(xdata[dx]), ydata[dx], type = 'p', pch = 20, cex = 0.5, col = 'palegreen3')
  textxy(-1.5, min(log10(cv), na.rm = TRUE) + 0.3,
         paste('alpha = ', round(coef(fit)[2], digit = 3), sep=''), cex = 1, offset = 0)
  textxy(-1.5, min(log10(cv), na.rm = TRUE) + 0.2,
         paste('beta = ', round(coef(fit)[1], digit = 3), sep=''), cex = 1, offset = 0)
  textxy(-1.5, min(log10(cv), na.rm = TRUE) + 0.1,
         paste('variable genes = ', length(dx), sep=''), cex = 1, offset = 0)
  dev.off()
  # plot CV difference histogram and fitted normal distribution
  filename <- paste(tmpName,"cvHist",sep="_")
  pdf_filename<-paste(filename,".pdf",sep="")
  pdf(pdf_filename)
#   eps_filename <- paste(filename,".eps",sep="")
#   setEPS()
#   postscript(eps_filename)
  hist(cvDist, freq = FALSE, breaks = 200, xlim = c(-1, 2), ylim = c(0, max(dor$y) + 1), xlab = "CV difference")
  du <- distFit$estimate[1]
  ds <- distFit$estimate[2]
  # curve(dnorm(x, mean = du, sd = ds)*max(dor$y)/dnorm(du, mean = du, sd = ds)
  #       , col = "red", add = TRUE)
  lines(dor$x, dnorm(dor$x, mean = du, sd = ds)*max(dor$y)/dnorm(du, mean = du, sd = ds)
        , col = "red")
  lines(dor, col = "blue")
  dev.off()
  sig <- data[dx, ]
  # curve(dgamma(x, shape = distFit$estimate[1], rate = distFit$estimate[2]), col = "red", add = TRUE)
}