PBM.glmnet.getCoefs <- function(e, m, alpha=0.05, randomize=F, center=FALSE, file="AP2-coefs.txt") {
   e.coef <- apply(e, 2, IDC.glmnet, m, mode="coefficients", alpha=alpha, randomize=randomize)
   if (dim(e)[2] > 1) {
     e.coef.s <- t(apply(e.coef, 1, scale, center=center))
   } else {
     e.coef.s <- e.coef
   }
   rownames(e.coef.s) <- colnames(m)
   colnames(e.coef.s) <- colnames(e)
   write.table(e.coef.s, file=file, sep="\t", quote=F, row.names=T, col.names=NA)
   e.coef.s
}

PBM.glmnet.getFit <- function(e, m, alpha=0.05, file="fit.eps") {
   e.fits.true <- apply(e, 2, IDC.glmnet, m, mode="predict", alpha=alpha, randomize=F)
   e.fits.rand <- apply(e, 2, IDC.glmnet, m, mode="predict", alpha=alpha, randomize=T)
   ymax <- max(c(e.fits.true, e.fits.rand))
   print(sprintf("avg c=%f", mean( e.fits.true )))
   setEPS()
   postscript(file)
   barplot(rbind(e.fits.true, e.fits.rand), las=2, beside=T, col=c("red", "black"), ylim=c(0,ymax+0.15), ylab="Fit (max = 1.0)")
   mtext(paste0("avg c=", mean( e.fits.true )), side = 3)
   legend(10, ymax+0.15, c("Actual expression                                ", "Randomized expression"), c("red", "black"))
   par(cex=0.8)
   dev.off()
}

PBM.writeSingleMotifSingleConditionTargets <- function(coefs, e, m, motif, cond, t=2, outdir=".")
{
  c <- coefs[motif, cond]
  tgt <- c()
  if (c > 0) {
    tgt <- rownames(e)[ m[,motif]>0 & e[,cond]>=t ]
  } else {
    tgt <- rownames(e)[ m[,motif]<0 & e[,cond]<=t ]
  }
  write.table(tgt, file=sprintf("%s/%s.%s.txt", outdir, motif, cond), sep="\t", quote=F, row.names=F, col.names=F)
  tgt
}

PBM.writeTargetFiles <- function(coefs, e, m, fdr=0.01, corr=0.5, outdir=".")
{
  for (f in rownames(coefs)) {
    co <- t(apply(e[m[,f] > 0,], 1, function(x) { co <- cor.test(x, coefs[f,]); c(co$p.value, co$estimate) }))
    tgt <- co[p.adjust(co[,1], method="BH")<fdr & co[,2] > corr,]
    write.table(rownames(tgt), file=sprintf("%s/%s_targets.txt", outdir, f), sep="\t", quote=F, row.names=F, col.names=F)
  }

}

PBM.drawCoefHeatMap <- function(coefs, file="AP2-coefs.txt", minmax=2, clustrows=0, sortrowsbymax=0, sortrowsusingfile=NA)
{
   #todo <- sprintf("draw_expression_heatmap.pl --matrix=%s --draw=open --minmax=%d --saverownames=%s.rownames.txt", file, minmax, file)
   #if (clustrows == 1) {
     #todo <- paste(todo, "--clustrows=1")
   #}
   #if (sortrowsbymax == 1) {
     #todo <- paste(todo, "--sortrowsbymax=1")
   #}
   #if (!is.na(sortrowsusingfile)) {
     #todo <- paste(todo, sprintf("--sortrowsusingfile=%s", sortrowsusingfile))
   #}
   #print(todo)
   #system(todo)

   # make a better one..
   require(tidyr)
   require(cowplot)
   require(magrittr)
   notheme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.position="none",
                    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),plot.background=element_blank())

   color <- scale_fill_gradient2(low="blue", high="yellow", mid="black", midpoint=0)

   coefs <- read.table(file, header=T, stringsAsFactors=F)
   coefs$motif <- rownames(coefs)
   rownames(coefs) <- NULL
   coefs %<>% gather(tp, exp, -motif)
   o <- c("PF14_0533",
          "PFD0985wD1",
          "PF13_0026",
          "PF10_0075D2",
          "PFF0670wD2ext",
          "PF11_0442",
          "PF10_0075D3",
          "PF13_0235D1",
          "PFL1085w",
          "PF11_0404D1",
          "PF14_0633",
          "MALP8P1153",
          "PF14_0079",
          "PFD0985wD2",
          "PFF0200cDsL",
          "PF10_0075D1",
          "PF13_0267",
          "PF07_0126DsL",
          "PFL1900wDsL",
          "PF13_0097",
          "PFL1075w",
          "PFF0670wD1",
          "PF11_0091",
          "PFE0840cD2")
   coefs <- transform(coefs, motif = factor(motif, levels = o))

   g <- ggplot(coefs, aes(x = tp, y = motif)) +
     geom_tile(aes(fill=exp)) +
     scale_fill_gradient2(low="blue", mid="black", high="yellow", space="rgb", midpoint = 0, limits=c(-2.5,2.5)) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     xlab("Time") +
     ylab("Motif") +
     theme_classic() +
     ggtitle(paste0("ApiAP2 Activity"))

   ggsave(plot = g, file = paste0(file, ".png"))

}

PBM.plotAP2CoefVSExp <- function(coefs, e, m, fdr=0.01, corr=0.5, file="coef-vs-exp.pdf") {
  genes <- rownames(e)

  old_ids <-
    c("PF14_0533",
      "PFD0985wD1",
      "PF13_0026",
      "PF10_0075D2",
      "PFF0670wD2ext",
      "PF11_0442",
      "PF10_0075D3",
      "PF13_0235D1",
      "PFL1085w",
      "PF11_0404D1",
      "PF14_0633",
      "MALP8P1153",
      "PF14_0079",
      "PFD0985wD2",
      "PFF0200cDsL",
      "PF10_0075D1",
      "PF13_0267",
      "PF07_0126DsL",
      "PFL1900wDsL",
      "PF13_0097",
      "PFL1075w",
      "PFF0670wD1",
      "PF11_0091",
      "PFE0840cD2")

  new_ids <-
    c("PF3D7_1456000",
      "PF3D7_0420300",
      "PF3D7_1305200",
      "PF3D7_1007700",
      "PF3D7_0613800",
      "PF3D7_1143100",
      "PF3D7_1007700",
      "PF3D7_1342900",
      "PF3D7_1222600",
      "PF3D7_1139300",
      "PF3D7_1466400",
      "PF3D7_0802100",
      "PF3D7_1408200",
      "PF3D7_0420300",
      "PF3D7_0604100",
      "PF3D7_1007700",
      "PF3D7_1350900",
      "PF3D7_0730300",
      "PF3D7_1239200",
      "PF3D7_1317200",
      "PF3D7_1222400",
      "PF3D7_0613800",
      "PF3D7_1107800",
      "PF3D7_0516800")
  names(new_ids) <- old_ids

  par(mfrow=c(4,6))
  par(cex=0.5)
  pdf(file, onefile=TRUE)
  for (mo in rownames(coefs)) {
    ms <- new_ids[[mo]]
    AP2 <- !is.na(lapply(genes, function(x) { grep( sprintf("%s", x), ms) }) > 0)
    if (length(which(AP2 == T)) == 1) {
      cat(paste0(mo, "\n"))
      eAP2 <- e[AP2,]
      cAP2 <- coefs[mo,]

      # from mot file pick genes for a specific motif with a score greater than one
      # then grab the expression profiles of those genes and find the correlation between them
      # and the activity profile for the motif/protein being looked at
      co <- t(apply(e[m[,mo] >= 1,], 1, function(x) {co <- cor.test(x, coefs[mo,]); c(co$p.value, co$estimate) }))
      # filter by bonferroni corrected p-value (false discovery rate) and correlation value
      tgt <- co[p.adjust(co[,1], method="BH") < fdr & co[,2] > corr,]
      eTGT <- colMeans(e[ rownames(tgt),], na.rm=T)

      # if there wasn't a good enough fit for the particular motif skip it
      if (sum(is.na(eTGT)) > 0) {
        message(paste0(mo," was not a good fit; plot won't contain targets."))

        correl1 <- cor(as.numeric(eAP2),as.numeric(cAP2), use="complete.obs")
        correl2 <- NA
        eAP2 <- as.list(eAP2)
        cAP2 <- as.list(cAP2)
        #eTGT <- as.list(eTGT)
        # red is the AP2 protein expression profile
        plot(1:7, eAP2, type="l", col="red", xlab="", yaxt='n', ann=FALSE)
        par(new=T)
        # black is the AP2 activity profile -correlation values are with the black
        plot(1:7, cAP2, type="l", col="black", xlab="Time", yaxt='n', ann=FALSE)
        par(new=T)
        # green is average expression of predicted targets
        #plot(1:7, eTGT, type="l", col="green", main="", xlab="", yaxt='n', ann=FALSE)
        axis(2, at=seq(-3, 3, 0.1))
        title(main=sprintf("%s\nAP2 Exp vs. AP2 Act (red vs. black) = %2.1f\nAP2 Exp vs. Avg Tgt Exp (red vs. green) = %2.1f", mo, correl1, correl2))
        next()
      }

      correl1 <- cor(as.numeric(eAP2),as.numeric(cAP2), use="complete.obs")
      correl2 <- cor(as.numeric(eAP2),as.numeric(eTGT), use="complete.obs")
      eAP2 <- as.list(eAP2)
      cAP2 <- as.list(cAP2)
      eTGT <- as.list(eTGT)
      # red is the AP2 protein expression profile
      plot(1:7, eAP2, type="l", col="red", xlab="", yaxt='n', ann=FALSE)
      par(new=T)
      # black is the AP2 activity profile -correlation values are with the black
      plot(1:7, cAP2, type="l", col="black", xlab="Time", yaxt='n', ann=FALSE)
      par(new=T)
      # green is average expression of predicted targets
      plot(1:7, eTGT, type="l", col="green", main="", xlab="", yaxt='n', ann=FALSE)
      axis(2, at=seq(-3, 3, 0.1))
      title(main=sprintf("%s\nAP2 Exp vs. AP2 Act (red vs. black) = %2.1f\nAP2 Exp vs. Avg Tgt Exp (red vs. green) = %2.1f", mo, correl1, correl2))

    }
  }
  dev.off()
}

IDC.glmnet <- function(e, m, mode="coef", randomize=F, alpha=0.5) {
  nona  <- !is.na(e);
  enona <- e[nona]
  mnona <- m[nona,]
  e.cv <- cv.glmnet( mnona, enona, nfolds=10)
  l <- e.cv$lambda.min
  #print(l)
  if (randomize == TRUE) {
    enona <- sample(enona)
  }
  e.fits <- glmnet( mnona, enona, family="gaussian", alpha=alpha, nlambda=100)
  if (mode == "predict") {
        #cor.test ( predict(lars( mnona, enona, type="lasso"), mnona, type="fit", mode="fraction", s=s)$fit, enona)$estimate
     #cor.test ( predict(e.fits, mnona, type="response", s=l)$fit, enona)$estimate
       cor.test( predict(e.fits, mnona, type="response", s=l), enona)$estimate
  } else {
     as.matrix(predict(e.fits, s=l, type="coefficients")[-1,])

     #as.matrix(coef(lars( mnona, enona, type="lasso"), mode="fraction", s=s))
  }
}



# should be same thing as IDC glmnet just using a different package called lars; glmnet is more up to date apparently
# http://stats.stackexchange.com/questions/58531/using-lasso-from-lars-or-glmnet-package-in-r-for-variable-selection
IDC.lasso <- function(e, m, mode="coef", randomize=F) {
  nona  <- !is.na(e);
  enona <- e[nona]
  mnona <- m[nona,]
  e.cv <- cv.lars( mnona, enona, type="lasso", trace=F, K=10)
  s <- e.cv$fraction[order(e.cv$cv)[1]]
  if (randomize == TRUE) {
    enona <- sample(enona)
  }
  if (mode == "predict") {
        cor.test ( predict(lars( mnona, enona, type="lasso"), mnona, type="fit", mode="fraction", s=s)$fit, enona)$estimate
  } else {
     as.matrix(coef(lars( mnona, enona, type="lasso"), mode="fraction", s=s))
  }
}

