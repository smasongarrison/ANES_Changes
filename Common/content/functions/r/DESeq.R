### R code from vignette source 'vignettes/DESeq/inst/doc/DESeq.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(digits=3)


###################################################
### code chunk number 2: pasilla2
###################################################
library( "DESeq" )
library( "pasilla" )
data( "pasillaGenes" )


###################################################
### code chunk number 3: look_into_pasillaGenes
###################################################
head( counts(pasillaGenes) )
pData( pasillaGenes )


###################################################
### code chunk number 4: pairedSamples
###################################################
pairedSamples <- pData(pasillaGenes)$type == "paired-end"
countsTable <- counts(pasillaGenes)[ , pairedSamples ]
conds <- pData(pasillaGenes)$condition[ pairedSamples ]


###################################################
### code chunk number 5: conds
###################################################
conds


###################################################
### code chunk number 6: conds (eval = FALSE)
###################################################
## #not run
## conds <- factor( c( "treated", "treated", "untreated", "untreated" ) )


###################################################
### code chunk number 7: instantiate
###################################################
cds <- newCountDataSet( countsTable, conds )


###################################################
### code chunk number 8: headcounts1
###################################################
head( counts(cds) )


###################################################
### code chunk number 9: estimateSizeFactors
###################################################
cds <- estimateSizeFactors( cds )
sizeFactors( cds )


###################################################
### code chunk number 10: headcounts2
###################################################
head( counts( cds, normalized=TRUE ) )


###################################################
### code chunk number 11: estimateDispersions
###################################################
cds <- estimateDispersions( cds )


###################################################
### code chunk number 12: str
###################################################
str( fitInfo(cds) )


###################################################
### code chunk number 13: fitUntr
###################################################
plotDispEsts <- function( cds )
{
   plot(
      rowMeans( counts( cds, normalized=TRUE ) ),
      fitInfo(cds)$perGeneDispEsts,
      pch = '.', log="xy" )
   xg <- 10^seq( -.5, 5, length.out=300 )
   lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}


###################################################
### code chunk number 14: figFit
###################################################
plotDispEsts( cds )


###################################################
### code chunk number 15: DESeq.Rnw:276-277
###################################################
all(table(conditions(cds))==2)


###################################################
### code chunk number 16: head
###################################################
head( fData(cds) )


###################################################
### code chunk number 17: str
###################################################
str( fitInfo(cds) )


###################################################
### code chunk number 18: nbt1
###################################################
res <- nbinomTest( cds, "untreated", "treated" )


###################################################
### code chunk number 19: nbt2
###################################################
head(res)


###################################################
### code chunk number 20: figDE
###################################################
plotDE <- function( res )
   plot(
      res$baseMean,
      res$log2FoldChange,
      log="x", pch=20, cex=.3,
      col = ifelse( res$padj < .1, "red", "black" ) )

plotDE( res )


###################################################
### code chunk number 21: histp
###################################################
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")


###################################################
### code chunk number 22: ressig1
###################################################
resSig <- res[ res$padj < 0.1, ]


###################################################
### code chunk number 23: ressig2
###################################################
head( resSig[ order(resSig$pval), ] )


###################################################
### code chunk number 24: ressig3
###################################################
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )


###################################################
### code chunk number 25: ressig4
###################################################
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )


###################################################
### code chunk number 26: writetable (eval = FALSE)
###################################################
## #not run
## write.table( res, file="results.txt" )


###################################################
### code chunk number 27: ncu
###################################################
ncu <- counts( cds, normalized=TRUE )[ , conditions(cds)=="untreated" ]


###################################################
### code chunk number 28: MArepl
###################################################
plot( rowMeans(ncu), log2( ncu[,2] / ncu[,1] ), pch=".", log="x" )


###################################################
### code chunk number 29: subset
###################################################
cdsTTU <- cds[ , 1:3]
pData( cdsTTU )


###################################################
### code chunk number 30: est123
###################################################
cdsTTU <- estimateSizeFactors( cdsTTU )
cdsTTU <- estimateDispersions( cdsTTU )
resTTU <- nbinomTest( cdsTTU, "untreated", "treated" )


###################################################
### code chunk number 31: figDE_Tb
###################################################
plotDE( resTTU )


###################################################
### code chunk number 32: subset2
###################################################
cds2 <- cds[ ,c(  "untreated3fb", "treated3fb"   ) ]


###################################################
### code chunk number 33: cds2
###################################################
cds2 <- estimateDispersions( cds2, method="blind", sharingMode="fit-only" )


###################################################
### code chunk number 34: res2
###################################################
res2 <- nbinomTest( cds2, "untreated", "treated" )


###################################################
### code chunk number 35: figDE2
###################################################
plotDE( res2 )


###################################################
### code chunk number 36: addmarg
###################################################
addmargins( table( res_sig = res$padj < .1, res2_sig = res2$padj < .1 ) )


###################################################
### code chunk number 37: design
###################################################
design <- pData(pasillaGenes)[ , c("condition","type") ]
design


###################################################
### code chunk number 38: fct
###################################################
fullCountsTable <- counts( pasillaGenes )

cdsFull <- newCountDataSet( fullCountsTable, design )


###################################################
### code chunk number 39: estsfdisp
###################################################
cdsFull <- estimateSizeFactors( cdsFull )
cdsFull <- estimateDispersions( cdsFull )


###################################################
### code chunk number 40: figFitPooled
###################################################
plotDispEsts( cdsFull )


###################################################
### code chunk number 41: fit1
###################################################
fit1 <- fitNbinomGLMs( cdsFull, count ~ type + condition )
fit0 <- fitNbinomGLMs( cdsFull, count ~ type  )


###################################################
### code chunk number 42: fitstr
###################################################
str(fit1)


###################################################
### code chunk number 43: pvalsGLM
###################################################
pvalsGLM <- nbinomGLMTest( fit1, fit0 )
padjGLM <- p.adjust( pvalsGLM, method="BH" )


###################################################
### code chunk number 44: addmarg2
###################################################
tab = table( "paired end only" = res$padj < .1, "all samples" = padjGLM < .1 )
addmargins( tab )


###################################################
### code chunk number 45: figPval
###################################################
bottom = function(x, theta=1e-12) pmax(x, theta)
plot( bottom(res$pval), bottom(pvalsGLM), log="xy", pch=20, cex=.3 )
abline(a=0, b=1, col="blue")


###################################################
### code chunk number 46: lookatfit1
###################################################
head(fit1)


###################################################
### code chunk number 47: fullAnalysisSimple
###################################################
cdsFullB <- newCountDataSet( fullCountsTable, design$condition )
cdsFullB <- estimateSizeFactors( cdsFullB )
cdsFullB <- estimateDispersions( cdsFullB )
resFullB <- nbinomTest( cdsFullB, "untreated", "treated" )


###################################################
### code chunk number 48: table
###################################################
addmargins(table(
   `all samples simple` = resFullB$padj < 0.1,
   `all samples GLM`    = padjGLM < 0.1 ))


###################################################
### code chunk number 49: rs
###################################################
rs <- rowSums ( counts ( cdsFull ))
use <- (rs > quantile(rs, 0.4))
table(use)
cdsFilt <- cdsFull[ use, ]


###################################################
### code chunk number 50: check
###################################################
stopifnot(!any(is.na(use)))


###################################################
### code chunk number 51: fitFilt
###################################################
fitFilt1  <- fitNbinomGLMs( cdsFilt, count ~ type + condition )
fitFilt0  <- fitNbinomGLMs( cdsFilt, count ~ type  )
pvalsFilt <- nbinomGLMTest( fitFilt1, fitFilt0 )
padjFilt  <- p.adjust(pvalsFilt, method="BH" )


###################################################
### code chunk number 52: doublecheck
###################################################
stopifnot(all.equal(pvalsFilt, pvalsGLM[use]))


###################################################
### code chunk number 53: tab
###################################################
padjFiltForComparison = rep(+Inf, length(padjGLM))
padjFiltForComparison[use] = padjFilt
tab = table( `no filtering`   = padjGLM < .1,
             `with filtering` = padjFiltForComparison < .1 )
addmargins(tab)


###################################################
### code chunk number 54: figscatterindepfilt
###################################################
plot(rank(rs)/length(rs), -log10(pvalsGLM), pch=".")


###################################################
### code chunk number 55: histindepfilt
###################################################
h1 = hist(pvalsGLM[!use], breaks=50, plot=FALSE)
h2 = hist(pvalsGLM[use], breaks=50, plot=FALSE)
colori = c(`do not pass`="khaki", `pass`="powderblue")


###################################################
### code chunk number 56: fighistindepfilt
###################################################
barplot(height = rbind(h1$counts, h2$counts),
        beside = FALSE, col = colori,
        space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0,
     label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


###################################################
### code chunk number 57: vsd
###################################################
cdsBlind <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cdsBlind )


###################################################
### code chunk number 58: modlr
###################################################
mod_lfc <- (rowMeans( vsd[, conditions(cds)=="treated", drop=FALSE] ) -
            rowMeans( vsd[, conditions(cds)=="untreated", drop=FALSE] ))


###################################################
### code chunk number 59: dah
###################################################
lfc <- res$log2FoldChange
finite <- is.finite(lfc)
table(as.character(lfc[!finite]), useNA="always")


###################################################
### code chunk number 60: repl
###################################################
largeNumber <- 10
lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)


###################################################
### code chunk number 61: figmodlr
###################################################
plot( lfc, mod_lfc, pch=20, cex=.3,
      col = ifelse( finite, "#80808040", "red" ) )
abline( a=0, b=1, col="#40404040" )


###################################################
### code chunk number 62: figHeatmap2a
###################################################
select <- order(res$pval)[1:40]
colors <- colorRampPalette(c("white","darkblue"))(100)
heatmap( vsd[select,],
         col = colors, scale = "none")


###################################################
### code chunk number 63: figHeatmap2b
###################################################
heatmap( counts(cds)[select,],
         col = colors, scale = "none")


###################################################
### code chunk number 64: sampleClust
###################################################
cdsFullBlind <- estimateDispersions( cdsFull, method = "blind" )
vsdFull <- getVarianceStabilizedData( cdsFullBlind )
dists <- dist( t( vsdFull ) )


###################################################
### code chunk number 65: figHeatmapSamples
###################################################
heatmap( as.matrix( dists ),
   symm=TRUE, scale="none", margins=c(10,10),
   col = colorRampPalette(c("darkblue","white"))(100),
   labRow = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type ) )


###################################################
### code chunk number 66: sessi
###################################################
sessionInfo()


