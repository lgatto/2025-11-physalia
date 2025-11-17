## installation

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("remotes")
BiocManager::install("tidyverse")
BiocManager::install("factoextra")
BiocManager::install("MsDataHub")
BiocManager::install("mzR")
BiocManager::install("rhdf5")
BiocManager::install("rpx")
BiocManager::install("MsCoreUtils")
BiocManager::install("QFeatures")
BiocManager::install("Spectra")
BiocManager::install("ProtGenerics")
BiocManager::install("PSMatch")
BiocManager::install("pheatmap")
BiocManager::install("limma")
BiocManager::install("MSnID")
BiocManager::install("Biostrings")
BiocManager::install("cleaver")
BiocManager::install("RforMassSpectrometry/SpectraVis")

## version

version
BiocManager::version()

## access data

library("rpx")

px <- PXDataset("PXD000001")
px
pxtax(px)
pxref(px)
pxurl(px)

px

pxfiles(px)

print(object.size(px), units = "Mb")

fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"

mzf <- pxget(px, fn)

mzf

fa <- pxget(px, grep("fasta", pxfiles(px), value = TRUE))

fa

library("MsDataHub")

MsDataHub()
ko15.CDF()

MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()

## library("msdata")