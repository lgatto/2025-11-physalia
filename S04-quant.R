library(QFeatures)
library(tidyverse)
library(PSMatch)

data(feat1)

feat1

colData(feat1)

colData(feat1)$X <- c("X1", "X2")

feat1$Y <- c("Y1", "Y2")

feat1

feat1[[1]]

feat1[["psms"]]

assay(feat1[[1]])

rowData(feat1[[1]])

feat1 <- aggregateFeatures(
    feat1,
    i = "psms",
    fcol = "Sequence",
    name = "peptides",
    fun = colMeans)

assay(feat1[[2]])

rowData(feat1[[2]])

## Ex: aggregate peptides into proteins

feat1 <- aggregateFeatures(
    feat1,
    i = "peptides",
    fcol = "Protein",
    name = "proteins",
    fun = colMedians)

feat1

assay(feat1[[3]])

rowData(feat1[[3]])

feat1["ProtA", , ]

filterFeatures(feat1, ~ pval < 0.05)

filterFeatures(feat1, ~ pval < 0.05)[[1]] |> rowData()

filterFeatures(feat1, ~ pval < 0.05, keep = TRUE)

filterFeatures(feat1, ~ location == "Mitochondrion")

## Ex: filter features and keep those that are not localised in the MT.

filterFeatures(feat1, ~ location != "Mitochondrion")

## Create a QFeatures object

data(hlpsms)

class(hlpsms)
dim(hlpsms)

hlpsms[1:5, 1:5]

hlpsms[1:2, ]

names(hlpsms)

hl <- readQFeatures(hlpsms, quantCols = 1:10, name = "psms")

hl

assay(hl[[1]])

rowData(hl[[1]])

se <- readSummarizedExperiment(hlpsms, quantCols = 1:10)

se

hl[[1]]


QFeatures(list(psms = se))

QFeatures(list(psms = se, peptides = se, proteins = se))

## CPTAC

f <- MsDataHub::cptac_a_b_peptides.txt()

basename(f)

## Ex: create a SE using readSummarizedExperiment() and/or a QFeatures object
## using readQFeatures().

head(readLines(f))

tab <- read.delim(f)

tab
dim(tab)
names(tab)

grep("Intensity\\.", names(tab), value = TRUE)

i <- grep("Intensity\\.", names(tab))

readQFeatures(tab, quantCols = i, name = "peptides")

## readSummarizedExperiment(f, quantCols = i,
##                          sep = "\t",
##                          fname = "Sequence")

se <- readSummarizedExperiment(tab, quantCols = i, fname = "Sequence")
se

## Use unique rownames!
## tmp <- readQFeatures(tab, quantCols = i, fname = "Proteins")

colnames(se) <- sub("Intensity\\.", "", colnames(se))

se

colData(se)

## Ex: Create an appropriate colData containing the 2 groups
## and the replicate numbers from the sample names.

colnames(se)


se$group <- rep(c("6A", "6B"), each = 3)
se$rep <- rep(7:9, times = 2)


## cd <- read.csv("....")
## identical(rownames(cd), colnames(se))
## colData(se) <- cd

se
colData(se)

names(rowData(se))

keep_var <- c("Sequence", "Proteins", "Leading.razor.protein", "PEP",
              "Score", "Reverse", "Potential.contaminant")

rowData(se) <- rowData(se)[, keep_var]

rowData(se)


## Missing values: random or not at random

head(assay(se))

se <- zeroIsNA(se)

names(nNA(se))

nNA(se)[[1]]

nNA(se)[["nNAcols"]]

barplot(nNA(se)[["nNAcols"]]$nNA)

nNA(se)[["nNArows"]]

table(nNA(se)[["nNArows"]]$nNA)

se <- filterNA(se, pNA = 4/6)

table(nNA(se)[["nNArows"]]$nNA)

MsCoreUtils::imputeMethods()

impute(se, method = "knn")

impute(se, method = "MinDet")

## Id QC/filtering

names(rowData(se))

table(rowData(se)$Reverse)

table(rowData(se)$Potential.contaminant)

## Ex: compare the Score distribution between Reverse and Fwd peptides

rowData(se) |>
    ggplot(aes(x = Score, colour = Reverse)) +
    geom_density()

rowData(se) |>
    ggplot(aes(x = PEP, colour = Reverse)) +
    geom_density()

prots <- rowData(se)$Proteins
names(prots) <- rowData(se)$Sequence
head(prots)
adj <- makeAdjacencyMatrix(prots, split = ";")
adj[1:10, 1:5]

## Create a QFeatures object

qf <- QFeatures(list(peptides = se))

qf

colData(qf)

colData(qf[[1]])

colData(qf) <- colData(se)

## Ex: Use filterFeatures() to remove reverse hits, contaminants and peptides
## with a PEP > 0.05.