library(QFeatures)
library(tidyverse)
library(PSMatch)
library(limma)

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
## with a PEP > 0.05 in the QFeatures object.

qf <- qf |>
    filterFeatures(~ Reverse != "+") |>
    filterFeatures(~ Potential.contaminant != "+") |>
    filterFeatures(~ PEP < 0.05)

qf

## Transform and normalise

qf <- logTransform(qf, i = "peptides",
                   name = "log_peptides")

qf

qf <- normalize(qf, i = "log_peptides",
                name = "lognorm_peptides",
                method = "center.mean")

qf

plotDensities(assay(qf[[1]]))
plotDensities(assay(qf[[2]]))
plotDensities(assay(qf[[3]]))

## Ex: aggregate the date into proteins (using Leading.razor.protein) and
## median.

qf <- aggregateFeatures(qf,
                        i = "lognorm_peptides",
                        name = "proteins_med",
                        fcol = "Leading.razor.protein",
                        fun = colMedians,
                        na.rm = TRUE)

qf

table(rowData(qf[["proteins_med"]])$.n)

## PCA

library(factoextra)

qf[["proteins_med"]] |>
    filterNA() |>
    assay() |>
    t() |>
    prcomp() |>
    fviz_pca_ind(habillage = qf$group)

qf[[3]] |>
    filterNA() |>
    assay() |>
    t() |>
    prcomp() |>
    fviz_pca_ind(habillage = qf$group)

qf[["proteins_med"]] |>
    impute(method = "knn") |>
    assay() |>
    t() |>
    prcomp() |>
    fviz_pca_ind(habillage = qf$group)

## Visualisation

qf["P02787ups|TRFE_HUMAN_UPS", , 3:4] |>
    longForm() |>
    ggplot(aes(x = colname,
               y = value,
               colour = rowname)) +
    geom_line(aes(group = rowname)) +
    geom_point() +
    facet_wrap(~ assay)

plot(qf)

## CAREFUL!

##                              6A_7       6A_8      6A_9        6B_7        6B_8         6B_9
## APNHAVVTR                      NA         NA        NA  0.41079503 -2.12340198 -1.991160146
## ASYLDCIR                       NA         NA        NA -1.01813156 -0.30856517 -0.304856517
## FDEFFSEGCAPGSKK        -2.2898637 -1.8918720 -2.617213          NA           NA          NA
## IECVSAETTEDCIAK                NA         NA        NA -2.52842170 -1.79980353 -1.506790602
## KASYLDCIR                      NA         NA        NA -0.52486909 -1.45305770 -1.849688082
## SASDLTWDNLK            -1.3057464 -1.5910897 -1.574269          NA           NA          NA


## Ex: apply a different normalisation (quantiles.robust) on the log_peptides
## set and then aggregate that newly created set using colMedians().

qf |>
    normalize(i = "log_peptides",
              name = "logquantiles_peptides",
              method = "quantiles.robust") |>
    aggregateFeatures(i = "logquantiles_peptides",
                      name = "proteins_med2",
                      fcol = "Leading.razor.protein",
                      fun = colMedians,
                      na.rm = TRUE) |>
    plot()

## statistical analyses

pr <- getWithColData(qf, "proteins_med")

colData(pr)

design <- model.matrix(~ pr$group)
fit <- lmFit(assay(pr), design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "pr$group6B", number = Inf) |>
    rownames_to_column("protein") |>
    as_tibble() |>
    mutate(TP = grepl("ups", protein)) |>
    select(protein, logFC, adj.P.Val, TP)

res |>
    ggplot(aes(x = logFC,
               y = -log10(adj.P.Val),
               colour = TP)) +
    geom_point() +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = c(-1, 1))