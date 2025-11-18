library(Biostrings)
library(cleaver)
library(MsDataHub)
library(PSMatch)
library(ggplot2)
library(tidyverse)

packageVersion("MsDataHub")
BiocManager::version()
version

## identifcations

idf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
idf
id <- PSM(idf)

class(id)

id
dim(id)

names(id)

cnt <- readAAStringSet(MsDataHub::crap_gpm.fasta())
cnt
rev(cnt)
pep <- cleave(cnt)
cleave(rev(cnt))

id$sequence

## Ex: How many decoy hits?

table(id$isDecoy)

## Ex: Verify that this table contains 5802 matches for 5343 scans and 4938
## peptides.

nrow(id)

length(unique(id$spectrumID))

table(id$rank)

length(unique(id$sequence))

table(table(id$spectrumID))


i <- which(id$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

data.frame(id[i, ])[, 1:5]

id2 <- reducePSMs(id, id$spectrumID)

id2

j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

id2[j, ]

id2[j, "spectrumID"]

id2[j, "DatabaseAccess"]

## Filtering id results
## - keep rank 1

table(id$rank == 1)

## - remove decoys

table(id$isDecoy)

dim(id)

describeProteins(id)

describePeptides(id)

id_filtered <- filterPSMs(id)

describePeptides(id_filtered)

## Ex: compare the distribution of the MS.GF.RawScore between the decoy and
## non-decoy/forward matches.

ggplot(id,
       aes(x = MS.GF.RawScore,
           colour = isDecoy)) +
    geom_density()

## Peptide-protein relations


id2 <- id |>
    filterPsmDecoy() |>
    filterPsmRank()

id2

data.frame(id2[1:10, c("sequence", "DatabaseAccess")])

adj <- makeAdjacencyMatrix(id2)

dim(adj)

adj[1:5, 1:5]

cc <- ConnectedComponents(adj)

cc

connectedComponents(cc, 1)

connectedComponents(cc, 527)

connectedComponents(cc, 38)

connectedComponents(cc, 920)

length(cc)

(i <- which(nrows(cc) > 2 & ncols(cc) > 2))

dims(cc)[i, ]

cx <- connectedComponents(cc, 1082)

cx

plotAdjacencyMatrix(cx)

plotAdjacencyMatrix(cx, 1)

a2 <- adj[1:10, c("ECA3389", "ECA3399")]

cc

## Adding id to raw data

id_filtered

sp

table(table(id_filtered$spectrumID))

id4 <- id_filtered[which(id_filtered$spectrumID == "controllerType=0 controllerNumber=1 scan=5490"), ]

data.frame(id4[, c("spectrumID", "sequence", "modName", "modLocation")])

id_filtered <- reducePSMs(id_filtered, id_filtered$spectrumID)

spectraVariables(sp)

sp <- joinSpectraData(sp, id_filtered,
                      by.x = "spectrumId",
                      by.y = "spectrumID")

sp

spectraVariables(sp)

## Ex: verify that the identiciation data has been added to the correct spectra.

## - check that MS1 scans have no id info.

all(is.na(filterMsLevel(sp, 1)$sequence))

## - check how many MS2 scans have an id.

table(!is.na(filterMsLevel(sp, 2)$sequence))

## - Compare the precursor MZ to the peptide/theoretical MZ

sp2 <- filterMsLevel(sp, 2)

summary(sp2$precursorMz - sp2$experimentalMassToCharge)

sp2$sequence[407]

sp2$precursorMz[407]

## id-annotated chromatogram

filterMsLevel(sp, 1) |>
    spectraData() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()

sp <- countIdentifications(sp)

table(sp$msLevel)

table(msLevel(sp), sp$countIdentifications)

filterMsLevel(sp, 1) |>
    spectraData() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line() +
    geom_point(
        aes(colour = ifelse(countIdentifications == 0,
                            NA, countIdentifications)),
        size = 3)


table(sp2$isDecoy)

## Visualising PSMs

i <- which(sp$MS.GF.RawScore > 100)[1]

plotSpectra(sp[i])

calculateFragments(sp[i]$sequence)

labelFragments(sp[i])

plotSpectra(sp[i], labels = labelFragments,
            labelPos = 3, labelCol = "steelblue")


plotSpectraPTM(sp[i])

## Comparing spectra with compareSpectra()

## Create a new Spectra object with the MS2 scans with sequences
## "SQILQQAGTSVLSQANQVPQTVLSLLR" and "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR".

i <- which(sp$sequence %in% c("SQILQQAGTSVLSQANQVPQTVLSLLR", "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR"))
spx <- sp[i]

spx$sequence

## Compare these with compareSpectra() and visualise the results using a heatmap
## using for example pheatmap::pheatmap().

mat <- compareSpectra(spx)
colnames(mat) <- rownames(mat) <- strtrim(spx$sequence, 2)

pheatmap::pheatmap(mat)

## Visualise/compare these scans: plotSpectra, plotSpectraOverlay,
## plotSpectraMirror, plotSpectraPTM

## Bonus: becore plotting with plotSprectra(), keep only peaks that have an
## intensity > 1e3 (see filterIntensity())

filterIntensity(spx, 1e3) |>
    plotSpectra()

plotSpectraMirror(spx[1], spx[2])

spx$sequence

plotSpectraMirror(spx[3], spx[4])

plotSpectraOverlay(spx[1:2], col = c("red", "blue"))

plotSpectraOverlay(spx[c(5, 3, 4)], col = c("red", "blue", "green"))

plotSpectraPTM(spx[1:2])

plotSpectraPTM(spx[3:5])


## Ex 6.1

## Get the data

px <- PXDataset("PXD022816")

mzf <- pxget(px, grep("mzML", pxfiles(px), value = TRUE)[1:3])
idf <- pxget(px, grep("mzID", pxfiles(px), value = TRUE)[1:3])

## Generate a Spectra object and a table of filtered PSMs. Visualise the total
## ion chromatograms and check the quality of the identification data by
## comparing the density of the decoy and target PSMs id scores for each file.

sp <- Spectra(mzf)

sp

length(sp)

table(msLevel(sp))

table(basename(sp$dataOrigin))
table(basename(sp$dataOrigin), msLevel(sp))

filterMsLevel(sp, 1) |>
    spectraData() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent,
               colour = basename(dataOrigin))) +
    geom_line()

id <- PSM(idf)

id

table(id$idFile)

names(id)

data.frame(id) |>
    ggplot(aes(x = MetaMorpheus.score,
               colour = isDecoy)) +
    geom_density() +
    facet_wrap(~ idFile)

max(id$PSM.level.q.value)

table(id$idFile, id$isDecoy)

id_filtered <- filterPSMs(id)


grep("scan=18216", id$spectrumID, value = TRUE)