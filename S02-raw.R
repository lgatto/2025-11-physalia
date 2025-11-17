library("Spectra")

spd <- DataFrame(msLevel = c(1L, 2L),
                 rtime = c(1.1, 1.2))

spd

spd$mz <- list(c(100, 103.2, 104.4, 111), c(45.6, 120.4, 190.1))
spd$intensity <- list(c(200, 400, 100, 13.2), c(111, 23, 370))

spd

sp <- Spectra(spd)

class(sp)

sp

spectraVariables(sp)

spectraData(sp)

peaksData(sp)[[1]]
peaksData(sp)[[2]]

sp[1]

sp[1:2]

mzf
sp <- Spectra(mzf)
sp

print(object.size(sp), units = "Mb")

## Ex 1
## - Repeat the data manipulations above.
## - Check the number of scans in the object with length().
## - Note the difference in the first line when showing the object in the
##   console. We will get back to this idea of backend later.
## - Extract the ms levels by using the msLevel column from the spectraData()
##   and my using the msLevel() accessor.
##
## See also https://rformassspectrometry.github.io/Spectra/
## and lecture notes.

sp

spectraVariables(sp)

spectraData(sp)

peaksData(sp)

peaksData(sp[1])[[1]]

length(sp)

spectraData(sp)$msLevel
msLevel(sp)
sp$msLevel

peaksData(sp[1:5])
intensity(sp[1:5])
mz(sp[1:5])

range(rtime(sp))
range(rtime(sp)/60)

## Ex 2
## - How many MS levels are there? How many scans of each level?

unique(msLevel(sp))
table(msLevel(sp))

## - Extract the index of the MS2 spectrum with the highest base peak intensity.

spectraVariables(sp)

table(msLevel(sp) == 2)

sp[msLevel(sp) == 2]

sp2 <- filterMsLevel(sp, 2)

sp_i <- sp2[which.max(sp2$basePeakIntensity)]

plotSpectra(sp_i)

## Backends

x <- filterMsLevel(sp, 2)[1:100]
x

peaksData(x)

print(object.size(x), units = "Kb")

xm <- setBackend(x, MsBackendMemory())

print(object.size(xm), units = "Kb")

xh <- setBackend(x, MsBackendHdf5Peaks(), hd5path = tempdir())
xh

dataStorage(xh)[1]

## Visualisation

## Ex chrom: The chromatogram can be created by extracting the totIonCurrent and
## rtime variables for all MS1 spectra. Annotate the spectrum of interest
## (2807th scan).

spectraVariables(sp)

plot(rtime(sp), tic(sp), type = "l")

plot(rtime(filterMsLevel(sp, 1)), tic(filterMsLevel(sp, 1)),
     type = "l")

sp1 <- filterMsLevel(sp, 1)
plot(rtime(sp1), tic(sp1), type = "l")
abline(v = rtime(sp)[2807], col = "red")

library(tidyverse)

spectraData(sp1) |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()


## The filterPrecursorScan() function can be used to retain a set parent (MS1)
## and children scans (MS2), as defined by an acquisition number. Use it to
## extract the MS1 scan of interest and all its MS2 children.

sp[2807]

spx <- filterPrecursorScan(sp, 2807)
spx

## Plot the MS1 spectrum of interest and highlight all the peaks that will be
## selected for MS2 analysis.

range(mz(spx[1])[[1]])
range(mz(spx[2])[[1]])

plotSpectra(spx[1])

plotSpectra(spx[1], xlim = c(400, 1000))

sort(precursorMz(spx))

abline(v = precursorMz(spx)[-1], col = "grey")
abline(v = precursorMz(spx)[2], col = "red")

range(precursorMz(sp), na.rm = TRUE)
range(precursorCharge(sp), na.rm = TRUE)
table(precursorCharge(sp))

## Zoom in mz values 521.1 and 522.5 to reveal the isotopic envelope of that
## peak.

plotSpectra(spx[1], xlim = c(521.1, 522.5))

abline(v = precursorMz(spx)[-1], col = "grey")
abline(v = precursorMz(spx)[2], col = "red")

## The plotSpectra() function is used to plot all 10 MS2 spectra in one call.

spx[-1]

plotSpectra(spx[-1])

plotSpectra(spx[7],
            xlim = c(126, 132))

mzLabel <- function(z) {
    z <- peaksData(z)[[1L]]
    lbls <- format(z[, "mz"], digits = 4)
    lbls[z[, "intensity"] < 1e5] <- ""
    lbls
}

plotSpectra(spx[7],
            xlim = c(126, 132),
            labels = mzLabel,
            labelOffset = 0.1,
            labelPos = 2)

## ?plotSpectra

## Filter MS2 level spectra and find any 2 MS2 spectra that have matching
## precursor peaks based on the precursor m/z values.

sp2 <- filterMsLevel(sp, 2L)

sp2

anyDuplicated(precursorMz(sp2))


i <- which(precursorMz(sp2) == precursorMz(sp2)[37])

precursorMz(sp2)[i]

sp2i <- sp2[i]


## Visualise the matching pair using the plotSpectraOverlay() and
## plotSpectraMirror() functions.

plotSpectra(sp2i)

plotSpectraOverlay(sp2i, col = c("red", "steelblue"),
                   xlim = c(100, 400))

plotSpectraMirror(sp2i[1], sp2i[2])