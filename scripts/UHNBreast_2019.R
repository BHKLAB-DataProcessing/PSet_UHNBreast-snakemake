options(stringsAsFactors = FALSE)
library(PharmacoGx)
library(readr)
library(rhdf5)
# library(gdata)
# library(readxl)
# library(openxlsx)
library(data.table)
library(reshape2)
library(Biobase)
library(CoreGx)
library(SummarizedExperiment)
library(parallel)
# library(biocompute)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- paste0(args[[1]], "download")
output_dir <- args[[1]]
filename <- args[[2]]
tools <- args[[3]]
transcriptome <- args[[4]]

tools <- gsub("-", "_", tools)
untar(file.path(input_dir, paste0(tools, ".tar.gz")), exdir = file.path(input_dir))
tool_path <- expand.grid(a = tools, b = transcriptome)
tool_path <- paste0(tool_path$a, "_", tool_path$b)

print(tool_path)

print(args)

standardize <- args[grep("filtered", args)]

standardizeRawDataConcRange <- function(sens.info, sens.raw) {
  unq.drugs <- unique(sens.info$drugid)

  conc.m <- data.table(melt(sens.raw[, , 1], as.is = TRUE))
  conc.m[, drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[, .(l = min(value, na.rm = T), r = max(value, na.rm = T)), c("drugid", "Var1")]
  conc.ranges[, Var1 := NULL]
  conc.ranges <- conc.ranges[, unique(.SD), drugid]
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[
    ,
    {
      sq <- sort(unique(c(l, r)))
      l <- sq[seq(1, length(sq) - 1)]
      r <- sq[seq(2, length(sq))]
      .(l = l, r = r)
    },
    drugid
  ]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N > 0)
    unlist(sapply(seq(1, N), function(ii) {
      return(sapply(seq(ii, N), function(jj) {
        return(seq(ii, jj))
      }))
    }), recursive = FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl) {
    if (length(indicies) == 1) {
      return(TRUE)
    }
    sq <- seq(indicies[1], indicies[length(indicies)] - 1)
    all(lr.tbl[["l"]][sq + 1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[, .N, drugid][, N], returnConsInts)

  names(per.drug.range.indicies) <- conc.ranges.disj[, unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]


  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug) {
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug) {
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) {
      return(c(lr.tbl[x[1], l], lr.tbl[x[length(x)], r]))
    }))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify = FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)

  conc.ranges <- conc.m[, .(l = min(value, na.rm = T), r = max(value, na.rm = T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)
  # tic()
  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug) {
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid == drug, .(l, r)], 1, function(rng) {
      conc.m[drugid == drug][conc.ranges[drugid == drug][l <= rng["l"]][r >= rng["r"], Var1], on = "Var1"][value >= rng["l"]][value <= rng["r"], .N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid == drug][which(num.points.in.range == max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })
  # toc()
  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug) {
    rng <- unlist(chosen.drug.ranges[[drug]][, .(l, r)])
    exp.out.range <- conc.ranges[drugid == drug][l > rng["l"] | r < rng["r"], Var1]
    return(exp.out.range)
  }))

  sens.raw[removed.experiments, , ] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]

  for (drug in unq.drugs) {
    rng <- unlist(chosen.drug.ranges[[drug]][, .(l, r)])
    myx <- conc.ranges.kept[drugid == drug, Var1]
    doses <- sens.raw[myx, , "Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, , "Dose"][which(which.remove, arr.ind = TRUE)] <- NA_real_
    sens.raw[myx, , "Viability"][which(which.remove, arr.ind = TRUE)] <- NA_real_

    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid == drug, "chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid == drug, "chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments, "rm.by.conc.range"] <- TRUE

  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}


# filter noisy curves from PSet (modified function to take into account standardized conc range)
filterNoisyCurves2 <- function(pSet, epsilon = 25, positive.cutoff.percent = .80, mean.viablity = 200, nthread = 1) {
  acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
    # for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp, , ], 2, as.numeric), stringsAsFactors = FALSE)
    if (!all(is.na(drug.responses))) {
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[, "delta"] <- .computeDelta(drug.responses$Viability)

      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)

      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)

      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
        (delta.sum < epsilon) &
        (max.cum.sum < (2 * epsilon)) &
        (mean(drug.responses$Viability) < mean.viablity)) {
        return(xp)
      }
    }
  }, mc.cores = nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy" = noisy, "ok" = acceptable))
}

.computeDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if (trunc) {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx) - 1]), 0))
  } else {
    return(c(xx[2:length(xx)] - xx[1:length(xx) - 1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if (trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2, simplify = TRUE))
  tt <- tt[which(((tt[, 2] - tt[, 1]) >= 2) == TRUE), ]
  if (is.null(nrow(tt))) {
    tt <- matrix(tt, ncol = 2)
  }
  cum.sum <- unlist(lapply(1:nrow(tt), function(x) {
    xx[tt[x, 2]] - xx[tt[x, 1]]
  }))
  return(max(cum.sum))
}

matchToIDTable <- function(ids, tbl, column, returnColumn = "unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)", Hmisc::escapeRegex(x), "((///)|$)"), tbl[, column])
    if (length(myx) > 1) {
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if (length(myx) == 0) {
      return(NA_character_)
    }
    return(tbl[myx, returnColumn])
  })
}

### curation cell/drug/tissue###

# cell
cell_all <- read.csv(file = file.path(input_dir, "cell_annotation_all.csv"), na.strings = c("", " ", "NA"))
curationCell <- cell_all[which(!(is.na(cell_all[, "Ben_Neel.cellid"]) & is.na(cell_all[, "Cescon.cellid"]))), ]
curationCell <- curationCell[, c("unique.cellid", "Ben_Neel.cellid", "Cescon.cellid")]
rownames(curationCell) <- curationCell[, "unique.cellid"]

# drug

drug_all <- read.csv(file = file.path(input_dir, "drugs_with_ids.csv"), na.strings = c("", " ", "NA"))
drug_all <- drug_all[, c("unique.drugid", "UHNBreast.drugid", "smiles", "inchikey", "cid", "FDA")]
# rownames(drug_all) <- drug_all[, "unique.drugid"]
curationDrug <- drug_all[which(!is.na(drug_all[, "UHNBreast.drugid"])), ]
curationDrug <- curationDrug[, c("unique.drugid", "UHNBreast.drugid")]
rownames(curationDrug) <- curationDrug[, "unique.drugid"]


# tissue
curationTissue <- cell_all[which(!(is.na(cell_all[, "Ben_Neel.cellid"]) & is.na(cell_all[, "Cescon.cellid"]))), ]
curationTissue <- curationTissue[, c("unique.tissueid", "Ben_Neel.tissueid", "Cescon.cellid")]
rownames(curationTissue) <- curationCell[, "unique.cellid"]


### cell line info###
cellline_info <- read.csv(file = file.path(input_dir, "bc_cellines_neel_subtypes.csv"), na.strings = c("", " ", "NA"))
cellline_info$cellid <- matchToIDTable(ids = cellline_info$cellid, tbl = curationCell, column = "Ben_Neel.cellid", returnColumn = "unique.cellid")
dupCell <- cellline_info$cellid[duplicated(cellline_info$cellid)]
for (cell in dupCell) {
  myx <- which(cellline_info$cellid == cell)
  cellline_info[myx[1], ] <- apply(cellline_info[myx, ], 2, function(x) {
    return(unique(na.omit(x)))
  })
  cellline_info <- cellline_info[-myx[2], ]
}
rownames(cellline_info) <- cellline_info$cellid
cellline_info$tissueid <- curationTissue$unique.tissueid[match(cellline_info$cellid, rownames(curationTissue))]
# cellline_info["HCT 116","tissueid"] <- "large_intestine"
# cellline_info["DU145","tissueid"] <- "prostate"
# cellline_info[83,"cellid"] <- "HCT 116"
# cellline_info[84,"cellid"] <- "DU145"

### drug info###
drug_info <- data.frame("drugid" = curationDrug$unique.drugid)
rownames(drug_info) <- drug_info$drugid
drug_all <- drug_all[rownames(drug_info), ]
drug_info[, c("smiles", "inchikey", "cid", "FDA")] <- drug_all[, c("smiles", "inchikey", "cid", "FDA")]


# import recomputed sensitivity
load(file.path(input_dir, "UHN_recomputed.RData"))

sensitivityData_final <- lapply(sensitivityData, function(x) {
  ibx <- rownames(x$info)[x$info$drug == "nothing"]

  if (length(ibx) != 0) {
    x$info <- x$info[-match(ibx, rownames(x$info)), , drop = F]
    x$profile <- x$profile[-match(ibx, rownames(x$profile)), , drop = F]
    x$raw <- x$raw[-match(ibx, rownames(x$raw)), , , drop = F]
    x$raw.complete <- x$raw.complete[-match(ibx, rownames(x$raw.complete)), , , drop = F]
  }
  return(x)
})

sensitivityData <- sensitivityData_final


rm(sensitivityData_final)

metaData <- lapply(sensitivityData, function(x) {
  xx <- x$info
  xx[, "exp"] <- rownames(xx)
  rownames(xx) <- NULL


  return(xx)
})

metaData <- do.call(rbind, metaData)


sum(duplicated(metaData$exp)) ## 84 duplicates

table(metaData$cell)
table(metaData$drug)

cell_unannot <- unique(metaData$cell)

cell_annot <- gsub(x = cell_unannot, pattern = "^[0-9]+-", replacement = "")
cell_annot <- gsub(cell_annot, pattern = "_2", replacement = "")
cell_annot <- toupper(cell_annot)

cell_unannot <- cell_unannot[order(cell_annot)]
cell_annot <- cell_annot[order(cell_annot)]

cellMatch <- data.frame("annotated" = cell_annot, "unannotated" = cell_unannot)


badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

closeMatches <- lapply(cellMatch[, "annotated"], function(x) {
  myx <- cell_all$unique.cellid[which(x == toupper(cell_all[, "unique.cellid"]))]
  if (length(myx) == 0) {
    if (grepl(pattern = "^MDA", x) && x != "MDA361-20%" && x != "MDA143") { ## different way of writing these cell lines
      print(x)
      x <- paste(strsplit(x, split = "A")[[1]][1], "A", "MB", strsplit(x, split = "A")[[1]][2])
      #      print(x)
    }
    if (grepl(pattern = "^UAC[1-9]+", x)) { # common typo
      print(x)
      x <- paste(strsplit(x, split = "C")[[1]][1], "CC", strsplit(x, split = "C")[[1]][2])
    }
    myx <- grep(pattern = toupper(gsub(badchars, "", x)), x = toupper(gsub(badchars, "", cell_all$unique.cellid)))
    myx <- cell_all$unique.cellid[myx]
  }


  if (x == "MFN223") {
    myx <- "MFM-223"
  } ## typo in data
  if (x == "AU655") {
    myx <- "AU565"
    print(x)
  } ## typo in data
  if (x == "HCL70") {
    myx <- "HCC70"
  } ## typo in data
  if (x == "MPE600") {
    myx <- "600MPE"
  } ## typo in data
  if (x == "MDA361-20%") {
    print(x)
    myx <- "MDA-MB-361"
  } ## typo in data
  if (x == "436") {
    print(x)
    myx <- "MDA-MB-436"
  } ## typo in data
  if (x == "MDA143") {
    myx <- "MDA-MB-134-VI"
    print(x)
    #    print(myx)
  } ## typo in data
  if (x == "OCUB-1" | x == "OCUB1") {
    myx <- "OCUB-M"
  }
  if (x == "ACC202") {
    print(x)
    myx <- "HCC202"
  } ## Closest match (is child of cell line)
  return(ifelse(length(myx) > 0, myx, NA))
})
sum(!sapply(closeMatches, is.na))


cellMatch[, "closeMatches"] <- unlist(closeMatches)


warning("All matches are now OK, but this may change if other cells are added in the future")
cellMatch$annotated[!is.na(cellMatch$closeMatches)] <- na.omit(cellMatch$closeMatches)
metaData$newcell <- cellMatch$annotated[match(metaData$cell, cellMatch$unannotated)]
metaData$newexp <- paste(metaData$newcell, metaData$drug, sep = "_")
for (exp in unique(metaData$newexp)) {
  myx <- which(metaData$newexp == exp)
  metaData$newexp[myx] <- paste(metaData$newexp[myx], rep("rep", length(myx)), 1:length(myx), sep = "_")
}
profilesData <- lapply(sensitivityData, function(x) {
  xx <- x$profile
  rownames(xx) <- NULL
  return(xx)
})
profilesData <- do.call(rbind, profilesData)
rownames(profilesData) <- metaData$newexp
require(abind)
rawData <- lapply(sensitivityData, function(x) {
  xx <- x["raw"]
  rownames(xx) <- NULL
  return(xx)
})

maxConc <- max(unlist(lapply(rawData, function(x) {
  return(dim(x[[1]])[2])
})))

for (i in 1:length(rawData)) {
  currentDim <- dim(rawData[[i]][[1]])[2]
  if (currentDim < maxConc) {
    remainingDim <- maxConc - currentDim
    tmpRaw <- array(NA, dim = c(dim(rawData[[i]][[1]])[1], remainingDim, dim(rawData[[i]][[1]])[3]), dimnames = list(NULL, NULL, unlist(dimnames(rawData[[i]][[1]])[[3]])))
    finalRaw <- array(NA, dim = c(dim(rawData[[i]][[1]])[1], maxConc, dim(rawData[[i]][[1]])[3]), dimnames = list(NULL, NULL, unlist(dimnames(rawData[[i]][[1]])[[3]])))
    for (j in 1:dim(rawData[[i]][[1]])[1]) {
      finalRaw[j, , ] <- rbind(rawData[[i]][[1]][j, , ], tmpRaw[j, , ])
    }

    rawData[[i]] <- finalRaw
  } else {
    rawData[[i]] <- rawData[[i]][[1]]
  }
}

lapply(rawData, function(x) {
  print(dim(x))
})

rawData <- abind(rawData, along = 1)
dimnames(rawData)[[1]] <- metaData$newexp
rownames(metaData) <- metaData$newexp
####### bind everything dtogether
metaData[, c("cell", "exp", "newexp")] <- NULL
colnames(metaData) <- gsub(pattern = "new", x = colnames(metaData), replace = "")
colnames(metaData) <- c("drugid", "file", "cellid")
metaData <- metaData[, c("drugid", "cellid", "file")]
sensitivity <- list(info = metaData, profiles = profilesData, raw = rawData)


## sensitivity.info##

x2 <- as.character(matchToIDTable(ids = sensitivity$info[, "drugid"], tbl = curationDrug, column = "UHNBreast.drugid", returnColumn = "unique.drugid"))
sensitivity$info[, "drugid"] <- x2

x3 <- as.character(matchToIDTable(ids = tolower(gsub(badchars, "", sensitivity$info[, "cellid"])), tbl = curationCell, column = "Ben_Neel.cellid", returnColumn = "unique.cellid"))
sensitivity$info[, "cellid"] <- x3


sensitivity$profiles$auc_recomputed <- sensitivity$profiles$auc_recomputed / 100
colnames(sensitivity$profiles)[2] <- "aac_recomputed"
colnames(sensitivity$profiles)[3] <- "slope_recomputed"


# RNA-seq

summarizeRnaSeq <- function(dir,
                            features_annotation,
                            samples_annotation,
                            method) {
  library(Biobase)
  library(readr)
  library(tximport)

  load(features_annotation)

  tx2gene <- as.data.frame(cbind("transcript" = tx2gene$transcripts, "gene" = tx2gene$genes))

  files <- list.files(dir, recursive = TRUE, full.names = T)
  if (method == "kallisto") {
    resFiles <- grep("abundance.h5", files)
  } else {
    resFiles <- grep("quant.sf", files)
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))

  if (features_annotation == file.path(input_dir, "Ensembl.v99.annotation.RData")) {
    txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else {
    txi <- tximport(resFiles, type = method, tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)
  }

  head(txi$counts[, 1:5])
  dim(txi$counts)

  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp), ]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp), ]
  annotation(gene.exp) <- "rnaseq"

  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count), ]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count), ]
  annotation(gene.count) <- "rnaseq"

  txii <- tximport(resFiles, type = method, txOut = T)

  if (features_annotation == file.path("Ensembl.v99.annotation.RData")) {
    # remove non-coding transcripts in ensembl
    rownames(txii$abundance) <- gsub("\\..*", "", rownames(txii$abundance))
    txii$abundance[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
    missing_transcript <- rownames(txii$abundance)[which(!rownames(txii$abundance) %in% features_transcript$transcript_id)]
    txii$abundance <- txii$abundance[-which(rownames(txii$abundance) %in% missing_transcript), ]
  }

  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 0.001))
  if (features_annotation == file.path(input_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(input_dir, "Gencode.v33lift37.annotation.RData")) {
    featureNames(transcript.exp) <- gsub("\\|.*", "", featureNames(transcript.exp))
    fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
  } else {
    fData(transcript.exp) <- features_transcript[featureNames(transcript.exp), ]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp), ]
  annotation(transcript.exp) <- "isoform"


  if (features_annotation == file.path("Ensembl.v99.annotation.RData")) {
    # remove non-coding transcripts in ensembl
    rownames(txii$counts) <- gsub("\\..*", "", rownames(txii$counts))
    txii$counts <- txii$counts[-which(rownames(txii$counts) %in% missing_transcript), ]
  }
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[, 1:length(resFiles)] + 1))
  if (features_annotation == file.path(input_dir, "Gencode.v33.annotation.RData") || features_annotation == file.path(input_dir, "Gencode.v33lift37.annotation.RData")) {
    featureNames(transcript.count) <- gsub("\\|.*", "", featureNames(transcript.count))
    fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
  } else {
    fData(transcript.count) <- features_transcript[featureNames(transcript.count), ]
  }
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count), ]
  annotation(transcript.count) <- "isoform"

  return(list(
    "rnaseq" = gene.exp,
    "rnaseq.counts" = gene.count,
    "isoforms" = transcript.exp,
    "isoforms.counts" = transcript.count
  ))
}

rnaseq.sampleinfo <- read.csv(file.path(input_dir, "uhn_metadata_new.csv"), stringsAsFactors = FALSE, row.names = 2)
rnaseq.sampleinfo[, "cellid"] <- as.character(matchToIDTable(ids = tolower(gsub(badchars, "", rnaseq.sampleinfo[, "cell.id"])), tbl = curationCell, column = "Ben_Neel.cellid", returnColumn = "unique.cellid"))
rnaseq.sampleinfo$cell.line <- rnaseq.sampleinfo$cell.id
rnaseq.sampleinfo$cell.id <- NULL

rnaseq_results <- c()
for (r in 1:length(tool_path)) {
  print(tool_path[r])
  if (length(grep(pattern = "Kallisto", x = tool_path[r])) > 0) {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir <- paste0("uhn_rnaseq_", gsub(".", "_", tolower(tool), fixed = T), "/", tool, "/")
    rnatool <- "kallisto"
  } else {
    tool <- sub("(_[^_]+)_.*", "\\1", tool_path[r])
    # tdir <- paste0("uhn_rnaseq_", gsub(".", "_", tolower(tool), fixed = T), "/", tool, "/")
    rnatool <- "salmon"
  }


  if (length(grep(pattern = "lift37", x = tool_path[r])) > 0) {
    annot <- file.path(input_dir, "Gencode.v33lift37.annotation.RData")
  } else if (length(grep(pattern = "v33", x = tool_path[r])) > 0) {
    annot <- file.path(input_dir, "Gencode.v33.annotation.RData")
  } else {
    annot <- file.path(input_dir, "Ensembl.v99.annotation.RData")
  }
  print(annot)

  rnaseq <- summarizeRnaSeq(
    dir = file.path(input_dir, tool, tool_path[r]),
    features_annotation = annot,
    samples_annotation = rnaseq.sampleinfo,
    method = rnatool
  )
  rnaseq_results <- c(rnaseq_results, c(
    rnaseq <- setNames(rnaseq, paste0(tool, ".", names(rnaseq)))
  ))
}


z <- list()

z <- c(z, c(
  rnaseq_results
))


# add missing cells to cell_info
rnaseq_cellid_all <- pData(rnaseq_results[[1]])[, "cellid"]
cellnall <- CoreGx::.unionList(rownames(cellline_info), rnaseq_cellid_all, sensitivity$info$cellid)
newcells <- setdiff(cellnall, rownames(cellline_info))
newRows <- matrix(NA_character_, nrow = length(newcells), ncol = ncol(cellline_info))

rownames(newRows) <- newcells
colnames(newRows) <- colnames(cellline_info)
newRows[, "cellid"] <- newcells

cellline_info <- rbind(cellline_info, newRows)

cellsPresent <- sort(CoreGx::.unionList(sensitivity$info$cellid, rnaseq_cellid_all))
cellline_info <- cellline_info[cellsPresent, ]

cellline_info$tissueid <- curationTissue[rownames(cellline_info), "unique.tissueid"]
cellline_info$cellid <- rownames(cellline_info)


# add cellosaurus disease type to cell-info
disease <- cell_all$Cellosaurus.Disease.Type[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$Cellosaurus.Disease.Type <- disease

# add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$Cellosaurus.Accession.id <- assession

# add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$PharmacoDB.id <- pdb

# add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$unique.tissueid.fromstudies <- study_tissue

# add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$CellLine.Type <- cell_type

# add metastatic info to cell_info
metastatic <- cell_all$Metastatic[match(cellline_info$cellid, cell_all$unique.cellid)]
cellline_info$Metastatic <- metastatic

curationCell <- curationCell[rownames(cellline_info), ]

drug_info <- drug_info[unique(sensitivity$info$drugid), ]

curationDrug <- curationDrug[rownames(drug_info), ]
curationTissue <- curationTissue[rownames(cellline_info), ]


.converteSetToSE <- function(eSets) {
  SEfinal <- lapply(
    eSets,
    function(eSet) {
      # Change rownames from probes to EnsemblGeneId for rna data type
      if (grepl("^rna$", Biobase::annotation(eSet))) {
        rownames(eSet) <- Biobase::fData(eSet)$EnsemblGeneId
      }

      # Build summarized experiment from eSet
      SE <- SummarizedExperiment::SummarizedExperiment(
        ## TODO:: Do we want to pass an environment for better memory efficiency?
        assays = S4Vectors::SimpleList(as.list(Biobase::assayData(eSet))),
        # Switch rearrange columns so that IDs are first, probes second
        rowData = S4Vectors::DataFrame(Biobase::fData(eSet),
          rownames = rownames(Biobase::fData(eSet))
        ),
        colData = S4Vectors::DataFrame(Biobase::pData(eSet),
          rownames = rownames(Biobase::pData(eSet))
        ),
        metadata = list(
          "experimentData" = eSet@experimentData,
          "annotation" = Biobase::annotation(eSet),
          "protocolData" = Biobase::protocolData(eSet)
        )
      )
      ## TODO:: Determine if this can be done in the SE constructor?
      # Extract names from expression set
      SummarizedExperiment::assayNames(SE) <- Biobase::assayDataElementNames(eSet)
      mDataType <- Biobase::annotation(eSet)
      eSets[[mDataType]] <- SE
    }
  )
  # setNames(pSet@molecularProfiles, names(eSets))
  return(SEfinal)
}

z <- .converteSetToSE(z)

if (length(standardize) > 0) {

  # standardize <- standardizeRawDataConcRange(sens.info = sensitivity$info, sens.raw = sensitivity$raw)
  # sensitivity$info <- standardize$sens.info
  # sensitivity$raw <- standardize$sens.raw
} else {
  print("unfiltered PSet")
}


UHNBreast2019 <- PharmacoSet(
  name = "UHNBreast",
  molecularProfiles = z,
  cell = cellline_info,
  drug = drug_info,
  sensitivityInfo = sensitivity$info,
  sensitivityRaw = sensitivity$raw,
  sensitivityProfiles = sensitivity$profiles,
  curationCell = curationCell,
  curationDrug = curationDrug,
  curationTissue = curationTissue,
  datasetType = "both"
)


if (length(standardize) > 0) {
  noisy_out <- filterNoisyCurves2(UHNBreast2019)
  print("filter done")
  UHNBreast2019@sensitivity$profiles[noisy_out$noisy, ] <- NA
} else {
  print("unfiltered PSet")
}

UHNBreast2019@annotation$version <- 2
saveRDS(UHNBreast2019, file = file.path(output_dir, filename))

# dataset <- "UHNBreast"
# #output ORCESTRA_ID and Pachyderm commit id
# write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
# write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
# pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
# write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)


# ###CREATE BIOCOMPUTE OBJECT###
# tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
# tools <- rnaseq_select[tools]
# tools <- gsub("-", "_", tools)

# tool <- gsub("\\_.*","", tools)
# version <- gsub(".*_","", tools)

# print(tool)
# print(version)

# if (length(tools) > 1){

#   if (isTRUE(any(duplicated(tool)))){

#     tool_name <- tool[1]
#     version <- c(version[1],version[2])

#     if (tool[1] == "Kallisto" || tool[2] == "Kallisto"){

#       parameter <- "-t"
#       parameter_value <- "2"
#       uri <- rep("https://pachterlab.github.io/kallisto/",2)

#     } else{
#       parameter <- c("-l", "--validateMappings")
#       parameter_value <- c("A", "")

#       uri <- rep("https://combine-lab.github.io/salmon/",2)
#     }

#   } else {

#     tool <- c(tool[1], tool[2])
#     v <- version
#     #version <- c(version[1],version[2])
#     if (v[1] == "0.8.2" || v[2] == "0.8.2"){
#       parameter <- c("-t","-l")
#       parameter_value <- c("2", "A")
#     } else{
#       parameter <- c("-t","-l", "--validateMappings")
#       parameter_value <- c("2", "A","")
#     }
#     uri <- c("https://pachterlab.github.io/kallisto/","https://combine-lab.github.io/salmon/")
#   }

# } else {

#   if (tool == "Kallisto"){

#     parameter <- "-t"
#     parameter_value <- c("2")
#     uri <- "https://pachterlab.github.io/kallisto/"
#   } else{
#     if (version == "0.8.2"){
#       parameter <- paste0("-l")
#       parameter_value <- c("A")
#     } else{
#       parameter <- c("-l", "--validateMappings")
#       parameter_value <- c("A","")
#     }
#     uri <- "https://combine-lab.github.io/salmon/"
#   }
# }

# ###Selection Transcriptome###

# if (transcriptome == "Gencode_v33"){

#   transcriptome_link <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz"

# } else if (transcriptome == "Gencode_v33lift37"){

#   transcriptome_link <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/GRCh37_mapping/gencode.v23lift37.annotation.gtf.gz"
# } else {

#   transcriptome_link <- "ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz"
# }

# print(tool)
# print(version)
# print(parameter)
# print(parameter_value)
# print(uri)
# print(transcriptome)
# print(transcriptome_link)


# ###########################
# #####Provenance Domain#####
# ###########################

# #Created and modified dates
# #Sys.setenv(TZ = "EST")
# created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")

# #Contributions
# contributors <- data.frame(
#   "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
#   "affiliation" = c(rep("University Health Network", 3)),
#   "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
#   "contribution" = c("createdBy","createdBy","authoredBy"),
#   "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
#   stringsAsFactors = FALSE
# )

# #License
# license <- "https://opensource.org/licenses/Apache-2.0"

# #Name of biocompute object
# name <- "UHNBreast_2019"

# #Version of biocompute object
# bio_version <- "1.0.0"

# #Embargo (none)
# embargo <- c()

# #Derived from and obsolete after (none)
# derived_from <- c()
# obsolete_after <- c()

# #reviewers (none)
# review <- c()

# #compile domain
# provenance <- compose_provenance_v1.3.0(
#   name, bio_version, review, derived_from, obsolete_after,
#   embargo, created, modified, contributors, license
# )
# provenance %>% convert_json()


# ############################
# #####Description Domain#####
# ############################
# times_rnaseq <- as.POSIXct("2020-01-20T1:04:10", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
# #Keywords and platform info
# keywords <- c("Biomedical", "Pharmacogenomics", "Cellline", "Drug")
# platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")

# #Metadata for each pipeline step
# pipeline_meta <- data.frame(
#   "step_number" = c("1","2","3","4"),
#   "name" = c("Expression processing",
#              "Curated Sample and treatment identifier compilation",
#              "Drug sensitivity processing",
#              "Build data object"),
#   "description" = c("Pseudoalignment of reads",
#                     "Download of appropriate sample and treatment identifiers from GitHub (curations performed by BHK Lab - http://bhklab.ca)",
#                     "Process sensitivity data",
#                     "Building of ORCESTRA data object"),
#   "version" = c(1.0,1.0,1.0,1.0),
#   stringsAsFactors = FALSE
# )

# #Inputs for each pipeline step
# pipeline_input <- data.frame(
#   "step_number" = c("1","1","2","2","3","4"),
#   "filename" = c("Raw RNA-seq data",
#                  "Transcriptome",
#                  "Sample annotation data",
#                  "Treatment annotations",
#                  "Raw sensitivity data",
#                  "Script for data object generation"),
#   "uri" = c(
#     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73526",
#     transcriptome_link,
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cell_annotation_all.csv",
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/drugs_with_ids.csv",
#     " https://codeocean.com/capsule/6718332/",
#     "https://github.com/BHKLAB-Pachyderm/getUHNBreast_2019/UHNBreast_2019.R"
#   ),
#   "access_time" = c(times_rnaseq,times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )


# #Outputs for each pipeline step
# pipeline_output <- data.frame(
#   "step_number" = c("1","2","2","3","4"),
#   "filename" = c("Processed expression data",
#                  "Downloaded sample annotations",
#                  "Downloaded treatment annotations",
#                  "Process and compile sensitivity data in parallel",
#                  "Data object"),
#   "uri" = c(
#     "https://orcestradata.blob.core.windows.net/uhn/2019/RNA-seq/",
#     "/pfs/downAnnotations/cell_annotation_all.csv",
#     "/pfs/downAnnotations/drugs_with_ids.csv",
#     "/pfs/calculateuhn2019raw/UHNRecomputed_2019.RData",
#     "/pfs/out/UHNBreast.rds"
#   ),
#   "access_time" = c(times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )

# #xref (none)
# xref <- c()

# #pipeline prereq (none)
# pipeline_prerequisite <- c()

# #compile domain
# description <- compose_description_v1.3.0(
#   keywords, xref, platform,
#   pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
# )
# description %>% convert_json()


# ############################
# ######Execution Domain######
# ############################

# script <- c()
# script_driver <- c()

# #software/tools and its versions used for data object creation
# software_prerequisites <- data.frame(
#   "name" = c(tool, "Pachyderm", "Docker Image"),
#   "version" = c(version, "1.9.3", "v3"),
#   "uri" = c(
#     uri,
#     "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/pharmacogx2.0"
#   ),
#   stringsAsFactors = FALSE
# )

# software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
# software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))

# external_data_endpoints <- c()
# environment_variables <- c()

# execution <- compose_execution_v1.3.0(
#   script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
# )
# execution %>% convert_json()


# ############################
# ######Extension Domain######
# ############################

# #repo of scripts/data used
# scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
# scm_type <- "git"
# scm_commit <- c()
# scm_path <- c()
# scm_preview <- c()

# scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
# scm %>% convert_json()

# extension <- compose_extension_v1.3.0(scm)
# extension %>% convert_json()

# ############################
# ######Parametric Domain#####
# ############################

# df_parametric <- data.frame(
#   "param" = parameter,
#   "value" = parameter_value,
#   "step" = c(1),
#   stringsAsFactors = FALSE
# )

# parametric <- compose_parametric_v1.3.0(df_parametric)
# parametric %>% convert_json()



# ############################
# ######Usability Domain######
# ############################

# #usability of our data objects
# text <- c(
#   "Pipeline for creating UHNBreast 2019 data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data."
# )

# usability <- compose_usability_v1.3.0(text)
# usability %>% convert_json()


# ######################
# ######I/O Domain######
# ######################

# input_subdomain <- data.frame(
#   "filename" = c("Raw RNA-seq data",
#                  "Transcriptome",
#                  "Sample annotation data",
#                  "Treatment annotations",
#                  "Raw sensitivity data",
#                  "Script for data object generation"),
#   "uri" = c(
#     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73526",
#     transcriptome_link,
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/cell_annotation_all.csv",
#     "https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/drugs_with_ids.csv",
#     " https://codeocean.com/capsule/6718332/",
#     "https://github.com/BHKLAB-Pachyderm/getUHNBreast_2019/UHNBreast_2019.R"
#   ),
#   "access_time" = c(times_rnaseq,times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )

# output_subdomain <- data.frame(
#   "mediatype" = c("tar.gz", "csv", "csv", "RData", "RDS"),
#   "uri" = c(
#     "https://orcestradata.blob.core.windows.net/uhn/2019/RNA-seq/",
#     "/pfs/downAnnotations/cell_annotation_all.csv",
#     "/pfs/downAnnotations/drugs_with_ids.csv",
#     "/pfs/calculateuhn2019raw/UHNRecomputed_2019.RData",
#     "/pfs/out/UHNBreast.rds"
#   ),
#   "access_time" = c(times_rnaseq,created,created,created,created),
#   stringsAsFactors = FALSE
# )

# io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
# io %>% convert_json()


# ########################
# ######Error Domain######
# ########################

# empirical <- c()
# algorithmic <- c()

# error <- compose_error(empirical, algorithmic)
# error %>% convert_json()


# ####Retrieve Top Level Fields####
# tlf <- compose_tlf_v1.3.0(
#   provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# tlf %>% convert_json()


# ####Complete BCO####

# bco <- biocompute::compose_v1.3.0(
#   tlf, provenance, usability, extension, description,
#   execution, parametric, io, error
# )
# bco %>% convert_json() %>% export_json("/pfs/out/UHN_2019_BCO.json") %>% validate_checksum()
