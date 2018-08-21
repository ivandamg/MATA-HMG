########################################################################################
#' Import transcript-level abundances and estimated counts for gene-level analysis packages
########################################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
#install.packages("gplots")
library("gplots")
library("RColorBrewer")
library(Rgraphviz)
library(topGO)
library(edgeR)
library(limma)
#libraries
library('nlme')
#options("scipen"=100,digits=3) 
library('lme4')
library('ggplot2')
library(reshape)
library("multcomp")
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/WGCNA_1.49.tgz",type = "source", repos = NULL, lib=.Library) 
library(WGCNA)
library(pheatmap)
########################################################################################
# ANALYSIS AT GENE LEVEL 
# input data: Kallisto pseudo-alignment. 

tximport <- function(files,
                     type=c("none","kallisto","salmon","sailfish","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     reader=read.delim,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer,
                     collatedFiles,
                     ignoreTxVersion=FALSE) {
  
  type <- match.arg(type, c("none","kallisto","salmon","sailfish","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment='#') 
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    importer <- reader
  }
  
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      
      #####################################################################
      # some temporary code for detecting older fishes
      if ((i == 1) &
          (type %in% c("salmon","sailfish")) &
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length" 
        # because the comment lines have the same comment character
        # as the header, need to name the column names
        importer <- function(x) {
          tmp <- reader(x, comment="#")
          names(tmp) <- c("Name","Length","TPM","NumReads")
          tmp
        }
        # re-read the first file
        raw <- as.data.frame(importer(files[i]))
      }
      #####################################################################
      
      # does the table contain gene association or was an external tx2gene table provided?
      if (is.null(tx2gene) & !txOut) {
        # e.g. Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    message("")
    
    txi <- list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx,
                countsFromAbundance="no")
    
    # if the user requested just the transcript-level data:
    if (txOut) {
      return(txi)
    }
    
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
    # e.g. RSEM already has gene-level summaries
    # just combine the gene-level summaries across files
  } else {
    # stating the obvious:
    if (txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
    
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}

# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize tx-level matrices to gene-level
#' @export
summarizeToGene <- function(txi,
                            tx2gene,
                            ignoreTxVersion=FALSE,
                            countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
) {
  
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- txi$abundance
  countsMatTx <- txi$counts
  lengthMatTx <- txi$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {
    colnames(tx2gene) <- c("tx","gene")
    if (ignoreTxVersion) {
      txId <- sapply(strsplit(as.character(txId), "\\."), "[[", 1)
    }
    tx2gene$gene <- factor(tx2gene$gene)
    tx2gene$tx <- factor(tx2gene$tx)
    # remove transcripts (and genes) not in the abundances
    tx2gene <- tx2gene[tx2gene$tx %in% txId,]
    tx2gene$gene <- droplevels(tx2gene$gene)
    ntxmissing <- sum(!txId %in% tx2gene$tx)
    if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)
    sub.idx <- txId %in% tx2gene$tx
    abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
    countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
    lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
    txId <- txId[sub.idx]
    geneId <- tx2gene$gene[match(txId, tx2gene$tx)]
  }
  
  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
  message("summarizing counts")
  countsMat <- fastby(countsMatTx, geneId, colSums)
  message("summarizing length")
  
  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
  lengthMat <- weightedLength / abundanceMat   
  
  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)
  
  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))
  
  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)
  
  if (countsFromAbundance != "no") {
    countsSum <- colSums(countsMat)
    if (countsFromAbundance == "lengthScaledTPM") {
      newCounts <- abundanceMat * rowMeans(lengthMat)
    } else {
      newCounts <- abundanceMat
    }
    newSum <- colSums(newCounts)
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  if (ncol(m) > 1) {
    t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
  } else {
    matrix(sapply(idx, function(i) fun(m[i,,drop=FALSE])),
           dimnames=list(levels(f), colnames(m)))
  }
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
        # otherwise use the geometric mean of the lengths from the other samples
        idx <- is.nan(lengthMat[i,])
        lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
      }
    }
  }
  lengthMat
}

#################################################################################
#################################################################################
# transfrom transcripts to genes.
####################################################################

#Rirregularis
setwd('~/Documents/Co-inoculation_manuscript/data/Results/kallisto_Umap_Rirregularis')   
filesToProcess <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Documents/Co-inoculation_manuscript/data/Results/kallisto_Umap_Rirregularis/UmapAMF_V1_10_B1_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')

txi_AMF <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_AMF<-txi_AMF[[2]]#[, c(-15,-45,-46,-40,-41)] # exclude repetead lib
colnames(FOUR_VARS_AMF)<-gsub("UmapAMF_"," ",colnames(FOUR_VARS_AMF),perl=T)
colnames(FOUR_VARS_AMF)

#annotation infos
########################################################################################################################################
#R. irregularis
rirregularis_go<-read.delim('~/Documents/Co-inoculation_manuscript/data/Bioinformatics/Annotations/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
head(mesculenta_go)

# AMF
Mercator_Rirregularis<-read.table("~/Documents/Co-inoculation_manuscript/data/Bioinformatics/Annotations/Mercator_Rirregularis_database_v4",h=F,quote="\'")
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene","details","Type")
Mercator_Rirregularis$Type<-unlist( lapply(lapply(strsplit(as.character(Mercator_Rirregularis[,2],20),split="[.]"),function (x) x[1:2]),function (x) paste(x[1],x[2],sep ="_" )  ))


########################################################################################################################################
####################################################################
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CTRL"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_3"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V1_7"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_11"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_15"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_6"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_9"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_6"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_19"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V6_14"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_19"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_7"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_17"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_13"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_8"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_16b"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_8"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V6_11"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_4"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_8"), colnames(FOUR_VARS_AMF),invert =T)]


colnames(FOUR_VARS_AMF[,grep("V3", colnames(FOUR_VARS_AMF),invert=F)])

# exclude varieties temporal 
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V1"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8"), colnames(FOUR_VARS_AMF),invert =T)]


####################################
# DESIGN DEFINITION


VAR<-gsub( "m","",gsub( "m2","",gsub("_\\w+$"," ",colnames(FOUR_VARS_AMF),perl=T)    ) )
TREAT<-gsub("^ (\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CANB1","CAN","B1"))
#TREAT_AMF<-factor(TREAT,levels=c("CAN","B1","CANB1"))

design_AMF <- model.matrix(~TREAT_AMF*VAR)
#design_AMF <- model.matrix(~TREAT_AMF)


####################################################################
#################################################################################
# DATA FILTERING AND NORMALIZATION
DGE_AMF <- DGEList(FOUR_VARS_AMF)

#FILTERING  


keep_A <- rowSums(DGE_AMF$counts>50) >= 3
DGE_AMF<- DGE_AMF[keep_A,]

#NORMALIZATION
DGE_AMF <- calcNormFactors(DGE_AMF)

DGE_AMF_N <- voom(DGE_AMF, design_AMF,plot=TRUE)


### PLOT each gene
tratamiento<-unlist(lapply(strsplit(colnames(DGE_AMF$counts),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(DGE_AMF$counts),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(7:9,1:3,4:6)]

boxplot(DGE_AMF_N$E[grep("g7952.t1",rownames(DGE_AMF_N$E)),]~interacts,las=2)

#####################################################################
# Plot all the HMG-box transcription results

# only HMg-box genes wich transcription value is higher than sum of 3 samples is higher than 50.
# increase confidence of data

# detected 75 of 76 on RNAseq data, but only 25 genes displayed a confident read number per samples.
# The ones not detected could mean that were not active, or very low active

HMG_box<-read.table("~/Documents/Co-inoculation_manuscript/v7_NewPhytol/MATA_HMG_GENOMES/HMG-BOX_n6/HMG-box_Predictedn6_names_vf.txt")
HMG_box<-as.vector(HMG_box[,1])
# Make barplot
genes_amf_compe2<-DGE_AMF$counts[rownames(DGE_AMF$counts) %in% HMG_box, ] 
#DGE_AMF$counts
rownames(genes_amf_compe2)

genes_amf_compe2<-genes_amf_compe2[c(24,1:23,25),]


tratamiento<-unlist(lapply(strsplit(colnames(genes_amf_compe2),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_amf_compe2),"_"),function (x) x[1]))   ))


pdf("~/Documents/Co-inoculation_manuscript/v7_NewPhytol/MATA_HMG_GENOMES/HMG-BOX_n6/SupFig_allHMG-box_transcription.pdf",width=18, height=18,useDingbats = F)
par(mfrow=c(5,5),mar=c(15,9,3,1),mgp=c(5, 1, 0))

for (i in 1:length(genes_amf_compe2[,1])) {
stDevs <-tapply(genes_amf_compe2[i,],interaction(tratamiento,cultivar),sd)
means<-tapply(genes_amf_compe2[i,],interaction(tratamiento,cultivar),mean)
mp<-barplot(tapply(genes_amf_compe2[i,],interaction(tratamiento,cultivar),mean),
            beside=T,las =2,ylim=c(0,max(genes_amf_compe2[i,])*1.2), col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 2,cex.names = 2,cex.lab=2,
            ylab=paste("Gene-transcription ", rownames(genes_amf_compe2)[i],sep = "\n") )
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

