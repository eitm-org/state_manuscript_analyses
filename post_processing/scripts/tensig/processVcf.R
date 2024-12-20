#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(VariantAnnotation)
library(rhdf5)
library(rtracklayer)
# PLEASE adjust paths
source("post_processing/scripts/tensig/mutations.R")
load("post_processing/scripts/tensig/constants.RData") #download constants.Rdata from here: https://drive.google.com/file/d/1jZVpvFOP8lOLKY1pTt1m8AUK9kyD-3u3/view

if(!require(reticulate)){
    install.packages("reticulate")
    library(reticulate)
}
source_python('post_processing/scripts/constants.py')

genomePath <- file.path(refs_dir, "GRCh38/GRCh38.primary_assembly.genome_X.fa") #"genome.fa.gz" # path to the reference genome

liftOverAnnots <- function(gr){
  chain_path = file.path(refs_dir, '/gh38_granges/hg19ToHg38.over.chain')
  ch = rtracklayer::import.chain(chain_path)
  seqlevelsStyle(gr) = "UCSC"  # necessary
  gr = liftOver(gr, ch)
  unlist(gr)
}


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}
readVcfSave <- function(path) {
    # may need some user modifcation
    vcf <- VariantAnnotation::readVcf(path)
    vcf <- vcf[seqnames(vcf) %in% paste0('chr', c(1:22, "X", "Y"))] # filt(vcf) == "PASS" &
    elementLengths <- elementNROWS(alt(vcf))
    vcf <- vcf[elementLengths==1]
}

processVcf <- function(vcf, ts, epi, nuc, rt) {
    
    # extracts the trinucleotide context of each single base substitution
    tnc <- getTrinucleotideContext(vcf, genomePath)
    # annotates the mutation substitution
    sub <- getTrinucleotideSubs(vcf, tnc)
    # extracts transcription directionality for each substitution
    ts <- getStrandOrientation(vcf, ts)
    # extracts replication directionality for each substitution
    rs <- getStrandOrientation(vcf, rt)
    # extracts epigenetic state for each substitution
    ep <- getChromatinState(vcf, epi)
    # extracts the nucleosome directionality for each substitution
    nu <- getNucleosomeState(vcf, nuc)
    # extracts clustering state
    cl <- getClustering(vcf)
    t <- table(ts=ts, rs=rs, ep=ep, nu=nu, cl=cl, sub=sub)
    t[,,,,,SUB]
}

# output file is specified in the last argument
outputPath <- args[length(args)]
inputDir = args[1]

inputPaths = sort(list.files(inputDir, pattern=args[2], full.names = T, recursive = TRUE))

chunk_vcfs = args[3]
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

TS38 = liftOverAnnots(TS)
EPI38 = liftOverAnnots(EPI)
NUC38 = liftOverAnnots(NUC)
RT38 = liftOverAnnots(RT)

if(chunk_vcfs == FALSE){
  vcf <- lapply(inputPaths, readVcfSave)
  snvTensor <- sapply(vcf, function(x) processVcf(x[isSNV(x)], TS38, EPI38, NUC38, RT38), simplify="array")
  indelTable <- sapply(vcf, function(x) getIndels(x[isIndel(x)]), simplify="array")
  mnvTable <- sapply(vcf, function(x) getMNV(x), simplify="array")

  print("Trying to save ...")
  print("Writing")
  h5write(snvTensor, "SNV", file=paste0(outputPath, '.h5'))
  h5write(indelTable, "INDELS", file=paste0(outputPath, '.h5'))
  h5write(mnvTable, "MNV", file=paste0(outputPath, '.h5'))
} else { 
  c = 1
  for (inputPaths_chunk in chunk(inputPaths, 3)){
    print(inputPaths_chunk)
    vcf <- lapply(inputPaths_chunk, readVcfSave)
    snvTensor <- sapply(vcf, function(x) processVcf(x[isSNV(x)], TS38, EPI38, NUC38, RT38), simplify="array")
    indelTable <- sapply(vcf, function(x) getIndels(x[isIndel(x)]), simplify="array")
    mnvTable <- sapply(vcf, function(x) getMNV(x), simplify="array")

    print("Trying to save ...")
    print("Writing")
    h5write(snvTensor, "SNV", file=paste0(outputPath, '_chunk', c, '.h5'))
    h5write(indelTable, "INDELS", file=paste0(outputPath, '_chunk', c, '.h5'))
    h5write(mnvTable, "MNV", file=paste0(outputPath, '_chunk', c, '.h5'))
    c = c + 1
  }
}