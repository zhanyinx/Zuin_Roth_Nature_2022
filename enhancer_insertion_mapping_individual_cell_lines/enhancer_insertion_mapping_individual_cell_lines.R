# Author: Pia Mach

library(BSgenome)
library(Biostrings)
library("BSgenome.Mmusculus.UCSC.mm9")
library(parallel)
library(stringr)

###GLOBAL OPTIONS#######################################################################################################################
# path to the directory with list of fasta files from Sanger sequencing (Microsynth AG) 
path = "/tungstenfs/scratch/ggiorget/Jessica/Mapping/3144766_Splinkerette_F07_dTAG_ei+eii_dCTCF_DDCTCF_08.03.21_3PRIME_PLATE_2/"
maplen = 24 # Length of read after restriction site/adaptor sequence that you would like to map to the genome
itr = 3 #If itr is set to 5, analysis for the 5'ITR splinkerette experiment is run, if set to 3, analysis for 3' ITR splinkerette is done

####TO NOT CHANGE ANYTHING AFTER HERE!!!################################################################################################
#Read in files with fasta format
files = list.files(path = path, pattern = ".fasta", full.names = TRUE)
reads = readDNAStringSet(filepath = files)
rev.com = reverseComplement(reads)
names(rev.com) = paste(names(rev.com), "rev_com", sep = "_")

#Read in SplitGFP and transgene sequences
splitGFP = readDNAStringSet("annotation/splitGFP.fa", format = "fasta")
cassette = readDNAStringSet("annotation/transgene.fa", format = "fasta")

if (itr==5) {
  
  #Define splinkerette adaptor sequences (end of it)
  restr.site = DNAString("CTAGGGTTAA")
  adaptor = DNAString("TAGTGGGATC")
  
  #Find adaptor sequences in reads
  matches_1 = vmatchPattern(restr.site, reads)
  matches_2 = vmatchPattern(adaptor, rev.com)
  
  #Trim reads to adaptor sequence
  for (p in 1:length(reads)){
    if (!is.null(unlist(startIndex(matches_1[p])))) {
      a = DNAStringSet(narrow(reads[p], start = unlist(start(matches_1[p]))[1]))
    reads[p] = a
    }
  }
  
  for (q in 1:length(rev.com)){
    if (!is.null(unlist(startIndex(matches_2[q])))) {
      b = DNAStringSet(narrow(rev.com[q], start = unlist(start(matches_2[q]))[1]))
      rev.com[q] = b
    }
  }
  
  #Filter on read-length (consider only reads that are longer than maplen+6 bp)
  # 6 is the number of bp in the adaptor excluding TTAA
  readlen = maplen+6
  reads = reads[width(reads)>readlen]
  
  #Removing all reads that do not start with the adaptor sequence and take maplen after the restriction site
  begin = narrow(reads, end=4) #depending on the length of the adaptor sequences
  keep = NULL
  for (i in 1:length(begin)) {
    x = ifelse(begin[i]==subseq(restr.site, start = 1, end = 4), TRUE, FALSE)
    keep = rbind(keep, x)
  }
  reads = reads[as.vector(keep)]
  reads = narrow(reads, start = 7, end = readlen) 
  
  # 10 is the length of the 
  readlen2 = maplen+10
  rev.com = rev.com[width(rev.com)>readlen2]
  begin2 = narrow(rev.com, end = 4) #depending on the length of the adaptor sequences
  keep2 = NULL
  for (i in 1:length(begin2)) {
    y = ifelse(begin2[i]==subseq(adaptor, start = 1, end = 4), TRUE, FALSE)
    keep2 = rbind(keep2, y)
  }
  rev.com = rev.com[as.vector(keep2)]
  rev.com = narrow(rev.com, start = 11, end = readlen2)
  
  
  #build the reverse complement of rev.com again, so all will be mapped on the same strand
  rev.com = reverseComplement(rev.com)
  #combine reverse complement and reads into one object
  seqs = c(reads, rev.com)

}

if(itr==3){
  #Define splinkerette adaptor sequences (end of it)
  restr.site = DNAString("CTTTCTAGGG")
  adaptor = DNAString("TAGTGGGATC")
  
  #Find adaptor sequences in reads
  matches_1 = vmatchPattern(restr.site, reads)
  matches_2 = vmatchPattern(adaptor, rev.com)
  
  #Trim reads to adaptor sequence
  for (p in 1:length(reads)){
    if (!is.null(unlist(startIndex(matches_1[p])))) {
      a = DNAStringSet(narrow(reads[p], start = unlist(start(matches_1[p]))[1]))
      reads[p] = a
    }
  }
  
  for (q in 1:length(rev.com)){
    if (!is.null(unlist(startIndex(matches_2[q])))) {
      b = DNAStringSet(narrow(rev.com[q], start = unlist(start(matches_2[q]))[1]))
      rev.com[q] = b
    }
  }
  
  #Filter on read-length (consider only reads that are longer than maplen+6 bp)
  readlen = maplen+10
  reads = reads[width(reads)>readlen]
  #Removing all reads that cannot be mapped and reducing the readlength
  begin = narrow(reads, end=4) #depending on the length of the adaptor sequences
  keep = NULL
  for (i in 1:length(begin)) {
    x = ifelse(begin[i]==subseq(restr.site, start = 1, end = 4), TRUE, FALSE)
    keep = rbind(keep, x)
  }
  reads = reads[as.vector(keep)]
  reads = narrow(reads, start = 11, end = readlen) 
  
  readlen2 = maplen+10
  rev.com = rev.com[width(rev.com)>readlen2]
  begin2 = narrow(rev.com, end = 4) #depending on the length of the adaptor sequences
  keep2 = NULL
  for (i in 1:length(begin2)) {
    y = ifelse(begin2[i]==subseq(adaptor, start = 1, end = 4), TRUE, FALSE)
    keep2 = rbind(keep2, y)
  }
  rev.com = rev.com[as.vector(keep2)]
  rev.com = narrow(rev.com, start = 11, end = readlen2)
  
  
  #build the reverse complement of rev.com again, so all will be mapped on the same strand
  rev.com = reverseComplement(rev.com)
  #combine reverse complement and reads into one object
  seqs = c(reads, rev.com)
  
}

#Map the sequences to the genome and write into GRanges object, if possible parallel process this
clObj = makeCluster(10)
#ranges = lapply(seqs, FUN = vmatchPattern, subject=BSgenome.Mmusculus.UCSC.mm9)
ranges = parLapply(cl = clObj, X=seqs, fun = vmatchPattern, subject=BSgenome.Mmusculus.UCSC.mm9, chunk.size = 10)
mappings = as(ranges, "GRangesList")
mappings = as.data.frame(mappings)

#Map for splitGFP
splitGFPmaps = parLapply(cl=clObj, X=reads, fun = vmatchPattern, subject = splitGFP, chunk.size =10)
splitGFPmaps = as(splitGFPmaps, "GRangesList")
splitGFPmaps = as.data.frame(splitGFPmaps)

#Map for cassette
cassettemaps = parLapply(cl=clObj, X=seqs, fun = vmatchPattern, subject = cassette, chunk.size =10)
cassettemaps = as(cassettemaps, "GRangesList")
cassettemaps = as.data.frame(cassettemaps)

#Bind mappings together to include splitGFP, cassette
maps = rbind(mappings, splitGFPmaps, cassettemaps)

#Make the table prettier
maps$PCR_end = maps$group_name
maps$clone= maps$group_name
maps$PCR_end[grep("_rev_com", maps$PCR_end)] <- rep(c("3PRIME"), length(grep("_rev_com", maps$PCR_end)))
maps$PCR_end[grep("_3PRIME", maps$PCR_end)] <- rep(c("5PRIME"), length(grep("_3PRIME", maps$PCR_end)))
maps$clone = str_sub(maps$group_name, 1, 3)
maps$splink_ITR = rep(paste(itr, "ITR", sep = ""), nrow(maps))

#Write mappings into a csv file
write.csv(maps, file = paste(path, "mappings_maplen=", maplen, "_", itr, "ITR_splitGFP_cassette.csv", sep = ""))
