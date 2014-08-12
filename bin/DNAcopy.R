library(DNAcopy)

data = read.delim("./DNAcopy.bed")

log2ratio = as.vector(as.matrix(data[6]))
chrom = as.vector(as.matrix(data[1]))
pos = as.vector(as.matrix(data[7]))

CNA.object <- CNA(log2ratio, chrom, pos, data.type="logratio", sampleid="ID", presorted=FALSE)
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1, min.width=5, alpha=0.001)
write.table(segment.smoothed.CNA.object$output[3:6], file = 'segments.DNAcopy', col.names = F, row.names = F, quote = F, sep = "\t")

