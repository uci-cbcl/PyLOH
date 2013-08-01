library(BICseq)

bicseq <- BICseq(sample = 'tumor_recal.bam', reference = 'normal_recal.bam')

segs <- getBICseg(object = bicseq, bin = 100, lambda = 2, winSize = 200, quant = 0.95, mult = 1)

seg.summary <- BICseq:::getSummary(segs, correction=TRUE)

write.table(seg.summary, file = './segments.BICseq', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)







