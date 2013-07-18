library(ExomeCNV)

normal = read.delim("http://genome.ucla.edu/~fah/ExomeCNV/data/normal.baf.txt", header=T)
tumor = read.delim("http://genome.ucla.edu/~fah/ExomeCNV/data/tumor.baf.txt", header=T)

eLOH = LOH.analyze(normal, tumor, alpha=0.05, method="two.sample.fisher")

eLOH_out = cbind(eLOH[1], eLOH[2], eLOH[8])

write.table(eLOH_out, file = "test.LOH", sep = "\t", row.names = FALSE, quote = FALSE)

