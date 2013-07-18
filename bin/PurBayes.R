library(PurBayes)

somatics = read.delim("./somatics.VarScan.somatic", header=F)

Y = as.vector(as.matrix(somatics[10]))
N = as.vector(as.matrix(somatics[9] + somatics[10]))
PB = PurBayes(N, Y, M=NULL, Z=NULL, pop.max=5, prior=NULL, burn.in=50000, n.post = 10000, fn.jags = "PB.jags", plot = FALSE)
results = summary(PB)
write.table(results$purity, file = 'results.PurBayes', col.names = T, quote = F)

