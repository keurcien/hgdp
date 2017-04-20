library(dplyr)
library(data.table)


### Extracting Uygur, French and Han ###
setwd("~/Data/hgdp/")
geno <- fread("HGDP_FinalReport_Forward.txt", data.table = FALSE, drop = 1)
map <- fread("HGDP_Map.txt", data.table = FALSE)
colnames(map) <- c("SNP", "Chr", "Position")

samp <- as.character(read.table("HGDP_SampleList.txt", header = FALSE)[,1])
ceph <- read.csv2("HGDP-CEPH-ID_populations.csv", header = TRUE)
pop.ceph <- as.character(ceph$population)
ceph <- ceph[pop.ceph %in% c("Han", "Uygur", "French"), ]
pop.ceph <- as.character(ceph$population)
ind.ceph <- as.character(ceph$CEPH.ID)

subs <- names(geno) %in% ind.ceph
geno <- geno[, subs]
geno <- cbind(map, geno)
chr.autho <- 1:22
geno <- geno[geno$Chr %in% chr.autho, ]
geno <- arrange(geno, Chr, Position)

write.table(geno, "HGDP_Uygur.txt", col.names = TRUE, row.names = FALSE, 
            quote = FALSE)
########################################

geno <- fread("HGDP_Uygur.txt")
info <- geno[, 1:3]
geno <- geno[, -(1:3)]

filename <- "hgdp.ped"
nIND <- ncol(geno)
nSNP <- nrow(geno)

x <- vector(mode = "numeric", length = (2 * nSNP + 6))
pop <- vector(mode = "character", length = nIND)
ind.i <- names(geno)[1]  
x[1] <- pop.ceph[which(ind.ceph == ind.i)]
x[2] <- ind.i
pop[1] <- pop.ceph[which(ind.ceph == ind.i)]
tmp <- paste(geno[, 1], collapse = "")
tmp.bis <- strsplit(tmp, split = "")[[1]]
tmp.bis[which(tmp.bis == "-")] <- "0"
x[-(1:6)] <- tmp.bis

cat(paste(x, collapse = " "), file = filename, sep = "\n", append = FALSE)

for (i in 2:nIND){
  ind.i <- names(geno)[i]  
  x[0] <- pop.ceph[which(ind.ceph == ind.i)]
  x[1] <- ind.i
  pop[i] <- pop.ceph[which(ind.ceph == ind.i)]
  tmp <- paste(geno[, i], collapse = "")
  tmp.bis <- strsplit(tmp, split = "")[[1]]
  tmp.bis[which(tmp.bis == "-")] <- "0"
  x[-(1:6)] <- tmp.bis
  cat(i, "\n")
  cat(paste(x, collapse = " "), file = filename, sep = "\n", append = TRUE)
}

x <- fread("hgdp.pcadapt", data.table = FALSE)
check.missing <- apply(x, MARGIN = 1, FUN = function(h){sum(h == 9)})
keep.snp <- check.missing < 83
x <- x[keep.snp, ]
info <- info[keep.snp, ]
gen.dist <- vector(mode = "numeric", length = length(info$Chr))
info.map <- data.frame(info$Chr, info$SNP, gen.dist, info$Position)

write.table(x, "hgdp_final.pcadapt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(info.map, "hgdp_final.map", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(pop, "hgdp_final.pop", col.names = FALSE, row.names = FALSE, quote = FALSE)


