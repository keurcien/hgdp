library(dplyr)
library(data.table)

setwd("~/Documents/thesis/datasets/hgdp")

geno <- fread("HGDP_FinalReport_Forward.txt")
map <- fread("HGDP_Map.txt")
colnames(map) <- c("SNP", "Chr", "Position")
geno <- cbind(map, geno)
arrange(geno, Chr, Position)
chr.autho <- 1:22


write.table(geno, "HGDP_Sorted.txt", col.names = TRUE, row.names = FALSE, 
            quote = FALSE)

df <- fread("HGDP_Sorted.txt")

df <- df[df$Chr %in% chr.autho, ]
info <- df[,1:4]
df <- df[,-(1:4)]
dstan <- as.character(read.table("HGDP_SampleList.txt", header = FALSE)[,1])
dceph <- read.csv2("HGDP-CEPH-ID_populations.csv", header = TRUE)
pop.ceph <- as.character(dceph$population)
dceph <- dceph[pop.ceph %in% c("Han", "Uygur", "French"), ]
pop.ceph <- as.character(dceph$population)
ind.ceph <- as.character(dceph$CEPH.ID)
coldf <- names(df) %in% ind.ceph
df <- df[,..coldf]
df <- as.data.frame(df)


filename <- "hgdp.ped"
nIND <- ncol(df)
nSNP <- nrow(df)
x <- vector(mode = "numeric", length = (2 * nSNP + 6))
ind.i <- names(df)[1]  
x[0] <- pop.ceph[which(ind.ceph == ind.i)]
x[1] <- ind.i
tmp <- paste(df[, 1], collapse = "")
tmp.bis <- strsplit(tmp, split = "")[[1]]
tmp.bis[which(tmp.bis == "-")] <- "0"
x[-(1:6)] <- tmp.bis
cat(1, "\n")
cat(paste(x, collapse = " "), file = filename, sep = "\n", append = FALSE)
for (i in 2:nIND){
  ind.i <- names(df)[i]  
  x[0] <- pop.ceph[which(ind.ceph == ind.i)]
  x[1] <- ind.i
  tmp <- paste(df[, i], collapse = "")
  tmp.bis <- strsplit(tmp, split = "")[[1]]
  tmp.bis[which(tmp.bis == "-")] <- "0"
  x[-(1:6)] <- tmp.bis
  cat(i, "\n")
  cat(paste(x, collapse = " "), file = filename, sep = "\n", append = TRUE)
}


