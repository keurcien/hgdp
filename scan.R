require(data.table)
require(pcadapt)
setwd("~/Data/hgdp/")

pop <- as.character(read.table("hgdp_final.pop")[, 1])
map <- fread("hgdp_final.map", data.table = FALSE)
imputed.geno <- impute.pcadapt("hgdp_final.pcadapt", pop = pop, skip.return = TRUE)
map <- map[-imputed.geno$ix, ]

write.table(imputed.geno$x, "hgdp_final_imputed.pcadapt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(map, "hgdp_final_imputed.map", col.names = FALSE, row.names = FALSE, quote = FALSE)

geno <- as.matrix(fread("hgdp_final_imputed.pcadapt", data.table = FALSE))
map <-  fread("hgdp_final_imputed.map", data.table = FALSE)
chr <- map$V1
stat <- scan.intro(geno, K = 1, ancstrl.1 = "Han", ancstrl.2 = "French", admxd = "Uygur", 
                   pop = pop, ploidy = 2, window.size = 100000, chr.info = chr, map = map$V4)
