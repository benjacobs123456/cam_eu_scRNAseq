args = commandArgs(TRUE)


cam = read.table(args[1])
tum = read.table(args[2])
sites = read.table(args[3])
outfile = args[4]


cam$Comb1 = paste0(cam$V1,"_",cam$V2)
tum$Comb1 = paste0(tum$V1,"_",tum$V2)
sites$Comb1 = paste0(sites$V1,"_",sites$V2)

# Only keep sites that are shared
cam = cam[cam$Comb1 %in% sites$Comb1 == T,]
tum = tum[tum$Comb1 %in% sites$Comb1 == T,]

# Remove indesl
cam = cam[nchar(cam$V3) == 1 & nchar(cam$V4) == T,]
tum = tum[nchar(tum$V3) == 1 & nchar(tum$V4) == T,]

cam$Comb2 = paste0(cam$Comb1,"_",cam$V3,"_",cam$V4)
cam$Comb3 = paste0(cam$Comb1,"_",cam$V4,"_",cam$V3)
tum$Comb2 = paste0(tum$Comb1,"_",tum$V3,"_",tum$V4)
tum$Comb3 = paste0(tum$Comb1,"_",tum$V4,"_",tum$V3)


print("SNPs with strand issues:")
print(nrow(tum[tum$Comb3 %in% cam$Comb2 == T,]))