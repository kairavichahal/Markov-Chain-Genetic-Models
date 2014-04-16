k1 = chr[which(chr$k == 1),c(1,2,3,12,13,14)]
k2 = chr[which(chr$k == 2),c(1,2,3,12,13,14)]
k3 = chr[which(chr$k == 3),c(1,2,3,12,13,14)]

chr.m = chr[which(chr$chromosome == "m"),c(1,2,3,12,13)]
chr.2f = chr[which(chr$chromosome == "2f"),c(1,2,3,12,13)]
chr.3f = chr[which(chr$chromosome == "3f"),c(1,2,3,12,13)]
chr.bf = chr[which(chr$chromosome == "bf"),c(1,2,3,12,13)]
chr.r = chr[which(chr$chromosome == "r"),c(1,2,3,12,13)]

chromosome = rep(c("m", "2f", "3f", "bf", "r"), each=30))
length = rep(c(55564, 161967, 16660, 75732, 85150), each=30)