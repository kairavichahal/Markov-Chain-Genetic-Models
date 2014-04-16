################################################################################
# ORGANIZE DATA ################################################################
################################################################################


# Only works if only simulation .Rdata files (nothing else) are in ~/sim.
to_load = list.files()
for (i in 1:length(to_load)){
  load(to_load[i])
}

# Create a data frame containing statistics to be analyzed.

# 2f_1
chr_2f_1_bw = matrix(data=unlist(sim_chr_2f_1[3,]), ncol=10, byrow=T)
colnames(chr_2f_1_bw) = names(unlist(sim_chr_2f_1[3,]))[1:10]
result_chr_2f_1 = cbind(unlist(sim_chr_2f_1[1,]),
                        chr_2f_1_bw,
                        unlist(sim_chr_2f_1[4,]),
                        unlist(sim_chr_2f_1[5,]))
colnames(result_chr_2f_1)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_2f_1) = paste('2f_1_', 1:10, sep='')

# 2f_2
chr_2f_2_bw = matrix(data=unlist(sim_chr_2f_2[3,]), ncol=10, byrow=T)
colnames(chr_2f_2_bw) = names(unlist(sim_chr_2f_2[3,]))[1:10]
result_chr_2f_2 = cbind(unlist(sim_chr_2f_2[1,]),
                        chr_2f_2_bw,
                        unlist(sim_chr_2f_2[4,]),
                        unlist(sim_chr_2f_2[5,]))
colnames(result_chr_2f_2)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_2f_2) = paste('2f_2_', 1:10, sep='')

# 2f_3
chr_2f_3_bw = matrix(data=unlist(sim_chr_2f_3[3,]), ncol=10, byrow=T)
colnames(chr_2f_3_bw) = names(unlist(sim_chr_2f_3[3,]))[1:10]
result_chr_2f_3 = cbind(unlist(sim_chr_2f_3[1,]),
                        chr_2f_3_bw,
                        unlist(sim_chr_2f_3[4,]),
                        unlist(sim_chr_2f_3[5,]))
colnames(result_chr_2f_3)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_2f_3) = paste('2f_3_', 1:10, sep='')

# m_1
chr_m_1_bw = matrix(data=unlist(sim_chr_m_1[3,]), ncol=10, byrow=T)
colnames(chr_m_1_bw) = names(unlist(sim_chr_m_1[3,]))[1:10]
result_chr_m_1 = cbind(unlist(sim_chr_m_1[1,]),
                       chr_m_1_bw,
                       unlist(sim_chr_m_1[4,]),
                       unlist(sim_chr_m_1[5,]))
colnames(result_chr_m_1)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_m_1) = paste('m_1_', 1:10, sep='')

# m_2
chr_m_2_bw = matrix(data=unlist(sim_chr_m_2[3,]), ncol=10, byrow=T)
colnames(chr_m_2_bw) = names(unlist(sim_chr_m_2[3,]))[1:10]
result_chr_m_2 = cbind(unlist(sim_chr_m_2[1,]),
                       chr_m_2_bw,
                       unlist(sim_chr_m_2[4,]),
                       unlist(sim_chr_m_2[5,]))
colnames(result_chr_m_2)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_m_2) = paste('m_2_', 1:10, sep='')

# m_3
chr_m_3_bw = matrix(data=unlist(sim_chr_m_3[3,]), ncol=10, byrow=T)
colnames(chr_m_3_bw) = names(unlist(sim_chr_m_3[3,]))[1:10]
result_chr_m_3 = cbind(unlist(sim_chr_m_3[1,]),
                       chr_m_3_bw,
                       unlist(sim_chr_m_3[4,]),
                       unlist(sim_chr_m_3[5,]))
colnames(result_chr_m_3)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_m_3) = paste('m_3_', 1:10, sep='')

# bf_1
chr_bf_1_bw = matrix(data=unlist(sim_chr_bf_k1[3,]), ncol=10, byrow=T)
colnames(chr_bf_1_bw) = names(unlist(sim_chr_bf_k1[3,]))[1:10]
result_chr_bf_1 = cbind(unlist(sim_chr_bf_k1[1,]),
                        chr_bf_1_bw,
                        unlist(sim_chr_bf_k1[4,]),
                        unlist(sim_chr_bf_k1[5,]))
colnames(result_chr_bf_1)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_bf_1) = paste('bf_1_', 1:10, sep='')

# bf_2
chr_bf_2_bw = matrix(data=unlist(sim_chr_bf_k2[3,]), ncol=10, byrow=T)
colnames(chr_bf_2_bw) = names(unlist(sim_chr_bf_k2[3,]))[1:10]
result_chr_bf_2 = cbind(unlist(sim_chr_bf_k2[1,]),
                        chr_bf_2_bw,
                        unlist(sim_chr_bf_k2[4,]),
                        unlist(sim_chr_bf_k2[5,]))
colnames(result_chr_bf_2)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_bf_2) = paste('bf_2_', 1:10, sep='')

# bf_3
chr_bf_3_bw = matrix(data=unlist(sim_chr_bf_k3[3,]), ncol=10, byrow=T)
colnames(chr_bf_3_bw) = names(unlist(sim_chr_bf_k3[3,]))[1:10]
result_chr_bf_3 = cbind(unlist(sim_chr_bf_k3[1,]),
                        chr_bf_3_bw,
                        unlist(sim_chr_bf_k3[4,]),
                        unlist(sim_chr_bf_k3[5,]))
colnames(result_chr_bf_3)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_bf_3) = paste('bf_3_', 1:10, sep='')

# r_1
chr_r_1_bw = matrix(data=unlist(sim_chr_r_k1[3,]), ncol=10, byrow=T)
colnames(chr_r_1_bw) = names(unlist(sim_chr_r_k1[3,]))[1:10]
result_chr_r_1 = cbind(unlist(sim_chr_r_k1[1,]),
                       chr_r_1_bw,
                       unlist(sim_chr_r_k1[4,]),
                       unlist(sim_chr_r_k1[5,]))
colnames(result_chr_r_1)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_r_1) = paste('r_1_', 1:10, sep='')

# r_2
chr_r_2_bw = matrix(data=unlist(sim_chr_r_k2[3,]), ncol=10, byrow=T)
colnames(chr_r_2_bw) = names(unlist(sim_chr_r_k2[3,]))[1:10]
result_chr_r_2 = cbind(unlist(sim_chr_r_k2[1,]),
                       chr_r_2_bw,
                       unlist(sim_chr_r_k2[4,]),
                       unlist(sim_chr_r_k2[5,]))
colnames(result_chr_r_2)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_r_2) = paste('r_2_', 1:10, sep='')

# r_3
chr_r_3_bw = matrix(data=unlist(sim_chr_r_k3[3,]), ncol=10, byrow=T)
colnames(chr_r_3_bw) = names(unlist(sim_chr_r_k3[3,]))[1:10]
result_chr_r_3 = cbind(unlist(sim_chr_r_k3[1,]),
                       chr_r_3_bw,
                       unlist(sim_chr_r_k3[4,]),
                       unlist(sim_chr_r_k3[5,]))
colnames(result_chr_r_3)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_r_3) = paste('r_3_', 1:10, sep='') 

# 3f_1
chr_3f_1_bw = matrix(data=unlist(sim_chr_3f_k1[3,]), ncol=10, byrow=T)
colnames(chr_3f_1_bw) = names(unlist(sim_chr_3f_k1[3,]))[1:10]
result_chr_3f_1 = cbind(unlist(sim_chr_3f_k1[1,]),
                        chr_3f_1_bw,
                        unlist(sim_chr_3f_k1[4,]),
                        unlist(sim_chr_3f_k1[5,]))
colnames(result_chr_3f_1)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_3f_1) = paste('3f_1_', 1:10, sep='') 

# 3f_2
chr_3f_2_bw = matrix(data=unlist(sim_chr_3f_k2[3,]), ncol=10, byrow=T)
colnames(chr_3f_2_bw) = names(unlist(sim_chr_3f_k2[3,]))[1:10]
result_chr_3f_2 = cbind(unlist(sim_chr_3f_k2[1,]),
                        chr_3f_2_bw,
                        unlist(sim_chr_3f_k2[4,]),
                        unlist(sim_chr_3f_k2[5,]))
colnames(result_chr_3f_2)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_3f_2) = paste('3f_2_', 1:10, sep='') 

# 3f_3
chr_3f_3_bw = matrix(data=unlist(sim_chr_3f_k3[3,]), ncol=10, byrow=T)
colnames(chr_3f_3_bw) = names(unlist(sim_chr_3f_k3[3,]))[1:10]
result_chr_3f_3 = cbind(unlist(sim_chr_3f_k3[1,]),
                        chr_3f_3_bw,
                        unlist(sim_chr_3f_k3[4,]),
                        unlist(sim_chr_3f_k3[5,]))
colnames(result_chr_3f_3)[c(1,12,13)] = c('k', 'mtx_cf', 'pw_cf')
rownames(result_chr_3f_3) = paste('3f_3_', 1:10, sep='') 

chr_half = rbind(result_chr_m_1,  result_chr_m_2,  result_chr_m_3, 
                 result_chr_2f_1, result_chr_2f_2, result_chr_2f_3, 
                 result_chr_3f_1, result_chr_3f_2, result_chr_3f_3, 
                 result_chr_bf_1, result_chr_bf_2, result_chr_bf_3, 
                 result_chr_r_1,  result_chr_r_2,  result_chr_r_3)
chr_half_df = as.data.frame(chr_half)
chr_half_df$chr = rep(c('m', '2f', '3f', 'bf', 'r'), each=30)
chr_half_df$leng = rep(chr_leng[c('M', '2F', '3F', 'BF', 'R')], each=30)
setwd("~/University/2013 Fall/36350/project")
save(chr_half_df, file='chr_half_df.Rdata')

################################################################################
# ANALYZE DATA #################################################################
################################################################################

load('chr_half_df.Rdata')

chr_half_df$orig_prop_single.a = chr_half_df$orig_prop_single.a / chr_half_df$leng
chr_half_df$orig_prop_single.c = chr_half_df$orig_prop_single.c / chr_half_df$leng
chr_half_df$orig_prop_single.g = chr_half_df$orig_prop_single.g / chr_half_df$leng
chr_half_df$orig_prop_single.t = chr_half_df$orig_prop_single.t / chr_half_df$leng

chr_half_df$sim_prop_single.a = chr_half_df$sim_prop_single.a / chr_half_df$leng
chr_half_df$sim_prop_single.c = chr_half_df$sim_prop_single.c / chr_half_df$leng
chr_half_df$sim_prop_single.g = chr_half_df$sim_prop_single.g / chr_half_df$leng
chr_half_df$sim_prop_single.t = chr_half_df$sim_prop_single.t / chr_half_df$leng

# Univariate EDA
png('plot_uni.png', width=3000, height=1000, res=150)
par(mfrow=c(1,3))
hist(chr_half_df$match_prop, xlab='Proportion of Exact Matches', freq=F,
     main='Histogram of Proportion of Exact Matches', col='skyblue',
     cex.main=2, cex.lab=2)
hist(chr_half_df$mtx_cf, xlab='Theta for Comparing Transition Matrices', freq=F,
     main='Histogram of Theta for Comparing Transition Matrices', col='tomato2',
     cex.main=2, cex.lab=2)
hist(chr_half_df$pw_cf, breaks=15, xlab='Pairwise Global Alignment Score', freq=F,
     main='Histogram of Pairwise Global Alignment Score', col='darkgoldenrod2',
     cex.main=2, cex.lab=2)
dev.off()

# Multivariate EDA
png('plot_prop.png', width=2000, height=2000, res=150)
par(mfrow=c(2,2))
boxplot(chr_half_df$sim_prop_single.a ~ chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Marginal Proportion of As', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Marginal Proportion of As in Simulated Sequences \n by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
points(1:15, chr_half_df$orig_prop_single.a[seq(1,150,by=10)[order(chr_half_df$chr[seq(1,150,by=10)])]], 
       col='red', pch=18)
legend('topleft', 'Original Sequence', pch=18, col='red')

boxplot(chr_half_df$sim_prop_single.c ~ chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Marginal Proportion of Cs', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Marginal Proportion of Cs in Simulated Sequences \n by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
points(1:15, chr_half_df$orig_prop_single.c[seq(1,150,by=10)[order(chr_half_df$chr[seq(1,150,by=10)])]], 
       col='red', pch=18)
legend('bottomleft', 'Original Sequence', pch=18, col='red')

boxplot(chr_half_df$sim_prop_single.g ~ chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Marginal Proportion of Gs', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Marginal Proportion of Gs in Simulated Sequences \n by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
points(1:15, chr_half_df$orig_prop_single.g[seq(1,150,by=10)[order(chr_half_df$chr[seq(1,150,by=10)])]], 
       col='red', pch=18)
legend('topleft', 'Original Sequence', pch=18, col='red')

boxplot(chr_half_df$sim_prop_single.t ~ chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Marginal Proportion of Ts', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Marginal Proportion of Ts in Simulated Sequences \n by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
points(1:15, chr_half_df$orig_prop_single.t[seq(1,150,by=10)[order(chr_half_df$chr[seq(1,150,by=10)])]], 
       col='red', pch=18)
legend('bottomleft', 'Original Sequence', pch=18, col='red')
dev.off()

png('plot_exact.png', width=1000, height=1000, res=150)
boxplot(chr_half_df$match_prop~chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Proportion of Base-wise Exact Matches', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Proportions of Base-wise Exact Matches in \n Simulated Sequences by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
dev.off()

png('plot_mtx.png', width=1000, height=1000, res=150)
boxplot(chr_half_df$mtx_cf~chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Theta for Comparing Transition Matrices', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Thetas for Comparing Transition Matrices in \n Simulated Sequences by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
dev.off()

png('plot_pw.png', width=1000, height=1000, res=150)
boxplot(chr_half_df$pw_cf~chr_half_df$k+chr_half_df$chr,
        xlab='Order k & Chromosome', ylab='Pairwise Global Alignment Scores', 
        border=rep(c('skyblue','hotpink','darkgoldenrod2','black', 'forestgreen'), each=3),
        main='Pairwise Global Alignment Scores in \n Simulated Sequences by Order K & Chromosome', 
        cex.main=1.5, cex.lab=1.5)
dev.off()

png('plot_pw_zoom.png', width=1000, height=1000, res=150)
boxplot(chr_half_df$pw_cf[which(chr_half_df$k==1&chr_half_df$chr=='2f')], 
        chr_half_df$pw_cf[which(chr_half_df$k==2&chr_half_df$chr=='2f')],
        chr_half_df$pw_cf[which(chr_half_df$k==3&chr_half_df$chr=='2f')],
        xlab='Order k', ylab='Pairwise Global Alignment Scores', 
        border='skyblue',
        main='Pairwise Global Alignment Scores in \n Simulated Sequences of Chromosome 2f by Order K', 
        cex.main=1.5, cex.lab=1.5)
axis(1, 1:3)
dev.off()

# Formal Analysis - Regression
lm.1 = lm(chr_half_df$match_prop ~ chr_half_df$chr + chr_half_df$k)
lm.2 = lm(chr_half_df$mtx_cf ~ chr_half_df$leng + chr_half_df$k)
lm.3 = lm(chr_half_df$mtx_cf ~ chr_half_df$chr + chr_half_df$k)
lm.4 = lm(chr_half_df$pw_cf ~ chr_half_df$leng + chr_half_df$k)
lm.5 = lm(chr_half_df$pw_cf ~ chr_half_df$chr + chr_half_df$k)
summary(lm.1)
summary(lm.2)
summary(lm.3)
summary(lm.4)
summary(lm.5)