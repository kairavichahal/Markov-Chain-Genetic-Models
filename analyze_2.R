load('chr_half_df.Rdata')  # this file in on google drive under 350_data

chr_half_df$orig_prop_single.a = chr_half_df$orig_prop_single.a / chr_half_df$leng
chr_half_df$orig_prop_single.c = chr_half_df$orig_prop_single.c / chr_half_df$leng
chr_half_df$orig_prop_single.g = chr_half_df$orig_prop_single.g / chr_half_df$leng
chr_half_df$orig_prop_single.t = chr_half_df$orig_prop_single.t / chr_half_df$leng

chr_half_df$sim_prop_single.a = chr_half_df$sim_prop_single.a / chr_half_df$leng
chr_half_df$sim_prop_single.c = chr_half_df$sim_prop_single.c / chr_half_df$leng
chr_half_df$sim_prop_single.g = chr_half_df$sim_prop_single.g / chr_half_df$leng
chr_half_df$sim_prop_single.t = chr_half_df$sim_prop_single.t / chr_half_df$leng

setwd("~/University/2013 Fall/36350/project/presentation")
# EDA
# univarite
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

# multivariate
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

# zoom in
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

# formal analysis
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


