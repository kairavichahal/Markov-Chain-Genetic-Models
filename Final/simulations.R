source('code.R')

# Run each simulation and save output to .Rdata file.
sim_chr_2f_1 = simulate_wrapper(chr.2f, k=1, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_2f_1, file='sim_chr_2f_1.Rdata')
rm(sim_chr_2f_1)

# (And so on for all combinations and chromosomes and ks...)