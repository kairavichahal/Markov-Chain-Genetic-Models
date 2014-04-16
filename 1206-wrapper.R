source('1206.R')

################
# CHROMOSOME 5 #
################
sim_chr_5_1 = simulate_wrapper(chr.5, k=1, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_5_1, file='sim_chr_5_1.Rdata'); print(paste("Saved sim_chr_5_1.Rdata!"))
rm(sim_chr_5_1)

sim_chr_5_2 = simulate_wrapper(chr.5, k=2, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_5_2, file='sim_chr_5_2.Rdata'); print(paste("Saved sim_chr_5_2.Rdata!"))
rm(sim_chr_5_2)

sim_chr_5_3 = simulate_wrapper(chr.5, k=3, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_5_3, file='sim_chr_5_3.Rdata'); print(paste("Saved sim_chr_5_3.Rdata!"))
rm(sim_chr_5_3)

################
# CHROMOSOME 6 #
################
sim_chr_6_1 = simulate_wrapper(chr.6, k=1, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_6_1, file='sim_chr_6_1.Rdata'); print(paste("Saved sim_chr_6_1.Rdata!"))
rm(sim_chr_6_1)

sim_chr_6_2 = simulate_wrapper(chr.6, k=2, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_6_2, file='sim_chr_6_2.Rdata'); print(paste("Saved sim_chr_6_2.Rdata!"))
rm(sim_chr_6_2)

sim_chr_6_3 = simulate_wrapper(chr.6, k=3, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
save(sim_chr_6_3, file='sim_chr_6_3.Rdata'); print(paste("Saved sim_chr_6_3.Rdata!"))
rm(sim_chr_6_3)