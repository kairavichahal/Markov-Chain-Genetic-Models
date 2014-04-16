test.chr = all_seq_wo_n[c(2, 8:10)]
for (chr in 1:4) {
  for (k in 1:3) {
    sim_chr = simulate_wrapper(test.chr[chr], k=k, m=+2, x=-2, g_open=-5, g_extend=-2, num_sim=10, incl_sim_seq=T)
    write(sim_chr, paste("chr_", chr, "_k_", k, ".txt", sep=""))
  }
}