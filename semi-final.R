### Don't forget to uncomment the chromosomes that you'll be working on 
### and comment out the ones that you won't

# 36-350 Final Project
# Markov Chain Genetic Models
# Kairavi Chahal, Tony Yang, Julian Zhou

########################################################
############### LOAD REQUIRED PACKAGES #################
########################################################
# install.packages('seqinr')
# source("http://bioconductor.org/biocLite.R")
# biocLite('Biostrings')
library(seqinr)
library(Biostrings)

########################################################
#################### READ SEQUENCE #####################
########################################################

# Read in FASTA file
raw_data = read.fasta(file = "dicty_chromosomal")

# use regex to retrieve names of chromosomes
retrieve_chr_name_length = function(fasta){
  # Input: data read in by read.fasta
  # Output: a character vector containing chromosome names and a vector of lengths
  chr_name_pattern = ":.+position"
  have_chr_name = grep(x=getAnnot(fasta), pattern=chr_name_pattern)
  chr_name_matches = gregexpr(pattern=chr_name_pattern, text=getAnnot(fasta)[have_chr_name])
  chr_name_pre = regmatches(x=unlist(getAnnot(fasta)[have_chr_name]), m=chr_name_matches)
  chr_name_pre = do.call(c, chr_name_pre)
  chr_name_pre = strsplit(chr_name_pre, ' ')
  chr_name_pre = do.call(rbind, chr_name_pre)
  chr_name = chr_name_pre[, 2] 
  
  chr_leng_pattern = "to.+"
  have_chr_leng = grep(x=getAnnot(fasta), pattern=chr_leng_pattern)
  chr_leng_matches = gregexpr(pattern=chr_leng_pattern, text=getAnnot(fasta)[have_chr_leng])
  chr_leng_pre = regmatches(x=unlist(getAnnot(fasta)[have_chr_leng]), m=chr_leng_matches)
  chr_leng_pre = do.call(c, chr_leng_pre)
  chr_leng_pre = strsplit(chr_leng_pre, ' ')
  chr_leng_pre = do.call(rbind, chr_leng_pre)
  chr_leng = as.numeric(chr_leng_pre[, 2])
  return(list(chr_name = chr_name, chr_leng = chr_leng))
}

chr_name = retrieve_chr_name_length(raw_data)$chr_name
num_chr = length(chr_name)
chr_leng = retrieve_chr_name_length(raw_data)$chr_leng

# verify lengths of chromosomes against function from package seqinr
all.equal(getLength(raw_data), chr_leng)

# Get sequences
# yields a list of length 11
all_seq = getSequence(raw_data)
rm(raw_data)

# Replace the Ns in sequence with a, c, g, t with prob. 0.25
replace_n = function(chr_vec){
  #Input:a vector chromosome sequence with possible unknown 'n' patterns
  #Output:a vector chromosome sequence with 'n' replaced by the four amino acids in uniform distr
  idx = which(chr_vec == "n")
  chr_vec[idx] = sample(c("a", "c", "g", "t"), length(idx), replace=T)
  return(chr_vec)
}


#all_seq_wo_n = lapply(all_seq, replace_n)
#save(all_seq_wo_n, file='all_seq_wo_n.Rdata')
load('all_seq_wo_n.Rdata')
rm(all_seq)

# If you want separate variables in each chromosome
# chr.m  = all_seq_wo_n[[1]]
chr.2f = all_seq_wo_n[[2]]
# chr.3f = all_seq_wo_n[[3]]
# chr.bf = all_seq_wo_n[[4]]
# chr.1  = all_seq_wo_n[[5]]
# chr.2  = all_seq_wo_n[[6]]
# chr.3  = all_seq_wo_n[[7]]
chr.4  = all_seq_wo_n[[8]]
chr.5  = all_seq_wo_n[[9]]
chr.6  = all_seq_wo_n[[10]]
# chr.r  = all_seq_wo_n[[11]]

# convert ATGC's in a sequence to 1234
seq_convert = function(seq){
  # Input: character vector
  # output: numeric vector
  stopifnot(is.vector(seq))
  seq[which(seq=='a')] = 1
  seq[which(seq=='c')] = 2
  seq[which(seq=='g')] = 3
  seq[which(seq=='t')] = 4
  return(as.numeric(seq))
}

# test seq_convert
test_seq_convert = function(){
  test = seq_convert(c('a','t','g','c','a','a','c'))
  right = c(1,4,3,2,1,1,2)
  return(all.equal(test, right))
}

# actual conversion
#all_seq_wo_n_num = lapply(all_seq_wo_n, seq_convert)
#save(all_seq_wo_n_num, file='all_seq_wo_n_num.Rdata')
load('all_seq_wo_n_num.Rdata')
#rm(all_seq_wo_n)

# doesn't make sense to concatenate all chromosomes to create a 'genome' and compute its transition matrix

########################################################
############## MISCELLANEOUS FUNCTIONS #################
########################################################

# character vector to string
char2str = function(char_vec){
  return(paste(char_vec, collapse=''))
}

# string to character vector
str2char = function(string){
  return(unlist(strsplit(string, '')))
}

# Function below simply return number of seconds taken
# to execute func with arguments ...
time_function = function(func, ...) {
  start_time = proc.time()
  func(...)
  stop_time = proc.time() - start_time
  return(stop_time[3])
}

########################################################
############## COMPUTE PROBABILITY MATRIX ##############
########################################################

# Convert raw count matrix into transition matrix
prob_of = function(mtx_row) {
  if (sum(mtx_row)>0){
    mtx_row = mtx_row/sum(mtx_row)
  } else {
    # Replace with 0.25 if this sub-sequence does not exist in the original
    mtx_row=rep(1/4, 4)
  }
  return(mtx_row)
}

# Helper function to create all combinations of ACGT
get_combinations = function(k) {
  l = rep(list(c('a','c','g','t')), k)
  grid = as.matrix(expand.grid(l))
  rownames = sort(apply(grid, 1, paste, collapse=""))
  stopifnot(length(rownames)==4^k)
  return(rownames)
}

# Helper function for compute_markov():
get_substring = function(i, k, character_vector) {
  # Input: index, order of Markov and character vector
  # Output: substring containing (i-k)th to ith characters of the sequence
  substring = char2str(character_vector[(i-k):i])
  return(substring)
}


compute_markov = function(sequence, k) {
  # Input: Sequence as a character vector and what order Markov we want.
  # Output: Markov probability matrix of dimension 4^k by 4.
  
  print("Computing Markov...")
  
  sub_seqs = sapply((k+1):length(sequence), get_substring, k=k, character_vector=sequence)
  sub_seqs_tab = table(sub_seqs)
  markov = matrix(0, ncol=4, nrow=4^k, dimnames=list(get_combinations(k), c("a", "c", "g", "t")))
  sub_seqs_rows = substr(rownames(sub_seqs_tab), 1, k)
  sub_seqs_cols = substr(rownames(sub_seqs_tab), k+1, k+1)
  mat = cbind(sub_seqs_rows, sub_seqs_cols, sub_seqs_tab)
  
  row_idx = match(mat[,1], rownames(markov))
  col_idx = match(mat[,2], colnames(markov))
  markov[cbind(row_idx, col_idx)] = as.numeric(mat[, 3])
  
  print("Done computing Markov!")
  return(t(apply(markov, 1, prob_of)))
}

########################################################
##################### SIMULATION #######################
########################################################

# Generate next base based on prior sequence
generate_next_base = function(prior_seq, markov, k){
  # Input: prior sequence (string), probability matrix, order k
  # Output: a single character, which will be the next base
  
  # find the matching row number given a prior sequence
  row_num = which(rownames(markov)==prior_seq)
  # get the matching probability vector
  prob_vec = markov[row_num, ]
  # randomly generate the next base depending on the probability vector
  next_base = sample(x=c('a', 'c', 'g', 't'), size=1, prob=prob_vec)
  return(next_base)
}

# simulate sequence
simulate_sequence = function(markov, orig_seq) {
  # Input: Markov probability matrix, original seq (character vector).
  # Output: Simulated sequence (character vector)
  
  print("Simulating sequence...")
  
  # Compute order of Markov Chain
  # nrow(markov) = 4^k
  k = log(dim(markov)[1], base=4)
  
  # Compute length of sequence
  exp_length = length(orig_seq)
  
  # Get initial prior seq from original sequence
  # Initialize simulate sequence
  prior_seq = orig_seq[1:k]
  new_seq = prior_seq
  new_length = length(new_seq)
  
  # keep growing simualted seq until it is as long as original seq
  while (new_length < exp_length){
    # Generate next base based on prior sequence
    next_base = generate_next_base(char2str(prior_seq), markov, k)
    # Append newly generated base to existing simualted sequence
    new_seq =  c(new_seq, next_base)
    # Update length of simulated sequence
    new_length = length(new_seq)
    print(new_length)
    # Update prior sequence
    prior_seq = c(prior_seq[-1], next_base)
    stopifnot(length(prior_seq)==k)
  }
  
  print("Done simulating sequence!")
  
  stopifnot(length(new_seq)==exp_length)
  return(new_seq)
}

# check how close probability matrices based on original and simulated sequences are
check_simulation = function(original_markov, simulated_sequence, threshold=0.001) { 
  # Input: Probability matrix based on original sequence, simulated sequence (char vectors)
  # Output: Proportion of entries in probability matrix based on simulated seq 
  # that are within threshold to entries in original_markov.
  
  print("Checking simulation...")
  
  # Compute order of Markov Chain
  k = log(dim(original_markov)[1], base=4)
  
  # Compute probability matrix based on simulated seq
  simulated_markov = compute_markov(simulated_seq, k)
  
  # Compute proportion of entries in simulated_markov that are close to original_markov
  comparison = abs(simulated_markov - original_markov) <= threshold
  
  print("Done checking simulation!")
  
  return(sum(comparison)/prod(dim(comparison)))
}

# base-wise comparison
base_wise_cf = function(orig_seq, sim_seq){
  # Input: orignal and simulated sequence (char vectors)
  # Output: returns a list of number of exact base-wise matches, proportion of single bases
  
  print("Doing base-wise comparison...")
  
  # number of exact base-wise matches
  match_exact = sum(orig_seq==sim_seq)
  
  # proportion of A, C, G, T in orig_seq and sim_seq
  orig_prop_single = table(orig_seq)
  sim_prop_single = table(sim_seq)
  
  print("Done with base-wise comparison!")
  
  return(list(match_abs = match_exact,
              match_prop = match_exact/length(orig_seq),
              orig_prop_single = orig_prop_single,
              sim_prop_single = sim_prop_single))
}

test_base_wise_cf = function(){
  orig_seq = str2char("cggaccgcgggatgatagatcaaggtaggcataccgtaatatataatattaccgtaccttagggatatcctagaatccgtcccggtagctgaaaagcggt")
  sim_seq = str2char("agagcaagatagactctcgctttgtatggggctaattcactagtcgagggcacctgctgactttgtgtatctcgttggcccatcactaacaagactctgc")
  orig_prop = c(a=30, c=21, g=26, t=23)
  sim_prop = c(a=23, c=25, g=24, t=28)
  match_abs = 19
  match_prop = .19
  test = base_wise_cf(orig_seq, sim_seq)
  return(test$match_abs==match_abs & test$match_prop==match_prop &
           sum(test$orig_prop_single==orig_prop)==length(orig_prop) &
           sum(test$sim_prop_single==sim_prop)==length(sim_prop))
}

# compare transition matrix
markov_cf = function(orig_markov, sim_markov){
  # Input: transition matrices built from original and simulated seqs, order K
  # Output: transition matrix comparison statistic
  
  # compute comparison statistic
  theta = sum(abs(orig_markov-sim_markov)) / prod(dim(orig_markov))
  return(theta)
}

test_markov_cf = function(){
  orig_markov = matrix(c(0,1,0,0, 0.1,0.2,0.6,0.1, 0.21,0.33,0.45, 0.01), nrow=3, ncol=4, byrow=T)
  sim_markov = matrix(c(1,0,0,0, 0.2,0.15,0.34,0.31, 0.5,0.4,0.03,0.07), nrow=3, ncol=4, byrow=T)
  theta = sum(1+1+.1+.05+.26+.21+.29+.07+.42+.06)/(3*4)
  return(identical(theta, markov_cf(orig_markov, sim_markov)))
}

# pairwise global alignment
pairwise_cf = function(orig_seq, sim_seq, m, x, g_open, g_extend){
  # Input: original and simulated seqs (char vectors); 
  # match score m, mismatch score x, gap opening score g_open, gap extension score g_extend
  # Output: pairwise globa alignment score based on pre-defined substitution matrix
  # significant test???
  
  print("Doing pair-wise global alignment...")
  
  # convert to string for use by globalAlign
  orig_seq = char2str(orig_seq)
  sim_seq = char2str(sim_seq)
  
  require(Biostrings)
  # define substitution matrix
  # alternative: use nucleotideSubstitutionMatrix function
  subst.mtx = matrix(data=x, nrow=4, ncol=4,
                     dimnames = rep(list(c('A','C','G','T')), 2))
  diag(subst.mtx) = m
  
  # pairwise global alignment
  globalAligns1s2 = pairwiseAlignment(toupper(orig_seq), toupper(sim_seq),
                                      substitutionMatrix = subst.mtx, 
                                      gapOpening = g_open, gapExtension = g_extend, 
                                      scoreOnly = T)
  
  print("Done with pair-wise global alignment!")
  
  return(globalAligns1s2)
}

# test pairwise_cf
test_pairwise_cf = function(){
  # s1: GAATTC
  # s2: GA-TTA
  # score: 2+2-2-8+2+2-1=-3
  
  s1 = c("g","a","a","t","t","c")
  s2 = c("g","a","t","t","a")
  return(pairwise_cf(s1, s2, m=2, x=-1, g_open=-2, g_extend=-8)==-3)
}


# actual simulation (single run)
simulate_core = function(orig_seq, orig_markov, k, m, x, g_open, g_extend, incl_sim_seq){
  # Input: original seq (char vector), transition matrix based on original seq, order K of MC
  # match score m, mismatch score x, gap opening score g_open, gap extension score g_extend
  # Output: k, sim_seq (string), base-wise cf, markov cf, alignment score
  
  print("Simulating core...")
  
  # generate simulated seq
  sim_seq = simulate_sequence(orig_markov, orig_seq)
  # compute transition matrix based on simulated seq
  sim_markov = compute_markov(sim_seq, k)
  # base-wise comparison
  base_wise_cf_list = base_wise_cf(orig_seq, sim_seq)
  # compare transition matrix
  markov_cf_theta = markov_cf(orig_markov, sim_markov)
  # pairwise global alignment
  pairwise_cf_score = pairwise_cf(orig_seq, sim_seq, m, x, g_open, g_extend)
  # return
  if (incl_sim_seq){
    return(list(k=k, sim_seq = char2str(sim_seq), bw_cf = unlist(base_wise_cf_list),
                mk_cf = markov_cf_theta, pw_cf = pairwise_cf_score))
  } else {
    
    print("Done simulating core!")
    
    return(list(k=k, bw_cf = unlist(base_wise_cf_list),
                mk_cf = markov_cf_theta, pw_cf = pairwise_cf_score))
  }
}

# wrapper
simulate_wrapper = function(orig_seq, k, m, x, g_open, g_extend, num_sim, incl_sim_seq=F){
  # Input: original seq (char vector), order,
  # match score m, mismatch score x, gap opening score g_open, gap extension score g_extend
  # number of simulations, 
  # Output: ?
  
  print("Starting simulation...")
  
  # compute transition matrix based on original sequence
  # before simluate_core so that this is only calculated once
  orig_markov = compute_markov(orig_seq, k)
  
  # actual simulation
  return(replicate(num_sim, simulate_core(orig_seq, orig_markov, k, m, x, g_open, g_extend, incl_sim_seq)))
} 

### for testing
### DO NOT DELETE
# simulate a short-length sequence with a fake markov using MC of fake.k-th order
# fake.k = 7
# fake.mtx = matrix(data=rep(c(0.25, 0.25, 0.37, 0.12), 4^fake.k), nrow=4^fake.k, ncol=4)
# rownames(fake.mtx) = get_combinations(fake.k)
# colnames(fake.mtx) = c('a','c', 'g', 't')
# m = all_seq[[1]]
# test = simulate_sequence(fake.mtx, m[1:1000])


########################################################
###################### COMPARE #########################
########################################################

# compare_simulated = function(original, simulated) {
#   # Input: Sequences to be compared.
#   # Output: Some measure of how similar the sequences are.
#   
# }

########################################################
###################### REFERENCES ######################
########################################################

# http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter1.html
# http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter4.html
# http://cran.r-project.org/web/packages/seqinr/seqinr.pdf
# http://tata-box-blog.blogspot.com/2012/04/introduction-to-markov-chains-and.html

########################################################

# thoughts
# do not delete

# dotPlot(m[1:200], test.split[1:200])
# abline(a=0, b=1, lwd=2, col=2, lty=2)
# # can compare dotPlot with diff k for the same DNA region


# assumption: markov transition matrix stays constant 

# hyptoheses? (same perf btw chr? perf change as k change?)
#compare between chromosome (ANOVA?)
#compare between k (plot)

# visualize & present result analysis

# dotplot to eyeball
# A sliding window plot of GC content (across simulations)

# base-wise comparison, across simulations, avg
# recal prob matrix from simulated seq, some sort of statistics, across simulations, avg?

# pairwise global alignment (N-W algorithm) [BLAST] (avg score [ANOVA?]/ SE across simulations, 
# across k, same substitution matrix, compare)
# statistical significance analysis

# multiple alignment (muscle)

# prolly a long time to run simulation

## TESTING
# test cases for compute_markov
# test cases for everything


### combine data
setwd("~/University/2013 Fall/36350/project/sim")

# only works if only simulation .Rdata files (nothing else) are in ~/sim
to_load = list.files()
for (i in 1:length(to_load)){
  load(to_load[i])
}

# create a data matrix/frame containing statistics to be analyzed

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

### analyze
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