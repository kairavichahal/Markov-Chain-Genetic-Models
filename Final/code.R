################################################################################
# 36-350 Final Project                   #######################################
# Markov Chain Genetic Models            #######################################
# Kairavi Chahal, Tony Yang, Julian Zhou #######################################
################################################################################

################################################################################
# LOAD REQUIRED PACKAGES #######################################################
################################################################################

# install.packages('seqinr')
# source("http://bioconductor.org/biocLite.R")
# biocLite('Biostrings')
library(seqinr)
library(Biostrings)

################################################################################
# READ IN DATA #################################################################
################################################################################

# Read in FASTA file
raw_data = read.fasta(file="dicty_chromosomal")
all_seq = getSequence(raw_data)
rm(raw_data)

# Take care of missing values. Replace the Ns in sequence with a, c, g, t 
# uniformly with prob. 0.25.
replace_n = function(chr_vec) {
  #Input: a character vector with possible unknown 'n' patterns.
  #Output: a character vector with 'n' replaced by a, c, g, t.
  idx = which(chr_vec == "n")
  chr_vec[idx] = sample(c('a', 'c', 'g', 't'), length(idx), replace=T)
  return(chr_vec)
}

# Save each chromosome in a separate variable.
chr.m  = replace_n(all_seq[[1]])
chr.2f = replace_n(all_seq[[2]])
chr.3f = replace_n(all_seq[[3]])
chr.bf = replace_n(all_seq[[4]])
chr.1  = replace_n(all_seq[[5]])
chr.2  = replace_n(all_seq[[6]])
chr.3  = replace_n(all_seq[[7]])
chr.4  = replace_n(all_seq[[8]])
chr.5  = replace_n(all_seq[[9]])
chr.6  = replace_n(all_seq[[10]])
chr.r  = replace_n(all_seq[[11]])

rm(all_seq)

################################################################################
# COMPUTE PROBABILITY MATRIX ###################################################
################################################################################

# Helper function to convert transition matrix to probabilities.
prob_of = function(mtx_row) {
  # Input: a row of the transition matrix with counts of each prior sequence.
  # Output: same row of matrix, with counts converted to probabilities.
  if (sum(mtx_row)>0) {
    mtx_row = mtx_row/sum(mtx_row)
  } else {
    # Replace with 0.25 if this sub-sequence does not exist in the original.
    mtx_row = rep(1/4, 4)
  }
  return(mtx_row)
}

# Helper function to create all combinations of ACGT in kth order transition 
# matrix.
get_combinations = function(k) {
  # Input: integer; order of the Markov chain to be computed.
  # Output: list of all possible prior sequences of length k.
  l = rep(list(c('a','c','g','t')), k)
  grid = as.matrix(expand.grid(l))
  rownames = sort(apply(grid, 1, paste, collapse=""))
  stopifnot(length(rownames)==4^k)
  return(rownames)
}

# Helper function to get subsequences of original chromosome.
get_substring = function(i, k, character_vector) {
  # Input: index, order of Markov and character vector
  # Output: substring containing (i-k)th to ith characters of the sequence
  substring = char2str(character_vector[(i-k):i])
  return(substring)
}

# Main function. Call this on a chromosome and specify order k.
compute_markov = function(sequence, k=1) {
  # Input: Sequence as a character vector and what order Markov we want.
  # Output: Markov probability matrix of dimension 4^k by 4.
  
  print("Computing Markov...")
  
  # Get subsequences of length k.
  sub_seqs = sapply((k+1):length(sequence), get_substring, k=k, 
                    character_vector=sequence)
  # Get counts of how many times each subsequences appears.
  sub_seqs_tab = table(sub_seqs)
  
  markov = matrix(0, ncol=4, nrow=4^k, dimnames=list(get_combinations(k), 
                                                     c('a', 'c', 'c', 't')))
  
  # Extract prior sequence and current base from rownames.
  sub_seqs_rows = substr(rownames(sub_seqs_tab), 1, k)
  sub_seqs_cols = substr(rownames(sub_seqs_tab), k+1, k+1)
  
  # Create matrix with columns: Prior sequences, current base, number of times
  # that combination appeared in the chromosome.
  mat = cbind(sub_seqs_rows, sub_seqs_cols, sub_seqs_tab)
  
  # Find each element of 'markov' in 'mat'.
  row_idx = match(mat[,1], rownames(markov))
  col_idx = match(mat[,2], colnames(markov))

  markov[cbind(row_idx, col_idx)] = as.numeric(mat[, 3])
  
  print("Done computing Markov!")
  return(t(apply(markov, 1, prob_of)))
}

################################################################################
# SIMULATE SEQUENCES ###########################################################
################################################################################

# Helper function that generates next base given prior.
generate_next_base = function(prior_seq, markov, k) {
  # Input: prior sequence (string), probability matrix, order k.
  # Output: a single character, which will be the next base.
  
  # Find the matching row number in 'markov' given a prior sequence.
  row_num = which(rownames(markov) == prior_seq)
  
  # Get the matching probability vector.
  prob_vec = markov[row_num,]
  
  # Randomly generate the next base depending on the probability vector.
  next_base = sample(x=c('a', 'c', 'g', 't'), size=1, prob=prob_vec)
  return(next_base)
}

# Main function. Call this on Markov matrix and original chromosome sequence.
simulate_sequence = function(markov, orig_seq) {
  # Input: Markov probability matrix, original sequence (character vector).
  # Output: Simulated sequence (character vector).
  
  print("Simulating sequence...")
  
  # Compute order of Markov chain.
  k = log(dim(markov)[1], base=4)

  exp_length = length(orig_seq)
  
  # Get first prior sequence from original sequence.
  prior_seq = paste(orig_seq[1:k], collapse='')
  new_seq = prior_seq
  new_length = nchar(new_seq)
  
  # Keep growing simulated sequence until it is as long as original sequence.
  while (new_length < exp_length) {
    # Generate next base based on prior sequence.
    next_base = generate_next_base(prior_seq, markov, k)
    
    # Append newly generated base to existing simulated sequence.
    new_seq =  paste(new_seq, next_base, sep='')
    
    # Update length of simulated sequence.
    new_length = nchar(new_seq)
    
    # Update prior sequence.
    prior_seq = paste(substring(prior_seq, 2), next_base, sep='')
    stopifnot(nchar(prior_seq) == k)
  }
  
  print("Done simulating sequence!")
  
  stopifnot(nchar(new_seq) == exp_length)
  return(str2char(new_seq))
}

################################################################################
# COMPARE SEQUENCES ############################################################
################################################################################

# Base-wise comparison.
base_wise_cf = function(orig_seq, sim_seq) {
  # Input: original and simulated sequence (character vectors).
  # Output: returns a list of number of exact base-wise matches and proportions 
  #         of each base.
  
  print("Doing base-wise comparison...")
  
  # Number of exact base-wise matches.
  match_exact = sum(orig_seq == sim_seq)
  
  # Proportion of a, c, g, t in orig_seq and sim_seq.
  orig_prop_single = table(orig_seq)
  sim_prop_single = table(sim_seq)
  
  print("Done with base-wise comparison!")
  return(list(match_abs = match_exact, match_prop = match_exact/length(orig_seq), 
              orig_prop_single = orig_prop_single, sim_prop_single = sim_prop_single))
}

# Compare transition matrices.
markov_cf = function(orig_markov, sim_markov) {
  # Input: transition matrices built from original and simulated sequences.
  # Output: transition matrix comparison statistic.
  
  theta = sum(abs(orig_markov - sim_markov))/prod(dim(orig_markov))
  return(theta)
}

# Pairwise Global Alignment
pairwise_cf = function(orig_seq, sim_seq, m, x, g_open, g_extend) {
  # Input: original and simulated sequences. Match score m, mismatch score x,
  #        gap opening score g_open, gap extension score g_extend.
  # Output: Pairwise Global Alignment score based on predefined substitution matrix.
  
  print("Doing pair-wise global alignment...")
  
  # Convert to string for use by globalAlign.
  orig_seq = char2str(orig_seq)
  sim_seq = char2str(sim_seq)
  
  require(Biostrings)
  # Define substitution matrix.
  # Alternative: use nucleotideSubstitutionMatrix function.
  subst.mtx = matrix(data=x, nrow=4, ncol=4, dimnames = rep(list(c('A','C','G','T')), 2))
  diag(subst.mtx) = m
  
  # pairwise global alignment
  globalAligns1s2 = pairwiseAlignment(toupper(orig_seq), toupper(sim_seq), 
                                      substitutionMatrix = subst.mtx, gapOpening = g_open, 
                                      gapExtension = g_extend, scoreOnly = T)
  
  print("Done with pair-wise global alignment!")
  return(globalAligns1s2)
}

################################################################################
# RUN SIMULATIONS ##############################################################
################################################################################

# Single run of simulation.
simulate_core = function(orig_seq, orig_markov, k, m, x, g_open, g_extend, incl_sim_seq) {
  # Input: original sequence, transition matrix based on original sequence,
  #        order k, match score m, mismatch score x, gap opening score g_open, 
  #        gap extension score g_extend.
  # Output: list of k, sim_seq (string), base-wise_cf output, markov_cf output, 
  #         Pairwise Global Alignment score.
  
  print("Running simulation...")
  
  # Generate simulated sequence.
  sim_seq = simulate_sequence(orig_markov, orig_seq)
  # Compute transition matrix based on simulated sequence.
  sim_markov = compute_markov(sim_seq, k)
  # Do base-wise comparison.
  base_wise_cf_list = base_wise_cf(orig_seq, sim_seq)
  # Compare transition matrices.
  markov_cf_theta = markov_cf(orig_markov, sim_markov)
  # Do Pairwise Global Alignment.
  pairwise_cf_score = pairwise_cf(orig_seq, sim_seq, m, x, g_open, g_extend)
  
  if (incl_sim_seq) {
    print("Done with simulation!")
    return(list(k=k, sim_seq = char2str(sim_seq), bw_cf = unlist(base_wise_cf_list), 
                mk_cf = markov_cf_theta, pw_cf = pairwise_cf_score))
  } else {
    print("Done with simulation!")
    return(list(k=k, bw_cf = unlist(base_wise_cf_list), mk_cf = markov_cf_theta, 
                pw_cf = pairwise_cf_score))
  }
}

# Final function that runs everything.
simulate_wrapper = function(orig_seq, k, m, x, g_open, g_extend, num_sim, incl_sim_seq=F) {
  # Input: original seq , order, match score m, mismatch score x, gap opening score g_open,
  #        gap extension score g_extend, number of simulations.
  # Output: 
  
  print("Starting simulation...")
  
  # Compute transition matrix based on original sequence before simluate_core so
  # that this is only calculated once.
  orig_markov = compute_markov(orig_seq, k)
  return(replicate(num_sim, simulate_core(orig_seq, orig_markov, k, m, x, g_open, 
                                          g_extend, incl_sim_seq)))
}

################################################################################
# TEST FUNCTIONS################################################################
################################################################################

test_base_wise_cf = function() {
  orig_seq = str2char("cggaccgcgggatgatagatcaaggtaggcataccgtaatatataatattaccgtaccttagggatatcctagaatccgtcccggtagctgaaaagcggt")
  sim_seq = str2char("agagcaagatagactctcgctttgtatggggctaattcactagtcgagggcacctgctgactttgtgtatctcgttggcccatcactaacaagactctgc")
  orig_prop = c(a=30, c=21, g=26, t=23)
  sim_prop = c(a=23, c=25, g=24, t=28)
  match_abs = 19
  match_prop = .19
  test = base_wise_cf(orig_seq, sim_seq)
  return(test$match_abs==match_abs & test$match_prop==match_prop & sum(test$orig_prop_single==orig_prop)==length(orig_prop) 
         & sum(test$sim_prop_single==sim_prop)==length(sim_prop))
}

test_markov_cf = function() {
  orig_markov = matrix(c(0,1,0,0, 0.1,0.2,0.6,0.1, 0.21,0.33,0.45, 0.01), nrow=3, ncol=4, byrow=T)
  sim_markov = matrix(c(1,0,0,0, 0.2,0.15,0.34,0.31, 0.5,0.4,0.03,0.07), nrow=3, ncol=4, byrow=T)
  theta = sum(1+1+.1+.05+.26+.21+.29+.07+.42+.06)/(3*4)
  return(identical(theta, markov_cf(orig_markov, sim_markov)))
}

test_pairwise_cf = function() {
  # s1: GAATTC
  # s2: GA-TTA
  # score: 2+2-2-8+2+2-1=-3
  
  s1 = c("g","a","a","t","t","c")
  s2 = c("g","a","t","t","a")
  return(pairwise_cf(s1, s2, m=2, x=-1, g_open=-2, g_extend=-8)==-3)
}

################################################################################
# MISCELLANEOUS FUNCTIONS ######################################################
################################################################################

# Character vector to string.
char2str = function(char_vec) {
  return(paste(char_vec, collapse=''))
}

# String to character vector.
str2char = function(string) {
  return(unlist(strsplit(string, '')))
}

# Function below simply return number of seconds taken to execute function and
# output of given function.
time_function = function(func, ...) {
  start_time = proc.time()
  output = func(...)
  stop_time = proc.time() - start_time
  return(list(stop_time[3], function_output=output))
}

################################################################################
# REFERENCES ###################################################################
################################################################################

# http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter1.html
# http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter4.html
# http://cran.r-project.org/web/packages/seqinr/seqinr.pdf
# http://tata-box-blog.blogspot.com/2012/04/introduction-to-markov-chains-and.html