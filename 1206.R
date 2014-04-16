# 36-350 Final Project
# Markov Chain Genetic Models
# Kairavi Chahal, Tony Yang, Julian Zhou

########################################################
############### LOAD REQUIRED PACKAGES #################
########################################################
# install.packages('seqinr')
# source("http://bioconductor.org/biocLite.R")
# library(BiocInstaller)
# biocLite('Biostrings')
library(seqinr)
library(Biostrings)

########################################################
#################### READ SEQUENCE #####################
########################################################

# Read in FASTA file
raw_data = read.fasta(file = "dicty_chromosomal")
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

# If you want separate variables in each chromosome
# chr.m  = replace_n(all_seq[[1]])
chr.2f = replace_n(all_seq[[2]])
# chr.3f = replace_n(all_seq[[3]])
# chr.bf = replace_n(all_seq[[4]])
# chr.1  = replace_n(all_seq[[5]])
# chr.2  = replace_n(all_seq[[6]])
# chr.3  = replace_n(all_seq[[7]])
chr.4  = replace_n(all_seq[[8]])
chr.5  = replace_n(all_seq[[9]])
chr.6  = replace_n(all_seq[[10]])
# chr.r  = replace_n(all_seq[[11]])

rm(all_seq)
#all_seq_wo_n = lapply(all_seq, replace_n)

# character vector to string
char2str = function(char_vec){
  return(paste(char_vec, collapse=''))
}

# string to character vector
str2char = function(string){
  return(unlist(strsplit(string, '')))
}

########################################################
############## COMPUTE PROBABILITY MATRIX ##############
########################################################

# Divide each row in markov by its sum
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
  prior_seq = paste(orig_seq[1:k], collapse='')
  new_seq = prior_seq
  new_length = nchar(new_seq)
  
  # keep growing simualted seq until it is as long as original seq
  while (new_length < exp_length){
    # Generate next base based on prior sequence
    next_base = generate_next_base(prior_seq, markov, k)
    # Append newly generated base to existing simualted sequence
    new_seq =  paste(new_seq, next_base, sep='')
    # Update length of simulated sequence
    new_length = nchar(new_seq)
    #print(new_length)
    # Update prior sequence
    prior_seq = paste(substring(prior_seq, 2), next_base, sep='')
    stopifnot(nchar(prior_seq)==k)
  }
  
  print("Done simulating sequence!")
  
  stopifnot(nchar(new_seq)==exp_length)
  return(str2char(new_seq))
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

# 
# 
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


########################################################
###################### TIME FUNCTION ###################
########################################################
# Function below simply return number of seconds taken
# to execute func with arguments ...
time_function = function(func, ...) {
  start_time = proc.time()
  func(...)
  stop_time = proc.time() - start_time
  return(stop_time[3])
}