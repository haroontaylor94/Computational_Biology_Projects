#Function to create markov chain from transition matrix and initial distribution
markov_chain <- function(A, mu0, N){
    A <- read.csv(A, check.names = FALSE)
    mu0 <- read.csv(mu0, check.names = FALSE)
    S <- as.numeric(colnames(A))
    mchain <- rep(NA, N)
    mchain[1] <- sample(S, 1, prob = mu0)
    for (i in 2:N){
        mchain[i] <- sample(S, 1, prob = A[as.character(mchain[i - 1]), ])
    }
    return(mchain)
}

#Function to infer maximum likelihood transition matrix and initial distribution
#from markov chain
ml_mu0_trans <- function(mchain, S){
    mchain <- as.numeric(scan(mchain, sep = ",", quiet = TRUE))
    N <- length(mchain)
    mu0 <- as.numeric(S == mchain[1])
    trans <- matrix(NA, nrow = length(S), ncol = length(S))
    dimnames(trans) <- list(S, S)
    for(i in S){
        ni <- sum(mchain[-N] == i)
        for (j in S){
            nij <- sum(mchain[which(mchain[-N] == i) + 1] == j)
            trans[as.character(i), as.character(j)] <- nij/ni
        }
    }
    return(list(mu0 = mu0, trans = trans))
}

#Function to implement hmm by adding emitted variable at each point in the chain
hmm <- function(A, mu0, B, N){
    hidden <- markov_chain(A, mu0, N)
    B <- read.csv(B, check.names = FALSE)
    V <- as.numeric(colnames(B))
    observed <- sapply(hidden, 
                       function(x) sample(V, 1, prob = B[as.character(x), ]))
    return(list(hid = hidden, obs = observed))
}

N <- 115
disc_hmm <- hmm("transition_matrix.csv", "initial_distribution.csv",
                "emission_matrix.csv", N)
#Plot emission and hidden states
par(xpd = TRUE)
plot(1:N, disc_hmm$obs, type = "l", col = "red", ylim = c(0, 5), 
     xlab = "Sequence position", ylab = "State")
lines(1:N, disc_hmm$hid, col = "blue")
legend(0, 6, title = "States", legend = c("hidden", "emitted"), 
       col = c("blue", "red"), lty = 1, cex = 0.75)

#Function to implement the scaled forward algorithm
forward <- function(seq, trans, emit, mu0){
    seq <- scan(seq, sep = ",", quiet = TRUE)
    S <- as.numeric(colnames(trans))
    V <- as.numeric(colnames(emit))
    c <- rep(NA, length(seq))
    alpha <- matrix(NA, ncol = length(S), nrow = length(seq))
    alpha[1, ] <- mu0 * emit[, as.character(seq[1])]
    c[1] <- sum(alpha[1, ])
    alpha[1, ] <- alpha[1, ]/c[1]
    for(i in 2:length(seq)){
        alpha[i, ] <- emit[,as.character(seq[i])] * 
            sapply(S, function(x) sum(alpha[i-1, ] * trans[,as.character(x)]))
        c[i] <- sum(alpha[i, ])
        alpha[i, ] <- alpha[i, ]/c[i]
    }
    return(list(alphas = alpha, scale_factors = c))
}

#Function to calculate the likelihood of an emitted sequence given the model
l_emit <- function(seq, trans, emit, mu0){
    forward_variables <- forward(seq, trans, emit, mu0)
    c <- forward_variables$scale_factors
    sum(log(c))
} 

#Calculate GC content of S.cerevisiae in 100bp windows
library(seqinr)
S.cerevisiae <- read.fasta("S.cerevisiae_chrom_III.fa")
seq <- S.cerevisiae$III
seq_index <- seq(1, length(seq), 100)
GC_vec <- rep(0, length(seq_index))
count <- 1
for (i in seq_index){
    if (i + 100 <= length(seq)){
        sub_seq <- seq[i:(99+i)]
    } else {
        sub_seq <- seq[i:length(seq)]
    }
    GC <- sum(sub_seq == "c" | sub_seq == "g")/length(sub_seq)
    GC_vec[count] <- GC
    count <- count + 1
}

GC_bin <- cut(GC_vec, 5, labels = 1:5)

#Function to make barplot of binning scheme
make_GC_barplot <- function(GC_vec, bins){
    GC_bin <- cut(GC_vec, bins)
    GC_barplot <- table(GC_bin)
    names(GC_barplot) <- gsub(",", "-", substring(names(GC_barplot), 2, 12))
    barplot(GC_barplot, xlab = "GC content", ylab = "Frequency", col = "blue") 
}

#Plot binning scheme barplot
make_GC_barplot(GC_vec, 5)
write.table(GC_bin, "GC_sequence.csv", row.names = FALSE, col.names = FALSE, 
            quote = FALSE)

#Set up model parameters
S <- c(0, 1)
V <- c(1, 2, 3, 4, 5)
A <- matrix(c(0.8, 0.1, 0.2, 0.9), ncol = 2)
B <- matrix(c(0.2, 0, 0.5, 0.1, 0.2, 0.4, 0.1, 0.4, 0, 0.1), ncol = 5)
mu0 <- c(0.5, 0.5)

dimnames(A) <- list(S, S)
dimnames(B) <- list(S, V)

#Calculate log likelihood
GC_l <- l_emit("GC_sequence.csv", A, B, mu0)
cat(paste("log likelihood =", round(GC_l, 3)))

backward <- function(seq, trans, emit, c){
    seq <- scan(seq, sep = ",", quiet = TRUE)
    S <- as.numeric(colnames(trans))
    V <- as.numeric(colnames(emit))
    n <- length(seq)
    beta <- matrix(0, ncol = length(S), nrow = n)
    beta[n, ] <- 1/c[n]
    for (i in (n-1):1){
        beta[i, ] <- sapply(S, 
                            function(x) sum(beta[i+1, ] * 
                                                trans[as.character(x),] * 
                                                emit[,as.character(seq[i+1])]))
        beta[i, ] <- beta[i, ]/c[i]
    }
    return(beta)
}

#Function to calculate gamma for specific timepoint
gamma <- function(alpha, beta){
    (alpha * beta)/sum(alpha * beta)
}

#Function to calculate xi for a specific timepoint
xi <- function(alpha, beta, trans, emit){
    t(t(alpha * trans) * (beta * emit))/sum(t(t(alpha * trans) * (beta * emit)))
}

#Function to perform Baum-Welch algorithm
baum_welch <- function(seq, trans, emit, mu0){
    dimnames(emit) <- list(S, V)
    forward_variables <- forward(seq, trans, emit, mu0)
    alpha <- forward_variables$alphas
    c <- forward_variables$scale_factors
    beta <- backward(seq, trans, emit, c)
    seq <- scan(seq, sep = ",", quiet = TRUE)
    mu0 <- gamma(alpha[1, ], beta[1, ])
    N <- length(seq)
    gammas <- t(sapply(1:N, function(x) gamma(alpha[x, ], beta[x, ])))
    xis <- lapply(1:(N-1), 
                  function(x) xi(alpha[x, ], beta[x+1, ], 
                                 trans, emit[,as.character(seq[x+1])]))
    trans <- Reduce("+", xis) / apply(gammas[1:(N-1), ], 2, sum)
    emit <- sapply(V, function(x) 
        apply(gammas * as.numeric(seq == x), 2, sum)/apply(gammas, 2, sum))
    return(list(mu0 = mu0, A = trans, B = emit))
}

bm_estimation <- function(seq, trans, emit, mu0, imax, ep){
    log_l <- l_emit(seq, trans, emit, mu0)
    for (i in 1:imax){
        params <- baum_welch(seq, trans, emit, mu0)
        mu0 <- params$mu0
        trans <- params$A
        emit <- params$B
        dimnames(trans) <- list(S, S)
        dimnames(emit) <- list(S, V)
        log_l_diff <- abs(log_l - l_emit(seq, trans, emit, mu0))
        if (log_l_diff < ep){
            break
        } else {
            log_l <- l_emit(seq, trans, emit, mu0)
        }
    }
    return(list(mu0 = mu0, A = trans, B = emit, i = i, log_l = log_l))
}

GC_params <- bm_estimation("GC_sequence.csv", A, B, mu0, 100, 10^-8)

cat(paste("log likelihood =", round(GC_params$log_l, 3)))

#Function to perform non-logarithmic viterbi
viterbi <- function(seq, trans, emit, u0){
    seq <- scan(seq, sep = ",", quiet = TRUE) 
    S <- as.numeric(colnames(trans))
    V <- as.numeric(colnames(emit))
    N <- length(seq)
    deltas <- matrix(NA, ncol = length(S), nrow = N)
    psis <- matrix(NA, ncol = length(S), nrow = N)
    hid <- rep(NA, N)
    deltas[1, ] <- u0 * emit[,as.character(seq[1])]
    deltas[1, ] <- deltas[1, ]/sum(deltas[1, ])
    for (i in 2:N){
        deltas[i, ] <- emit[,as.character(seq[i])] * 
            apply(deltas[i-1, ] * trans, 2, max)
        deltas[i, ] <- deltas[i, ]/sum(deltas[i, ])
        psis[i, ] <- apply(deltas[i-1, ] * trans, 2, 
                           function(x) which(x == max(x)))
    }
    hid[i] <- which(deltas[i, ] == max(deltas[i, ]))
    for (i in (N-1):1){
        hid[i] <- psis[i+1, hid[i+1]]
    }
    hid <- sapply(hid, function(x) S[x])
}

GC_hid <- viterbi("GC_sequence.csv", GC_params$A, GC_params$B, GC_params$mu0)
n <- 200
#Plot GC content and hidden states
par(xpd = TRUE)
plot(1:n, GC_vec[1:n], type = "l", col = "red", ylim = c(0.1, 0.7),
     xlab = "Genomic position (100 bp)", ylab = "GC content")
par(new = TRUE)
plot(1:n, GC_hid[1:n], type = "l", axes = FALSE, xlab = NA, ylab = NA, 
     col = "blue")
axis(side = 4, at=c(0,1), labels = c(0,1))
mtext(side = 4, line = 1, "Hidden State")
legend(0, 1.2, legend = c("GC content", "hidden states"), 
       col = c("red", "blue"), lty = 1, cex = 0.75)