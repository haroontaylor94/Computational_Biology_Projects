#Hopfield Networks
#Plotting function
plot_pattern <- function(x){
    image(t(matrix(x, nrow = 5, byrow = TRUE)), col = c("white", "black"), 
          xaxt = "n", yaxt = "n", ylim = c(1.1, -0.1))
}

#Function to calculate weight matrix
calc_w <- function(patterns){
    w <- matrix(0, ncol(patterns), ncol(patterns))
    for(i in 1:nrow(patterns)){
        wn <- outer(patterns[i, ], patterns[i, ])
        w <- w + wn
    }
    diag(w) <- 0
    return(w)
}

#Sigma function
theta <- function(x){
    if (x<0) {x = -1} else {x = 1}
    return(x)
}

#Asynchronous update function
async <- function(w, pattern){
    x <- sample.int(nrow(w), nrow(w), replace = FALSE)
    for (i in x){
        wi = w[i,]
        a <- sum(wi*pattern)
        pattern[i] <- theta(a) 
    }
    return(pattern)
}

#Patterns to store
H <- c(1,-1,-1,-1,1, 1,-1,-1,-1,1, 1,1,1,1,1, 1,-1,-1,-1,1, 1,-1,-1,-1,1)
M <- c(1,-1,-1,-1,1, 1,1,-1,1,1, 1,-1,1,-1,1, 1,-1,-1,-1,1, 1,-1,-1,-1,1)
T <- c(1,1,1,1,1, -1,-1,1,-1,-1, -1,-1,1,-1,-1, -1,-1,1,-1,-1, -1,-1,1,-1,-1)
one <- c(-1,-1,1,1,-1, -1,-1,-1,1,-1, -1,-1,-1,1,-1, -1,-1,-1,1,-1, -1,-1,1,1,1)

par(mfrow = c(2,2))
plot_pattern((H))
plot_pattern((M))
plot_pattern((T))
plot_pattern((one))

#Function to find storage capacity for a given spareness
find_capacity <- function(sp, new = FALSE){
    sp_vec <- c(rep(1, sp), rep(-1, I-sp))
    for(i in 2:I){
        patterns <- t(replicate(i, sample(sp_vec)))
        while(anyDuplicated.matrix(patterns) > 0){
            patterns <- t(replicate(i, sample(sp_vec)))
        }
        if(new){
            w <- calc_w_new(patterns)
        } else {
            w <- calc_w(patterns)
        }
        for(j in 1:nrow(patterns)){
            pattern <- patterns[j, ]
            new_pattern <- async(w, pattern)
            if(!isTRUE(all.equal(new_pattern, pattern))){
                return(i-1)
            }
        }
    }
}

#Function to calculate storage capacity for each successive spareness
sc_vs_sp <- function(new = FALSE){
    capacities <- rep(0, I-1)
    for(sp in 1:(I-1)){
        sp_cap <- rep(0, reps)
        for(i in 1:reps){
            capacity <- find_capacity(sp, new = new)
            sp_cap[i] <- capacity
        }
        av_cap <- mean(sp_cap)
        capacities[sp] <- av_cap
        print(sp)
    }
    return(capacities)
}

#Set number of neurons
I <- 200
reps <- 10
#Assess the effect of spareness on storage capcity
capacities <- sc_vs_sp()

png("Storage_capacity.png", width = 4, height = 4, units = "in", res = 300)
par(mfrow = c(1,1))
plot(1:(I-1)/I, capacities/I, type = "l", xlab = "Sparseness", 
     ylab = "Average storage capacity (N/I)")
dev.off()

#Assess how pattern recall is affected by random loss of weights
I <- 200
no_update <- 10
weights_lost <- seq(0, (I^2-I)/2, 50)
#Function to calculate recall for each successive weight lost
recall_vs_lw <- function(patterns, w){
    w_rows <- seq(I)
    lower_t_w_ind <- 
        cbind(row = unlist(lapply(2:I, function(x) x:I), use.names = FALSE),
              col = rep(w_rows[-length(w_rows)], times = rev(tail(w_rows, -1))-1))
    av_recalls <- rep(NA, length(weights_lost))
    for(i in 1:length(weights_lost)){
        lw <- weights_lost[i]
        lw_ind <- lower_t_w_ind[sample(nrow(lower_t_w_ind), lw), ]
        w2 <- w
        w2[lw_ind] <- 0
        w2[lw_ind[, c(2,1)]] <- 0
        recalls <- rep(NA, no_patterns)
        for(j in 1:no_patterns){
            pattern <- patterns[j, ]
            flip_bits <- sample.int(I, 1)
            pattern[flip_bits] <- pattern[flip_bits] * -1
            for(k in 1:no_update){
                pattern <- async(w2, pattern)
            }
            recall <- as.numeric(pattern %*% patterns[j, ]/I)
            recalls[j] <- recall 
        }
        av_recall <- mean(recalls)
        av_recalls[i] <- av_recall
        print(lw)
    }
    return(av_recalls)
}

#Create 10 unique patterns to store
patterns <- t(replicate(10, sample(c(1,-1), I, replace = TRUE)))
while(anyDuplicated.matrix(patterns) > 0){
    patterns <- t(replicate(10, sample(c(1,-1), I, replace = TRUE)))
}
no_patterns <- nrow(patterns)
w <- calc_w(patterns)

av_recalls <- recall_vs_lw(patterns, w)

png("Weights_lost.png", width = 4, height = 4, units = "in", res = 300)
par(mfrow = c(1,1))
plot((2*weights_lost)/(I^2-I), av_recalls, xlab = "Fraction of weights lost", 
     ylab = "Average Recall", type = "l", ylim = c(0,1))
dev.off()

#Alternative method for setting the weights
#Sigmoid function
sigmoid <- function(x){
    1/(1+exp(-x))
}

calc_w_new <- function(x, L=10, eta=0.1, alpha=0){
    w <- t(x) %*% x
    for(i in 1:L){
        diag(w) <- 0
        a <- x %*% w
        y <- sigmoid(a)
        t <- x
        t[which(t == -1)] <- 0
        e <- t - y
        gw <- t(x) %*% e
        gw <- gw + t(gw)
        w <- w + eta * (gw - alpha * w)
    }
    diag(w) <- 0
    return(w)
}

#Effect on storage capacity
png("Storage_capacity_new.png", width = 4, height = 4, units = "in", res = 300)
plot(1:(I-1)/I, capacities/I, type = "l", xlab = "Sparseness", 
     ylab = "Average storage capacity (N/I)", ylim = c(0,0.4))
capacities_new <- sc_vs_sp(new = TRUE)
lines(1:(I-1)/I, capacities_new/I, col = "red")
dev.off()

#Effect on robustness
png("Weights_lost_new.png", width = 4, height = 4, units = "in", res = 300)
plot((2*weights_lost)/(I^2-I), av_recalls, xlab = "Fraction of weights lost", 
     ylab = "Average Recall", type = "l", ylim = c(0,1))
w_new <- calc_w_new(patterns)
av_recalls_new <- recall_vs_lw(patterns, w_new)
lines((2*weights_lost)/(I^2-I), av_recalls_new, col = "red")
dev.off()





