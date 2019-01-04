#Theta function
theta <- function(k, imax, max = 0.9, min = 0.4){
    max - ((max - min)/imax) * k
}

#Function to update particle speeds
update_speeds <- function(speeds, k, pos, best_pos, global_best, imax, N){
    global <- apply(pos, 1, function(x) global_best - x)
    if (is.matrix(global)){
        global <- t(global)
    } else {
        global <- matrix(global)
    }
    theta(k, imax) * speeds + 2 * runif(N) * (best_pos - pos) + 
        2 * runif(N) * global
}

#Particle swarm optimisation function per dimension
PSO <- function(N, search_int, precision, imax, f, max = TRUE, verbose = TRUE){
    pos <- apply(search_int, 1, function(x) runif(N, x[1], x[2]))
    speeds <- matrix(0, nrow = N, ncol = nrow(search_int))
    best_pos <- pos
    for (k in 0:imax){
        if (all(apply(pos, 2, function(x) diff(range(x))) < precision)){
            if(verbose){
                cat(paste("PSO converged in", k, "iterations"))
            }
            return(global_best)
        } else {
            if (max){
                pos_index <- which(apply(pos, 1, f) > apply(best_pos, 1, f))
                best_pos[pos_index, ] <- pos[pos_index, ]
                global_best <- best_pos[which.max(apply(best_pos, 1, f)), ]
            } else {
                pos_index <- which(apply(pos, 1, f) < apply(best_pos, 1, f))
                best_pos[pos_index, ] <- pos[pos_index, ]
                global_best <- best_pos[which.min(apply(best_pos, 1, f)), ]
            }
            
            speeds <- update_speeds(speeds, k, pos, best_pos, 
                                    global_best, imax, N)
            pos <- pos + speeds
        }
    }
    return(global_best)
}

#Function 1 optimisation
f_1 <- function(x){
    -x^2 + 2*x + 11
}
f_1_PSO <- PSO(100, matrix(c(-2, 2), ncol = 2), 10^-4, 1000, f_1)

f_1_optim <- optimise(f_1, c(-2, 2), maximum = TRUE)$maximum
cat(paste("PSO maximum point =", signif(f_1_PSO, 4)))
cat(paste("Optim maximum point =", f_1_optim))

#Function 2 optimisation
f_2 <- function(x){
    if (x[1] >= 0 && x[2] >= 0 && x[1] + 4*x[2] <= 5 && 2*x[1] + 3*x[2] <= 6){
        result <- x[1]^2 + x[2]^2 - 2*x[1] - 4*x[2] 
    } else {
        result <- Inf
    }
    return(result)
}

f_2_PSO <- PSO(100, matrix(c(0, 3, 0, 1.25), ncol = 2, byrow = TRUE), 
               10^-3, 1000, f_2, max = FALSE)

f_2_optim <- optim(c(0,0), f_2)$par
cat(paste0("PSO minimum point = (", round(f_2_PSO[1], 3), " , ",
           round(f_2_optim[2], 3), ")"))

cat(paste0("Optim minimum point = (", round(f_2_optim[1], 3), " , ", 
           round(f_2_optim[2], 3), ")"))

f_3 <- function(x){
    sum(x^2)
}
f_3_accuracy <- matrix(NA, nrow = 3, ncol = 9)
for (i in 2:10){
    PSO_100 <- PSO(100, matrix(rep(c(-5.12, 5.12), i), ncol = 2, byrow = TRUE),
                   10^-3, 1000, f_3, max = FALSE, verbose = FALSE)
    acc_100 <- mean(abs(PSO_100))
    PSO_1000 <- PSO(1000, matrix(rep(c(-5.12, 5.12), i), ncol = 2, byrow = TRUE),
                    10^-3, 1000, f_3, max = FALSE, verbose = FALSE)
    acc_1000 <- mean(abs(PSO_1000))    
    f_3_optim <- optim(runif(i, -5.12, 5.12), f_3)$par
    acc_optim <- mean(abs(f_3_optim))
    f_3_accuracy[,(i-1)] <- c(acc_optim, acc_100, acc_1000)
}

f_3_accuracy <- -log10(f_3_accuracy)

par(xpd = TRUE)
plot(2:10, f_3_accuracy[1,], pch=0, cex = 0.75, col = "red", ylim = c(1,25), 
     xlab = "Number of dimensions", 
     ylab = "Accuracy (-log10 mean absolute difference from 0)")
lines(2:10, f_3_accuracy[1,], col = "red")
points(2:10, pch = 0, cex = 0.75, f_3_accuracy[2,], col = "blue")
lines(2:10, f_3_accuracy[2,], col = "blue")
points(2:10, pch = 0, cex = 0.75, f_3_accuracy[3,], col = "green3")
lines(2:10, f_3_accuracy[3,], col = "green3")
legend(2, 30, legend = c("Optim", "PSO with 100 particles", 
                         "PSO with 1000 particles"), 
       col = c("red", "blue", "green3"), lty = 1, cex = 0.75)