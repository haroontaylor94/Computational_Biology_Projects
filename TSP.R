setwd("TSP")

#Function phi
phi <- function(d, K){
    exp((-d ^ 2) / (2 * K ^ 2))
}

#Euclidean distance
euclid_dist <- function(x, y){
    sqrt(sum((x- y)^2))
}

#Force pulling point towards city
force_1 <- function(cities, points, K){
    #Create vector of distances between points and cities
    dist_mat <- apply(cities, 1, function(x) 
        apply(points, 1, function(y) euclid_dist(x,y)))
    phi_mat <- phi(dist_mat, K)
    t(sapply(1:nrow(phi_mat), 
             function(y) colSums(apply(phi_mat, 
                                       2, 
                                       function(x) 
                                           x[y]/sum(x))*t(apply(cities, 
                                                                1, function(z) 
                                                                    z - points[y,])))))
}

#Function to create a circulant matrix
circ<-function(x) { 
    n <- length(x)
    suppressWarnings(matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
}

force_2 <- function(points){
    mat <- circ(c(-2, 1, rep(0, nrow(points) - 3), 1))
    mat %*% points
}

#Function for vectorised change in y
change_y <- function(cities, points, K, alpha = 0.2, beta = 2){
    alpha * force_1(cities, points, K) + beta * K * force_2(points)
}

#Function to plot cities and path through points
plot_cities_points_path <- function(cities, points, main){
    plot(cities[,1], cities[,2], xlab = NA, ylab = NA, 
         main = main, pch = 0, cex = 0.5, xaxt = "n", 
         yaxt = "n", bty = "n")
    lines(c(points[,1], points[1,1]), c(points[,2], points[1,2]), 
          col = "red")
}

#Function to run elastic net model
elastic_net <- function(cities, K, plot, n = 25){
    N <- nrow(cities)
    center <- c(mean(cities[,1]), mean(cities[,2]))
    angles <- seq(0, 2*pi, length.out = 2.5 * N + 1)[1:(2.5*N)]
    points <- matrix(c(center[1] + r * cos(angles), 
                       center[2] + r * sin(angles)), ncol = 2)
    if(plot){
        main = paste("K =", round(K[1], 4))
        plot_cities_points_path(cities, points, main)
    }
    for(k in K){
        for(i in 1:n){
            change_coords <- change_y(cities, points, k)
            points <- points + change_coords
            if(any(is.nan(points))){
                stop("Let's break!")
            }
        }
        if(plot & k %in% K_plot[-1]){
            main = paste("K =", round(k, 4))
            plot_cities_points_path(cities, points, main)
        }
    }
    return(points)
}

#Function to find order of cities
find_order <- function(cities, points){
    #Create vector of distances between points and cities
    dist_mat <- apply(cities, 1, function(x) 
        apply(points, 1, function(y) euclid_dist(x,y)))
    city_points <- apply(dist_mat, 2, which.min)
    names(city_points) <- 1:length(city_points)
    city_order <- as.numeric(names(sort(city_points)))
    city_order <- c(city_order, city_order[1])
}

#Function to find shortest distance through cities
calculate_dist <- function(cities, city_order){
    dists <- rep(NA, length(city_order)-1)
    for(i in 1:length(dists)){
        city <- cities[city_order[i], ]
        next_city <- cities[city_order[i+1], ]
        dists[i] <- euclid_dist(city, next_city)
    }
    sum(dists)
}

#Function to plot cities and path
plot_city_path <- function(cities, order, main, cex_main = 1.2){
    plot(cities[,1], cities[,2], xlab = NA, ylab = NA, 
         main = main, cex.main = cex_main, pch = 0, cex = 0.5, xaxt = "n", 
         yaxt = "n", bty = "n")
    lines(cities[order, 1], cities[order, 2], col = "red")
}

#Function to read in TSP file and perform elastic net
Solve_TSP <- function(tsp_file){
    TSP <- read.delim(tsp_file, skip = 6, header = FALSE, sep = "")
    TSP_cities <- as.matrix(TSP[-nrow(TSP), 2:3])
    TSP_cities_unit <- apply(TSP_cities, 2, function(x) x/max(x))
    TSP_points <- elastic_net(TSP_cities_unit, K, plot = FALSE)
    TSP_order <- find_order(TSP_cities_unit, TSP_points)
    TSP_dist <- calculate_dist(TSP_cities, TSP_order)
    return(list(cities = TSP_cities, order = TSP_order, dist = TSP_dist))
}

#Set radius of initial circle
r <- 0.1

#Set up K values
K <- cumprod(c(0.2, rep(0.99, 280)))

#Set up values of K for plotting
K_plot <- K[seq(1, length(K), length.out = 6)]

#Set up random set of cities
cities <- matrix(runif(100), ncol = 2)

#Run elastic net on random cities to produce plots
par(mfrow = c(2, 3))
points <- elastic_net(cities, K, plot = TRUE)

#Find order of random cities
order <- find_order(cities, points)

#Plot approximate and actual path of random cities
par(mfrow = c(1, 2), pty = "s", mar = c(0,1,3,1), oma = c(0,0,0,0))
plot_cities_points_path(cities, points, "Approximate path")
plot_city_path(cities, order, "Actual path")

#Solve TSP for examples with known optimal paths
eil <- Solve_TSP("eil51.tsp")
berlin <- Solve_TSP("berlin52.tsp")
st <- Solve_TSP("st70.tsp")
rat <- Solve_TSP("rat99.tsp")
bier <- Solve_TSP("bier127.tsp")
u <- Solve_TSP("u159.tsp")

#Plot elastic net solutions and known solutions
par(mar=c(1.5,1.5,1,1))
m <- matrix(c(1,1,1,1,1,1,2,3,3,3,3,3,3,4,5,5,5,5,5,5,
              6,6,6,7,7,7,8,9,9,9,10,10,10,11,12,12,12,13,13,13,
              14,14,14,14,14,14,15,16,16,16,16,16,16,17,18,18,18,18,18,18,
              19,19,19,19,19,19,20,21,21,21,21,21,21,22,23,23,23,23,23,23), 
            ncol = 20, byrow = TRUE)
layout(m, heights = c(1,2,1,2))
plot.new()
text(0.5,0.5,"51 cities",cex=2)
plot.new()
plot.new()
text(0.5,0.5,"52 cities",cex=2)
plot.new()
plot.new()
text(0.5,0.5,"70 cities",cex=2)
plot_city_path(eil$cities, eil$order, "Elastic net path", 0.8)
eil_opt <- read.delim("eil51.opt.tour", skip = 5, header = FALSE, sep = "", 
                      stringsAsFactors = FALSE)
eil_opt <- abs(as.numeric(eil_opt[-nrow(eil_opt), 1]))
plot_city_path(eil$cities, eil_opt, "Optimal path", 0.8)
plot.new()
plot_city_path(berlin$cities, berlin$order, "Elastic net path", 0.8)
berlin_opt <- read.delim("berlin52.opt.tour", skip = 4, header = FALSE, sep = "",
                         stringsAsFactors = FALSE)
berlin_opt <- abs(as.numeric(berlin_opt[-nrow(berlin_opt), 1]))
plot_city_path(berlin$cities, berlin_opt, "Optimal path", 0.8)
plot.new()
plot_city_path(st$cities, st$order, "Elastic net path", 0.8)
st_opt <- read.delim("st70.opt.tour", skip = 5, header = FALSE, sep = "", 
                     stringsAsFactors = FALSE)
st_opt <- abs(as.numeric(st_opt[-nrow(st_opt), 1]))
plot_city_path(st$cities, st_opt, "Optimal path", 0.8)
plot.new()
text(0.5,0.5,"99 cities",cex=2)
plot.new()
plot.new()
text(0.5,0.5,"127 cities",cex=2)
plot.new()
plot.new()
text(0.5,0.5,"159 cities",cex=2)
plot_city_path(rat$cities, rat$order, "Elastic net path")
plot.new()
plot_city_path(bier$cities, bier$order, "Elastic net path")
symbols(14500, 11500, circles = 1500, add = TRUE, inches = FALSE, fg = "blue")
plot.new()
plot_city_path(u$cities, unique(u$order), "Elastic net path")
symbols(4200, 4500, circles = 400, add = TRUE, inches = FALSE, fg = "blue")
symbols(8200, 3100, circles = 400, add = TRUE, inches = FALSE, fg = "blue")

#Calculate distances of elastic net paths
eil_dist <- calculate_dist(eil$cities, eil$order)
berlin_dist <- calculate_dist(berlin$cities, berlin$order)
st_dist <- calculate_dist(st$cities, st$order)
rat_dist <- calculate_dist(rat$cities, rat$order)
bier_dist <- calculate_dist(bier$cities, bier$order)
u_dist <- calculate_dist(u$cities, u$order)
