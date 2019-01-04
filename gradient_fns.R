#Calculate gradients between consecutive points
grad <- function(x, y){
    diff(y) / diff(x)
}

#Calculate gradients between each point and all subsequent points
grad2 <- function(x, y){
    dist_x <- sapply(x, "-", x)
    dist_y <- sapply(y, "-", y)
    grads <- dist_y / dist_x
    grads[upper.tri(grads)]
}

#Calculate the gradient of the line of best fit between the points
lm_grad <- function(x, y){
    unname(lm(y~x)$coefficients[2])
}

#Calculate overlap between points
overlap <- function(x){
    range <- diff(range(x))
    total_length <- sum(abs(diff(x)))
    total_length - range
}

#Find outliers
outlier <- function(x, y, thresh = 0.15){
    fit <- lm(y ~ x)
    grad <- fit$coefficients[2]
    cooksd <- cooks.distance(fit)
    cook_bool <- cooksd < 4*mean(cooksd)
    x <- na.omit(x)[cook_bool]
    y <- na.omit(y)[cook_bool]
    lm_grad(x, y) < (grad / 2)
}

#Remove outliers
remove_outliers <- function(x, y){
    fit <- lm(y ~ x)
    cooksd <- cooks.distance(fit)
    cook_bool <- cooksd > 4*mean(cooksd)
    x[cook_bool] <- NA
    y[cook_bool] <- NA
    return(c(x, y))
}