library(deSolve)

#Hodgkin Huxley model
HH_ODE <- function(time, init, parms){
    with(as.list(c(init, parms)),{
        #alpha n
        an <- function(V) 0.01*(V+55)/(1-exp(-0.1*(V+55))) 
        #beta n
        bn <- function(V) 0.125*exp(-0.0125*(V+65))
        
        #alpha m
        am <- function(V) 0.1*(V+40)/(1-exp(-0.1*(V+40)))
        #beta m
        bm <- function(V) 4*exp(-0.0556*(V+65))
        
        #alpha h
        ah <- function(V) 0.07*exp(-0.05*(V+65))
        #beta h
        bh <- function(V) 1/(1+exp(-0.1*(V+35)))
        
        #dv/dt
        dv <- (I-(gl*(V-El)+gk*n^4*(V-Ek)+gna*m^3*h*(V-Ena)))/cm
        #dn/dt
        dn <- an(V)*(1-n)-bn(V)*n
        #dm/dt
        dm <- am(V)*(1-m)-bm(V)*m
        #dh/dt
        dh <- ah(V)*(1-h)-bh(V)*h
        return(list(c(dv, dm, dh, dn)))
    })
}

#Initial values
init <- c(V=-65, m=0.0529, h=0.5961, n=0.3177)
#Times
times <- seq(0,50,0.1)
#Parameters
parms <- c(Ena=50, Ek=-77, El=-54.402, gna=1200, gk=360, gl=3, cm=10, I=200)

HH <- lsoda(init, times, HH_ODE, parms)

png("Hodgkin_Huxley_200.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2))
plot(HH[,"time"], HH[,"V"], type = "l", xlab = "time (ms)", ylab = "V (mV)")
plot(HH[,"time"], HH[,"m"], type = "l", xlab = "time (ms)", ylab = "m", 
     ylim = c(0,1))
plot(HH[,"time"], HH[,"h"], type = "l", xlab = "time (ms)", ylab = "h", 
     ylim = c(0,1))
plot(HH[,"time"], HH[,"n"], type = "l", xlab = "time (ms)", ylab = "n",
     ylim = c(0,1))
dev.off()

#Calculate firing rate
firing_rate <- function(HH){
    x <- HH[,2]
    peaks <- which(diff(sign(diff(x * (x > 0)))) < 0) + 2
    peaks_s <- HH[peaks,1]/1000
    freq <- mean(diff(peaks_s))
    if (is.nan(freq)){
        rate <- 0
    } else {
        rate <- 1/freq
    }
    return(rate)
}

firing_rates <- rep(0, 500)
for(I in seq(0, 500)){
    parms <- c(Ena=50, Ek=-77, El=-54.402, gna=1200, gk=360, gl=3, cm=10, I=I)
    HH <- lsoda(init, times, HH_ODE, parms)
    rate <- firing_rate(HH)
    firing_rates[I] <- rate
    print(I)
    print(rate)
}

png("Firing_rate_vs_current.png", width = 4, height = 4, units = "in", 
    res = 300)
par(mfrow = c(1,1))
plot(1:500, firing_rates, type = "l", 
     xlab = expression("I"[e]*"/A (nA/mm"^2*")"), 
     ylab = "Firing rate (Hz)")
dev.off()

#Negative current of -50
#Initial values
init <- c(V=-65, m=0.0529, h=0.5961, n=0.3177)
#Times
times <- seq(0,5,0.1)
#Parameters
parms <- c(Ena=50, Ek=-77, El=-54.402, gna=1200, gk=360, gl=3, cm=10, I=-50)

HH_neg <- lsoda(init, times, HH_ODE, parms)

#Followed by current of 0
#Initial values
init <- HH_neg[nrow(HH_neg), -1]
#Times
times <- seq(5,50,0.1)
#Parameters
parms <- c(Ena=50, Ek=-77, El=-54.402, gna=1200, gk=360, gl=3, cm=10, I=0)

HH_0 <- lsoda(init, times, HH_ODE, parms)

HH <- rbind(HH_neg, HH_0[-1, ])
png("Hodgkin_Huxley_minus50_0.png", width = 8, height = 8, units = "in", 
    res = 300)
par(mfrow = c(2,2))
plot(HH[,"time"], HH[,"V"], type = "l", xlab = "time (ms)", ylab = "V (mV)")
plot(HH[,"time"], HH[,"m"], type = "l", xlab = "time (ms)", ylab = "m")
plot(HH[,"time"], HH[,"h"], type = "l", xlab = "time (ms)", ylab = "h")
plot(HH[,"time"], HH[,"n"], type = "l", xlab = "time (ms)", ylab = "n")
dev.off()
