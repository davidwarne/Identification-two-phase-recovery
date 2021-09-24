##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File LTMPModllingTools.R
# Summary: Functions for fitting growth models and subsequent analysis. 
#
# Current features:
#   - Rate per unit cover vs cover curve and uncertainty quantification.
#   - Change-point regression to identify delays.
#
#-------------------------------------------------------------------------------
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021
#
#-------------------------------------------------------------------------------
#
# Note: assumes processed recovery trajectories as per LTMPDataTools.R function
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Standard growth models
#-------------------------------------------------------------------------------
## 
# @brief Exponential Growth Model
# @details Solution to the exponential growth model,
#           dC/dt = alpha C(t)
# where C(t) is the population density, alpha is the intrinsic growth rate.
#
# @param t time to evalutate solution
# @param alpha intrisic growth rate
# @param c0 initial population density at t = 0
#
# @returns the population density, C(t), at time t > 0
# @note for applications to coral growth, the implicit assumption is that % coral
# cover is directly proportional to population density, furthermore that growth is
# unbounded, therefore this model is only approriate in early growth phases.
exponential_sol <- function(t,alpha,c0) {
    return(c0*exp(alpha*t))
}

## 
# @brief Logistic Growth Model
# @details Solution to the standard logistic growth model,
#           dC/dt = alpha C(t)(1-C(t)/K)
# where C(t) is the population density, alpha is the intrinsic growth rate, and 
# K is the carrying capacity density.
#
# @param t time to evalutate solution
# @param alpha intrisic growth rate
# @param K carrying capacity density
# @param c0 initial population density at t = 0
#
# @returns the population density, C(t), at time t > 0
# @note for applications to coral growth, the implicit assumption is that % coral
# cover is directly proportional to population density
logistic_sol <- function(t,alpha,K,c0) {
    return((K*c0)/((K-c0)*exp(-alpha*t) + c0))

}

## 
# @brief Gompertz Growth Model
# @details Solution to the Gompertz growth model,
#           dC/dt = -alpha C(t)log(C(t)/K)
# where C(t) is the population density, alpha is the intrinsic growth rate, and 
# K is the carrying capacity density.
#
# @param t time to evalutate solution
# @param alpha intrisic growth rate
# @param K carrying capacity density
# @param c0 initial population density at t = 0
#
# @returns the population density, C(t), at time t > 0
# @note for applications to coral growth, the implicit assumption is that % coral
# cover is directly proportional to population density
gompertz_sol <- function(t,alpha,K,c0) {
    return(K*exp(log(c0/K)*exp(-alpha*t)))
}

##
# @brief Delayed Logistic Growth Model
# Two phase model of the form,
#   dC/dt = alpha1 C(t) if t <= T and dC/dt = alpha2 C(t)(1-C(t)/K) otherwise,
# where C(t) is the population density, alpha1 and alpha2 are the intrinsic growth 
#rates in phases 1 and 2,0 K is the carrying capacity density, and T is the first
# phase duration. Furthermore, alpha1 < alpha2 thereby representing a form of 
# delay in logistic growth
#
# @param t time to evalutate solution
# @param T phase 1 duration (i.e., delay duration)
# @param alpha1 intrisic growth rate in phase 1
# @param alpha2 intrisic growth rate in phase 2
# @param K carrying capacity density
# @param c0 initial population density at t = 0
#
# @returns the population density, C(t), at time t > 0
delay_logistic_sol <- function(t,T,alpha1,alpha2,K,c0) {
    if (t <= T) {
        return(exponential_sol(t,alpha1,c0))
    } else {
        return(logistic_sol(t-T,alpha2,K,exponential_sol(T,alpha1,c0)))
    }
}

##
# @brief Delayed Gompertz Growth Model
# Two phase model of the form,
#   dC/dt = alpha1 C(t) if t <= T and dC/dt = -alpha2 log(C(t)/K) otherwise,
# where C(t) is the population density, alpha1 and alpha2 are the intrinsic growth 
# rates in phases 1 and 2,0 K is the carrying capacity density, and T is the first
# phase duration. Furthermore, alpha1 < alpha2 thereby representing a form of 
# delay in Gompertz growth
#
# @param t time to evalutate solution
# @param T phase 1 duration (i.e., delay duration)
# @param alpha1 intrisic growth rate in phase 1
# @param alpha2 intrisic growth rate in phase 2
# @param K carrying capacity density
# @param c0 initial population density at t = 0
#
# @returns the population density, C(t), at time t > 0
delay_gompertz_sol <- function(t,T,alpha1,alpha2,K,c0) {
    if (t <= T) {
        return(exponential_sol(t,alpha1,c0))
    } else {
        return(gompertz_sol(t-T,alpha2,K,exponential_sol(T,alpha1,c0)))
    }
}

#-------------------------------------------------------------------------------
# Finite difference functions 
#-------------------------------------------------------------------------------

## ***
# Summary: First order forward finite difference of first derivative
#
# param: f function to be differentiated
# param: t time point to evaluate derivative at
# param: h step-size of discretisation
#
# return: O(h) approximation to df/dt at t
finite_diff_1stF <- function(f,t,h) {
    return((f(t + h) - f(t))/h)    
}

##  ***
# Summary: First order backward finite difference of first derivative
#
# param: f function to be differentiated
# param: t time point to evaluate derivative at
# param: h step-size of discretisation
#
# return: O(h) approximation to df/dt at t
finite_diff_1stB <- function(f,t,h) {
    return((f(t) - f(t - h))/h)    
}

##  ***
# Summary: Second order central finite difference of first derivative
#
# param: f function to be differentiated
# param: t time point to evaluate derivative at
# param: h step-size of discretisation
#
# return: O(h^2) approximation to df/dt at t
finite_diff_2ndC <- function(f,t,h) {
    return((f(t + h) - f(t - h))/(2.0*h))    
}

##  ***
# Summary: differentiate function f at points ti
# Detail: uses first order for boundary nodes and second order on internal nodes
# 
# param: f function to be differentiated
# param: ti time points to evaluate derivative at
# param: h step-size of discretisation
finite_diff <- function(f,ti,h) {
    n <- length(ti)
    dfdt <- numeric(n)
    # computing boundaries with first order finite differences
    dfdt[1] <- finite_diff_1stF(f,ti[1],h)
    dfdt[n] <- finite_diff_1stB(f,ti[n],h)
    # 2nd order finite differences for internal nodes
    for (k in 2:(n-1)){
        dfdt[k] <- finite_diff_2ndC(f,ti[k],h)
    }
    return(dfdt) 
}

##  ***
# Summary: computes quantile function of the ratio of two correlated Gaussian
# Detail: using the approximation by Hinkley (1969) for the distribution of W = X_1/X_2
# with X_1 ~ N(mu1,sig1^2) X_2 ~ N(mu2,sig2^2) with correlation rho. 
# Approximation assumption is that mu2 >> sig2 > 0 i.e. Pr(X_2 < 0) -> 0
#
# See Hinkley D.V. (1969) Biometrika v56 
#
# param: u point in quantile function to evaluate
# param: mu1 mean of numerator
# param: mu2 mean of denominator
# param: sig1 standard deviation of numerator
# param: sig2 standard deviation of denominator
# param: rho correlation coefficient between X_1 and X_2
ratio_quantile <- function(u,mu1,mu2,sig1,sig2,rho) {
   
    if (mu2 - 2*sig2 < 0) {
        # treat X(t) as constant for really small cover
        w <- qnorm(u,mu1/mu2,sig1/mu2)
        if (is.nan(w))
            print(c(u,mu1,mu2,sig1))
        return(w)
    } else {
        f <- function(w) {
            pnorm((mu2*w-mu1)/(sig1*sig2*
                  sqrt((w/sig1)^2 - 2*(rho*w)/(sig1*sig2) + (1/sig2)^2 ))) - u
        }
        res <- uniroot(f,lower = -1.0, upper = 1.0, f.lower = -u, f.upper = 1)
        return(res$root)
    }
}

## ***
# Summary: computes derivatives and rate per unit cover uncertainty quantiles
#
# param: f interpolation of coral cover
# param: f_sig interpolation of standard error in cover estimates
# param: h step-size to use
# 
# return: rate per unit cover estimates along with 80% credible intervals
finite_diff_pc_uq <- function(f,f_sig,ti,h) {
    n <- length(ti)
    pc <- numeric(n)
    lb <- numeric(n)
    ub <- numeric(n)
    med <- numeric(n)

    # computing boundaries with first order finite differences
    mu1 <- finite_diff_1stF(f,ti[1],ti[2]-ti[1])
    sig1 <- sqrt((f_sig(ti[2])^2 + f_sig(ti[1])^2)/(ti[2]-ti[1])^2)
    
    mu2 <- f(ti[1])
    sig2 <- f_sig(ti[1])
    
    rho <-  -1*(f_sig(ti[1])/sqrt(f_sig(ti[2])^2 + f_sig(ti[1])^2))
    
    pc[1] <- mu1/mu2
    lb[1] <- ratio_quantile(0.1, mu1,mu2,sig1,sig2,rho)
    ub[1] <- ratio_quantile(0.9, mu1,mu2,sig1,sig2,rho)
    
    mu1 <- finite_diff_1stB(f,ti[n],ti[n]-ti[n-1])
    sig1 <- sqrt((f_sig(ti[n])^2 + f_sig(ti[n-1])^2)/(ti[n]-ti[n-1])^2)
    
    mu2 <- f(ti[n])
    sig2 <- f_sig(ti[n])
    
    rho <-  (f_sig(ti[n])/sqrt(f_sig(ti[n])^2 + f_sig(ti[n-1])^2))
    
    pc[n] <- mu1/mu2
    lb[n] <- ratio_quantile(0.1, mu1,mu2,sig1,sig2,rho)
    ub[n] <- ratio_quantile(0.9, mu1,mu2,sig1,sig2,rho)
    
    # 2nd order finite differences for internal nodes
    for (k in 2:(n-1)){
        h <- min(ti[k]-ti[k-1],ti[k+1]-ti[k])
        mu1 <- finite_diff_2ndC(f,ti[k],h)
        sig1 <- sqrt((f_sig(ti[k]+h)^2 + f_sig(ti[k]-h)^2)/(2*h)^2)
        
        mu2 <- f(ti[k])
        sig2 <- f_sig(ti[k])
        
        rho <-  0
        
        pc[k] <- mu1/mu2
        lb[k] <- ratio_quantile(0.1, mu1,mu2,sig1,sig2,rho)
        ub[k] <- ratio_quantile(0.9, mu1,mu2,sig1,sig2,rho)
    }
   
    return(list(pc,lb,med,ub))
}

## ***
# Summary: computes log likelihood for a linear regression model with Gaussian
# residuals
# 
# param: model fitted linear model (e.g., model <- lm(y ~ x))
# return: log likelihood using unbiased esimate of variance
lmLogLike <- function(model) {
    N <- length(model$residuals)
    SSE <- sum(model$residuals^2)
    return(-(N-2)/2 - N*log(2*pi)/2 - N*log(SSE/(N-2))/2)
}

## ***
# Summary: change point linear regression using supermum of the likelihood ratio
# Detail: computes change point regression for all possible data partitions
# selects the one with highest likelihood ratio with respect to single linear 
# regression.
#
# param: traj - recovery trajectory time series data for a single site recovery 
#               period.
# return: a list containing the following elements: the data index of the change
#         point, the two linear models, result of Chow test, R^2 and adjusted R^2
#         for change point regression, and single linear regression model. 
sup_Like_Ratio <- function(traj) {

    N <- length(traj$HC)
    ind <- 1:N
    traj <- traj %>% mutate(IND=ind) %>% arrange(HC)
    llr <- numeric(N)
    llr[1:N] <- -Inf 
    for (k in 3:(N-2)) {
        # 2. for each possible change point
        # 2.1 fit two linear models to split data [1,k],[k,N]
        traj1 <-traj %>% slice(1:k)
        traj2 <- traj %>% slice(k:N)
        y1 <- traj1$pcHC_lin_fd
        x1 <- traj1$HC
        mod1 <- lm(y1 ~ x1)
        y2 <- traj2$pcHC_lin_fd
        x2 <- traj2$HC
        mod2 <- lm(y2 ~ x2)
        # 1. fit linear model to entire data
        y <- c(y1,y2)
        x <- c(x1,x2)
        mod <- lm(y ~ x)
        # 2.2 compute log likelihood ratio (with possible penalty?)
        llr[k] <- lmLogLike(mod1) + lmLogLike(mod2) - lmLogLike(mod)
    }
    
    # 3. find K = argmax loglike ratio
    K <- which.max(llr)
    traj1 <-traj %>% slice(1:K)
    traj2 <- traj %>% slice(K:N)
    y1 <- traj1$pcHC_lin_fd
    x1 <- traj1$HC

    mod1 <- lm(y1 ~ x1)
    y2 <- traj2$pcHC_lin_fd
    x2 <- traj2$HC
    mod2 <- lm(y2 ~ x2)
    
    # 1. fit linear model to entire data
    y <- c(y1,y2)
    x <- c(x1,x2)
    mod <- lm(y ~ x)

    # coefficient of determination for each model
    R2_cp <- 1 - (sum(mod1$residuals^2) + 
                  sum(mod2$residuals^2))/(sum((y-mean(y))^2))
    adjR2_cp <- 1 - (1 - R2_cp)*(length(y)-1)/(length(y)-2) 

    # 4. apply Chow test to obtain F-stat and p-value for change-point
    cpt <- chow.test(y1,x1,y2,x2)
    
    # 5. return both change-point model and single linear model
    return(list(traj$IND[K],mod1,mod2,cpt,R2_cp,adjR2_cp,mod,
                traj$IND[K-1],traj$IND[K+1]))
}

# finds the intersection point of two models
compute_intersection <- function(model1,model2) {

    xx <- optimize(function(x) abs(predict(model1,list(x1 = x))
                                   - predict(model2,list(x2 = x))),
                   c(0.0, 100.0))$minimum
    return(xx)
}

## ***
# Summary: estimates the effect of the two-phase recovery
# Details: fits standard growth model to second phase of two-phase trajectory 
# then evaluate the % cover after a fixed interval and the number of years required 
# to reach given cover assuming no-delay. This can then be paired with the actual 
# delayed data. 
#
# param: traj - processed trajectory including group and benthos codes
# param: cp_ind - index of change point assuming data sorted by visit number
# param: TInterval - the length of time to project cover for
# param: CoverThresh - % cover to determine time to reach this target
quantify_delay_impact <- function(traj,cp_ind, TInterval,CoverThresh){

    traj <- traj %>% arrange(VISIT_NO)
    Nobs <- nrow(traj)
    
    # estimate carrying capacity K = 100 - AB (abiotic) + ST (silt transient)
    K <- 100.0 - traj$AB_G[Nobs] + traj$ST_B[Nobs]

    # filter trajectory for post-change point period
    y_obs <- traj$HC_G[cp_ind:Nobs]
    sigma <- traj$HC_se_G[cp_ind:Nobs]
    
    # reset time to cp and scale time to order of years
    t_obs <- (traj$T[cp_ind:Nobs] - traj$T[cp_ind])/365.0
    
    # Apply MLE using Gompertz and logistic models
    
    # neg-loglike for Gompertz model
    Gompertz_nll <- function(alpha) {
        hc_mu <- y_obs
        if (alpha > 0) {
            for (i in 2:length(t_obs)){
                hc_mu[i] <- gompertz_sol(t_obs[i],alpha,K,y_obs[1])
            }
            return(-sum(dnorm(y_obs,hc_mu,sigma,log = TRUE)))
        } else {
            return(NA)
        }
    }
    # neg-loglike for logistic model
    Logistic_nll <- function(alpha) {
        hc_mu <- y_obs
        if (alpha > 0) {
            for (i in 2:length(t_obs)){
                hc_mu[i] <- logistic_sol(t_obs[i],alpha,K,y_obs[1])
            }
            return(-sum(dnorm(y_obs,hc_mu,sigma,log = TRUE)))
        } else {
            return(NA)
        }
    }

    # compute MLE for both models
    fit_G <- mle(Gompertz_nll, start = list(alpha = 0.5),nobs = length(y_obs))
    fit_L <- mle(Logistic_nll, start = list(alpha = 0.5),nobs = length(y_obs))
    
    # choose best model
    if (AIC(fit_G) < AIC(fit_L)) {
        mle_alpha <- coef(fit_G)
        model_func <- gompertz_sol
    } else {
        mle_alpha <- coef(fit_L)
        model_func <- logistic_sol
    }

    cover_res <- c()
    time_res <- c()

    # return % cover after X years 
    cover_res[1] <- model_func(TInterval,mle_alpha,K,traj$HC_G[1])

    # return time to reach Y % cover
    f <- function(t) model_func(t,mle_alpha,K,traj$HC_G[1]) - CoverThresh
    res <- uniroot(f,lower = 0.0, upper = 20.0, f.lower = - CoverThresh, 
                   f.upper = K - CoverThresh)
    time_res[1] <- res$root

    # same thing, but taking the delay offset into account
    delay_cover_res <- c()
    delay_time_res <- c()
    if (TInterval > traj$T[cp_ind]/365.0) {
        delay_cover_res[1] <- model_func(TInterval - traj$T[cp_ind]/365.0,
                                         mle_alpha,K,traj$HC_G[cp_ind])
    } else {
        # get the closest data point to the TInterval
        I <- which.min(abs(traj$T[cp_ind]/365.0 - TInterval))
        delay_cover_res[1] <- traj$HC_G[I]
    }
    f <- function(t) model_func(t,mle_alpha,K,traj$HC_G[cp_ind]) - CoverThresh
    res <- uniroot(f,lower = 0.0, upper = 20.0, f.lower = - CoverThresh, 
                   f.upper = K - CoverThresh)
    delay_time_res[1] <- res$root + traj$T[cp_ind]/365.0

    return(list(cover_res,time_res,delay_cover_res,delay_time_res))
}
