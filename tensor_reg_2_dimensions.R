tensor.reg <- function(z.train, x.train, y.train, nsweep, rank, burn = 0.30, nskip = 3, a.lam, b.lam, phi.alpha, scale = TRUE) {
    library(plyr)
    library(GIGrvg)
    library(gtools)
    library(coda)
    library(fields)
        
    n <- length(y.train)            ## sample size
    p <- unique(dim(x.train)[-1])   ## tensor is of dimension p^d
    d <- length(dim(x.train)) - 1   
    pgamma <- ncol(z.train)         ## number of scalar predictors
    
    #### standarize ####
    mz <- colMeans(z.train)         
    sz <- apply(z.train, 2, function(z) diff(range(z)))
    Zt <- z.train ## scalar predictors
    
    my <- mean(y.train)  ## mean of responses
    sy <- ifelse(scale, sd(y.train), 1) ## sd of responses
    obs <- as.numeric(scale(y.train, center = my, scale = sy)) ## scaled response
    
    if(scale){
        ##centering & scale to have unity range
        Xt <- 0 * x.train
        mx <- apply(x.train, c(2:(d+1)), function(z) mean(z))
        sx <- apply(x.train, c(2:(d+1)), function(z) diff(range(z)))
        sx[sx == 0] <- 1
        for(jj in 1:nrow(x.train)) Xt[jj,,] <- (x.train[jj,,] - mx) / sx
    } else{
        ## do nothing;
        mx <- array(0, dim = dim(x.train)[-1])
        sx <- array(1, dim = dim(x.train)[-1])
        Xt <- x.train  ## tensor predictor
    }
        
    #### MCMC setup ####
    ZZ <- crossprod(Zt, Zt)

    require(glmnet)
    vecXt <- array(x.train, dim = c(length(y.train), p^d)) ## vectorize the tensor predictor
    vecXt <- cbind(vecXt, z.train) 

    #### Initialize the vector and the tensor coefficients ####
    #### lasso set up with response on scalar and vectorized tensor predictors ####
    las <- cv.glmnet(x = vecXt, y = y.train) 
    las <- glmnet(x = vecXt, y = y.train, lambda=las$lambda.min)
    beta.init <- as.numeric(las$bet)
    gam <- beta.init[1:pgamma]  ## initialized scalar predictor coefficients
    
    ## hyper-parameter initialize 
    if(missing(a.lam)) a.lam <- rep(3, rank)
    if(missing(b.lam)) b.lam <- (a.lam)**(1/(2 * d))
    if(missing(phi.alpha)) phi.alpha <- rep(1 / rank, rank)
    
    phi.a0 <- sum(phi.alpha)
    a.vphi <- phi.a0
    b.vphi <- phi.alpha[1] * rank^(1/d)
    
    c0 <- 0  
    s0 = 1; a.t = 2.5/2; b.t = 2.5/2 * s0^2
    tau2  <- 1 / rgamma(1, a.t, b.t) ## initialize tau2
    
    phi   <- rdirichlet(1, phi.alpha)    ## initialize phi
    varphi <- rgamma(1, a.vphi, b.vphi)  ## initialize varphi
    tau.r <- phi * varphi                ## initialize tau.r
    
    lambda <- matrix(rgamma(rank*d, a.lam[1], b.lam[1]), rank, d)            ## initialize lambda
    omega <- array(rexp(rank*p*d, 0.5*(a.lam[1]/b.lam[1])), dim=c(rank,d,p)) ## initialize omega

    hhhh <- list()
    hhhh1 <- list()
    svd.n <- svd(matrix(beta.init[-c(1:pgamma)],70,70))
    beta  <- replicate(rank, list(matrix(rnorm(p*d), p, d)))

    for(i in 1:rank){
        beta[[i]] <- cbind(svd.n$u[,i],svd.n$v[,i])*sqrt(svd.n$d[i])  ## initialize tensor coefficient margins 
    }

    lambda <- matrix(rgamma(rank*d, a.lam[1], b.lam[1]), rank, d)
    omega <- array(rexp(rank*p*d, 0.5*(a.lam[1]/b.lam[1])), dim=c(rank,d,p))
    beta  <- replicate(rank, list(matrix(rnorm(p*d), p, d)))   
    
    alpha.store <- rep(NA,nsweep)
    c0.store    <- rep(NA,nsweep)
    gam.store   <- array(data=NA,dim=c(nsweep,pgamma))
    tau2.store  <- rep(NA,nsweep)
    phi.store   <- array(data=NA,dim=c(nsweep,rank))
    varphi.store <- array(data=NA,dim=c(nsweep,1))
    beta.store  <- array(data=NA,dim=c(nsweep,rank,p,d))
    omega.store <- array(data=NA,dim=c(nsweep,rank,d,p))
    lambda.store   <- array(data=NA,dim=c(nsweep,rank,d))
    hyppar.store <-array(data=NA,dim=c(nsweep,rank,2)) 
    
    par.grid <- expand.grid(alam=seq(2.1,d+1,length.out=5), zeta=seq(0.5,ceiling(10 * rank**(1/(2*d))/2)/10,length.out=5))
    
    alpha.grid <- seq(rank**(-d), rank**(-0.1), length.out = 10); M <- 20   ## grid of alpha values
    score.store <- array(data=NA,dim=c(nsweep,length(alpha.grid)))
        
    #### MCMC run ####
    tt <- Sys.time()
    for(sweep in 1:nsweep){
        tens.mean <- getmean(Xt, beta, rank)
        
        ## update (a.lam, b.lam)
        Cjr <- sapply(1:rank, function(rr){bb <- sapply(1:d, function(jj) sum(abs(beta[[rr]][,jj]))); bb <- bb / sqrt(tau.r[rr]); return(bb)})
        mfun <- function(z,rank){o <- sapply(Cjr[,rank], function(val) return(lgamma(z[1]+p) - lgamma(z[1]) + z[1]*log(z[2]*z[1]) - (z[1]+p)*log(z[2]*z[1] + val))); return(sum(o))}
        ll <- sapply(1:rank, function(rr) apply(par.grid, 1, mfun, rank = rr))
        par.wt <- apply(ll, 2, function(z) return(exp(z - logsum(z))))
        ixx <- apply(par.wt, 2, sample, x = c(1:nrow(par.grid)), size = 1, replace = F)
        for(rr in 1:rank){
            a.lam[rr] <- par.grid[ixx[rr],1]
            b.lam[rr] <- par.grid[ixx[rr],2] * a.lam[rr]
        }
        
        ## update gamma (scalar predictor coefficients)
        Sig.g <- chol2inv(chol(diag(pgamma) + ZZ / tau2))
        mu.g <- Sig.g %*% (crossprod(Zt,obs-c0-tens.mean) / tau2)
        gam <- mu.g + chol(Sig.g) %*% rnorm(pgamma)
        
        ## update alpha (intercept)
        pred.mean <- Zt %*% gam
        mu.c0 <- mean(obs-pred.mean-tens.mean)
        c0 <- rnorm(1, mean=mu.c0, sd=sqrt(tau2 / n))
        
        ## update tau2 
        a.tau <- a.t + n / 2
        b.tau <- b.t + 0.5 * sum((obs-c0-pred.mean-tens.mean)^2)
        tau2 <- 1 / rgamma(1, a.tau, b.tau)
        
        ## update (alpha, phi, varphi)
        draw.phi_tau <- function(alpha){
            len <- length(alpha)
            
            m.phialpha <- rep(alpha[1], rank)
            m.phia0 <- sum(m.phialpha)
            m.avphi <- m.phia0
            ## assumes b.vphi const (use: alpha 1 / R)
            
            Cr <- sapply(1:rank, function(rr) {bb <- sapply(1:d, function(jj) crossprod(beta[[rr]][,jj], diag(1/omega[rr,jj,]) %*% beta[[rr]][,jj])); return(bb)})
            
            score.fn <- function(phi.alpha, phi.s, varphi.s, Cstat){
                ldirdens <- function(v, a){
                    c1 <- lgamma(sum(a))
                    c2 <- sum(lgamma(a))
                    return((c1-c2) + sum((a-1) * log(v)))
                }
                ldir <- apply(phi.s, 1, ldirdens, a = phi.alpha)
                
                lvarphi <- dgamma(varphi.s, sum(phi.alpha), b.vphi, log = T)
                
                dnorm.log <- -rowSums(Cstat) / (2 * varphi.s) -(p*d/2) * sapply(1:length(varphi.s), function(ii)return(sum(log(varphi.s[ii] * phi.s[ii,]))))
                return(dnorm.log + ldir + lvarphi)
            }
            
            phi <- NULL; varphi <- NULL; scores <- NULL
            if(len > 1){
                phi <- matrix(0, M*length(alpha.grid), rank)
                varphi <- matrix(0, M*length(alpha.grid), 1)
                Cstat <- matrix(0, M*length(alpha.grid), rank)
                scores <- list()
                
                ## get reference set
                for(jj in 1:len){  
                    m.phialpha <- rep(alpha[jj], rank)
                    m.phia0 <- sum(m.phialpha)
                    m.avphi <- m.phia0
                    
                    ## draw phi
                    Cr1 <- colSums(Cr)
                    phi.a <- sapply(1:rank, function(rr){rgig(M,m.phialpha[rr]-p*d/2,Cr1[rr],2*b.vphi)})
                    phi.a <- t(apply(phi.a, 1, function(z)return(z / sum(z)))) ## [M x rank]
                    
                    ## draw varphi ##colSums(Cr / t(replicate(d, z)))
                    Cr2 <- t(apply(phi.a, 1, function(z)return(Cr1 / z)))
                    varphi.a <- apply(Cr2, 1, function(z)return(rgig(1, m.avphi-rank*p*d/2, sum(z), 2*b.vphi)))
                    phi[seq((jj-1)*M+1, jj*M), ] <- phi.a
                    varphi[seq((jj-1)*M+1, jj*M)] <- varphi.a
                    Cstat[seq((jj-1)*M+1, jj*M), ] <- Cr2                
                }
                scores <- lapply(alpha.grid, function(z)return(score.fn(rep(z,rank), phi, varphi, Cstat)))
                lmax <- max(unlist(scores))
                scores <- sapply(scores, function(z)return(mean(exp(z - lmax))))
            } else{
                ## draw phi
                Cr1 <- colSums(Cr)
                phi <- sapply(1:rank, function(rr){rgig(1,m.phialpha[rr]-p*d/2,Cr1[rr],2*b.vphi)})
                phi <- phi / sum(phi)
                
                ## draw varphi
                Cr2 <- Cr1 / phi
                varphi <- rgig(1, m.avphi-rank*p*d/2, sum(Cr2), 2*b.vphi)
            }
            return(list(phi=phi, varphi=varphi, scores=scores))
        }
        
        ## sample astar
        o <- draw.phi_tau(alpha.grid)
        astar <- sample(alpha.grid, size = 1, prob = o$scores)
        cat(sprintf('scores: %s\n', paste(round(score <- o$scores/sum(o$scores),2),collapse = ', ')))
        score.store[sweep,] <- score
        
        ## sample (phi, varphi)
        o <- draw.phi_tau(astar)
        phi <- o$phi; varphi <- o$varphi
        tau.r <- varphi * phi
        phi.alpha <- rep(astar, rank); phi.a0 <- sum(phi.alpha); a.vphi <- phi.a0
                
        ## update rank specific params
        for(r in 1:rank){
            for(j in 1:d){
                tens.mu.r <- getmean(Xt, beta, rank, r)
                
                betj <- getouter(beta[[r]],j)
                H <- matrix(NA, n, p)
                for(i in 1:n) {
                    H[i, ] <- apply(Xt[i,,], j, function(x)return(sum(x * betj))) 
                }
                
                HH <- crossprod(H,H)
                K <- chol2inv(chol(HH / tau2 + diag(1/omega[r,j,]) / tau.r[r]))
                
                ## update betas
                mm <- (obs-c0-pred.mean-tens.mu.r)
                bet.mu.jr <- K %*% crossprod(H / tau2, mm)
                beta[[r]][,j] <- bet.mu.jr + chol(K) %*% rnorm(p)
                
                ## update lambda.jr
                lambda[r,j] <- rgamma(1, a.lam[r] + p, b.lam[r] + sum(abs(beta[[r]][,j])) / sqrt(tau.r[r]))
                
                ## update omega.jr
                omega[r,j,] <- sapply(1:p, function(kk) rgig(1, 1/2, beta[[r]][kk,j]^2 / tau.r[r], lambda[r,j]^2))
                #omega[r,j,] <- sapply(1:p, function(kk){a <- lambda[r,j]^2; b <- beta[[r]][kk,j]^2 / tau.r[r]; map <- besselK(sqrt(a*b),0.5 + 1) / besselK(sqrt(a*b), 0.5) * sqrt(b / a); return(map)})
            }
            
            beta.store[sweep,r,,] <- beta[[r]]
        }
        
        ## store params
        tau2.store[sweep] <- tau2
        c0.store[sweep]   <- c0
        gam.store[sweep,] <- gam
        alpha.store[sweep] <- astar ## not intercept
        phi.store[sweep,] <- phi
        varphi.store[sweep,] <- varphi
        omega.store[sweep,,,] <- omega
        lambda.store[sweep,,] <- lambda
        sapply(1:rank, function(rr) hyppar.store[sweep,rr,] <<- c(a.lam[rr], b.lam[rr]))
        
        if(sweep %% 5 ==0) cat(sprintf('%i, %s: %2.3f, %s: %2.3f %2.3f %2.3f\n', sweep, "tau2", tau2 * sy^2, "(alpha,a.lam,b.lam)", astar, a.lam[r], b.lam[r]))
    }
    
    tt <- abs(tt - Sys.time())
    cat('Time out: ', tt, '\n')
            
    #### finalize ####    
    out <- list(nsweep = nsweep, rank = rank, p = p, d = d, par.grid = par.grid, alpha.grid = alpha.grid, my = my, sy = sy, mz = mz, sz = sz, mx = mx, sx = sx, Zt = Zt, Xt = Xt, obs = obs, a.t = a.t, b.t = b.t, tau2.store = tau2.store, c0.store = c0.store, gam.store = gam.store, alpha.store = alpha.store, beta.store = beta.store, phi.store = phi.store, varphi.store = varphi.store, omega.store = omega.store, lambda.store = lambda.store, hyppar.store = hyppar.store, score.store = score.store, time = tt)
    
    class(out) <- "tensor.reg"
    return(out)
}

#### aux functions ####
getouter <- function(bet,j=NULL){
    idx <- setdiff(1:ncol(bet),j)
    out <- bet[,idx[1]]
    
    idx <- setdiff(idx,idx[1])
    for(k in idx) {
        out <- outer(out, bet[,k])
    }
    return(out)
}

getmean <- function(X,beta,rank,rank.exclude=NULL){    
    idx <- setdiff(1:rank,rank.exclude)
    B <- Reduce('+', lapply(beta[idx], getouter))
    mu.B <- apply(X, 1, function(xx, bb) sum(xx * bb), bb = B)
    return(mu.B)
}

# expand.grid <- function(x, y){
#     nx <- ifelse(is.null(dim(x)), length(x), nrow(x))
#     ny <- ifelse(is.null(dim(y)), length(y), nrow(y))
#     cbind(kronecker(rep(1, ny), x), kronecker(y, rep(1, nx)))
# }

uncollapse <- function(str, collapse = "", mode = "character"){
    a <- unlist(strsplit(str, collapse))
    mode(a) <- mode
    return(a)
}

rmse <- function(a, b) return(sqrt(mean((a - b)^2)))

logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

