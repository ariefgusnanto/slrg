# List of functions to perform sparse logistic regression
# with hierarchical likelihood
# This file contains functions for the next updated
# slrg package

# Here are the updates:
# The function log.sparse now implements multi-core
# and vectorised computation.

####################################################################
####################################################################
####################################################################
####################################################################

merge.fun <- function(a,b){
# merge the results of log.sparse()
# A.Gusnanto@leeds.ac.uk
names.a <- names(a)
res <- list()
for(j in 1:(length(a)-1)){
res[[j]] <- c(a[[j]], b[[j]])
}
estimates <- list()
names.est <- names(a$estimates)
for(j in 1:length(a$estimates)){
estimates[[j]] <- c(a$estimates[[j]], b$estimates[[j]])
}
names(estimates) <- names.est
res <- c(res, list(estimates))
names(res) <- names.a
return(res)
}# end of function



merge.fun.res <- function(a,b){
# merge the results cv.log.sparse()
# A.Gusnanto@leeds.ac.uk
AIC.a <- min(a$estimates$AIC.vector)
AIC.b <- min(b$estimates$AIC.vector)

names.a <- names(a)
res <- list()
for(j in 1:(length(a)-1)){
if(AIC.a <= AIC.b){
res[[j]] <- a[[j]]
} else {
res[[j]] <- b[[j]]
}}
estimates <- list()
names.est <- names(a$estimates)
for(j in 1:length(a$estimates)){
if("numeric" %in% is(a$estimates[[j]])){
estimates[[j]] <- c(a$estimates[[j]], b$estimates[[j]])
 } else {
 estimates[[j]] <- cbind(a$estimates[[j]], b$estimates[[j]])
 }
}
names(estimates) <- names.est
res <- c(res, list(estimates))
names(res) <- names.a
return(res)
}



plot.cv <- function(x,xlab=expression(paste("Log ", lambda)), ylab="Cross validation error", lwd=1.5, ylim=c(0, 0.5), las=1, ...){
n <- length(x$estimates$y.val.fold)
plot(x$log.lambda.seq, x$error.seq/n, xlab=xlab, ylab=ylab, ylim=ylim, type="n", las=las, ...)
abline(h=seq(0,0.5, by=0.1), v=seq(-10, 10, by=2), lty=3, col=2)
lines(x$log.lambda.seq, x$error.seq/n, lwd=lwd)
points(x$log.lambda.seq, x$error.seq/n, pch=19)
}# end of plot.cv


plot.sparse <- function(x, type="l", ...){
# OLD style, will be deprecated. See plot.lambda()
plot(x$estimates$log.lambda, x$estimates$AIC.vector, type="l", ...)
lines(x$estimates$log.lambda, x$estimates$BIC.vector, col=2, ...)
}# end of function


plot.lambda <- function(x, crit="BIC", type="l", las=1, ...){
if(crit=="AIC"){
plot(x$estimates$log.lambda, x$estimates$AIC.vector, type=type, las=las, ...)
} else {
plot(x$estimates$log.lambda, x$estimates$BIC.vector, type=type, las=las, ...)
}
}



#     plot b
plot.re <- function(res, chrom=NULL, crit="AIC", ylab="", cex.axis=0.8, ...){
if(class(res)[1]=="numeric"){
 bhat <- res
} else {
  if(crit=="BIC"){
   id <- which.min(res$estimates$BIC.vector)
   bhat <- res$estimates$b.matrix[,id] 
    } else {
   bhat <- res$bhat
    }
}

if(!is.null(chrom)){
   chrom <- as.numeric(gsub("chr", "", chrom))
   chrom.id <- unique(chrom)
   id <- 1:length(chrom)
   tick <- round(sapply(split(id, chrom), median))
   plot(id, bhat, type="h", ylab=ylab, axes=F, ...)
   axis(1, at=tick, chrom.id, cex.axis=cex.axis)
   axis(2, las=1, cex.axis=cex.axis)
   box()
   mtext(expression(paste(hat(b))), side = 2, line = 3, las = 1)
 } else {
   plot(bhat, type="h", ylab=ylab, cex.axis=cex.axis, ...)
 }
} # end of function plot.re



ann.out <- function(res, ann, crit="AIC", threshold=1e-4){
  if(crit=="BIC"){
   id <- which.min(res$estimates$BIC.vector)
   bhat <- res$estimates$b.matrix[,id] 
    } else {
   bhat <- res$bhat
   }
   bhat[abs(bhat)<=threshold] <- 0
   nonzero.id <- which(abs(bhat)>0)
   bhat.nonzero <- bhat[nonzero.id]
   ann.nonzero <- cbind(ann[nonzero.id,], bhat=bhat.nonzero)
   return(ann.nonzero[order(abs(bhat.nonzero), decreasing=T),])
}# end of function




###########################################################################################################
# Main function log.sparse() 
# This function relies on other functions to utilise multi-core
# vectorised computation.
# Use the RhpcBLASctl package(!)


log.sparse <- function(y,X=1,Z, penalty="HL", lambda.seq=exp(seq(-2,6,0.5)), tolerance=1e-4, write=NULL, plot.all=FALSE, opt.crit = "AIC", alpha.mix=0.5, alpha=1.00001, zero.threshold=NULL, fixed.b0=TRUE, epsilon=1e-3){
# Logistic regression with HL random effects
# WITH ELASTIC NET type penalty
# for sparse estimation
# utilising HL likelihood of Lee and Oh (2006)
# input: X, matrix of fixed predictors, default: 1 (fixed intercept)
# Z, matrix of genomic data
# penalty, either "L1", "SCAD", "ridge", "HL", "enet", or "HLnet" (strict)
# (addition) or "HLIG" or "HLIGnet" for inverse gamma (strict)
# lambda.seq, sequence of lambda to be run
# tolerance, convergence criterion
# write, if not NULL, the path to store the value of some key quantities
#      in the iteration -- leave it NULL, unless you know what they mean
# plot.all, logical, shall a plot of criterion be done?
# opt.crit, optimal criterion to estimate lambda, either "AIC" or "BIC" (strict)
# alpha.mix, mixing proportion for HL in HLnet
# (addition) alpha, parameter for HLIG, should be close to one for sparse solution

# Note: Here there is the addition of steps
# in which estimates of random effects less than zero.threshold is set to zero.
# The threshold does not have an effect when penalty="ridge".

# Contact: A.Gusnanto@leeds.ac.uk

    Yfun <- function (X = NULL, y = NULL, Z = NULL, beta = NULL, b = NULL, N = 1, epsilon=1e-3) {
        Z = as.matrix(Z)
        eta = X %*% beta + (Z %*% b)
        p = abs(exp(eta)/(1 + exp(eta))-epsilon)
        wt = (N * p * (1 - p))+(epsilon*0.1)
        Y = eta + ((y - N * p)/(N * p * (1 - p) + epsilon))
        return(list(eta = c(eta), p = c(p), wt = c(wt), Y = c(Y)))
    }# End of Yfun function
    
#solve1 <- function(Z,sinv,dinv){
#d <- 1/dinv
#s <- 1/sinv
#temp <- t(Z)%*%solve(Z%*%diag(d)%*%t(Z)+diag(s))%*%Z
#temp <- d*temp
#temp <- -1*t(t(temp)*d)
#diag(temp) <- d+diag(temp)
#return(temp)
#}

thresh = function(b, penalty, zero.threshold){
if(penalty=="ridge") zero.threshold = 0
if(penalty=="L1" | penalty=="HL" | penalty=="HLIG"){
 if(is.null(zero.threshold)) zero.threshold <- 1e-3}
if(penalty=="enet" | penalty=="HLnet" | penalty=="HLIGnet"){
 if(is.null(zero.threshold)) zero.threshold <- 1e-2}
b[abs(b) < zero.threshold] <- 0
return(b) 
}


    n=dim(Z)[1]  # number of samples
    q=dim(Z)[2]  # number of random effects
    if(length(c(X))==1) X=rep(X,n)
    X = as.matrix(X)
    nc=dim(X)[2]
    beta0=rep(mean(y),nc)  # fixed effects starting values
    #b0 = solve1(Z, sinv=rep(1,n), dinv=rep(median(lambda.seq), q)) %*% t(Z)%*%(y-X%*%beta0)
    # random effects starting values
    a<-3.7 # for SCAD
    w<-30  # for HL, sparse estimates
    b.matrix = NULL # matrix for random effects for different lambdas
    beta.matrix = NULL # matrix for fixed effects for different lambdas
                       # should be the same across columns.
                       # It is put here as a check.
    se.b.matrix = NULL # matrix for standard error of random effects
                       # for different lambdas
    se.beta.matrix = NULL # matrix for standard error of fixed effects
                          #for different lambdas
                          # should be different across columns (a function of lambda)
                          # It is put here as a check.
    AIC.vector = NULL # vector of AIC for each lambda
    BIC.vector = NULL # vector of BIC for each lambda
    df.vector = NULL # vector of df fit for Z (for each lambda)
    loglik.vector = NULL # vector of log likelihood (for each lambda)

    for(lamx in lambda.seq){ # Start iteration for different lambda
    tole=1e10
    iter=1
    beta= beta0
    if(fixed.b0){
     #b0 = solve1(Z, sinv=rep(1,n), dinv=rep(1, q)) %*% t(Z)%*%(y-X%*%beta0)
     b0 = fast_wb_solve(Z, rep(1,q), 1, rep(1,n), as.vector(y-X%*%beta0))
    } else {
     #b0 = solve1(Z, sinv=rep(1,n), dinv=rep(lamx, q)) %*% t(Z)%*%(y-X%*%beta0)
     b0 = fast_wb_solve(Z, rep(1,q), lamx, rep(1,n), as.vector(y-X%*%beta0))
    }# end if else fixed.b0

    b = b0
    sigx<-sd(b0)/sqrt(2)

    if (lamx==0 & penalty=="SCAD") {lamx<-1e-10}

    while(tole>tolerance & iter<=5000){  # Start iteration for a given lambda
     cat("Penalty:", penalty, "log(lambda)",log(lamx),"iter",iter,"rss", tole, "\n")
     old.b = b 
    if (penalty=="L1") {uux<-abs(b)+1e-08}
    if (penalty=="SCAD") {uux<-(abs(b)+1e-08)/ (as.numeric(abs(b)<=lamx)+as.numeric(abs(b)>lamx)*(a*lamx-abs(b))*as.numeric(a*lamx>abs(b))/((a-1)*lamx))}
    if (penalty=="HL"){    kax <- sqrt(4*b^2/(w*sigx^2)+((2/w)-1)^2) ;    uux <- 0.25*w*(((2/w)-1)+kax)+1e-08 }
    if (penalty=="ridge"){uux <- rep(1,q)}
    if (penalty=="enet"){uux <- 0.5*(abs(b)+1e-08)+0.5}
    if (penalty=="HLnet"){ kax <- sqrt(4*b^2/(w*sigx^2)+((2/w)-1)^2) ; uux <- alpha.mix*(0.25*w*(((2/w)-1)+kax)+1e-08)+(1-alpha.mix) }
    if (penalty=="HLIG"){co <-1; alpha1 <- 2*alpha+3; alpha2 <- 2*co*(alpha-1); uux <- (b^2/(2*sigx^2) + alpha2)/alpha1+1e-08}
    if (penalty=="HLIGnet"){co <-1; alpha1 <- 2*alpha+3; alpha2 <- 2*co*(alpha-1); uux <- alpha.mix*((b^2/(2*sigx^2) + alpha2)/alpha1+1e-08)+(1-alpha.mix)}
    WWx<-as.vector(1/uux)

    res.fun = Yfun(X, y, Z, beta, b, 1, epsilon)
            Y = res.fun$Y

     if(!is.null(write)){
     write.table(t(c(iter, log(lamx), b)), file=paste(write,"-b.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
     write.table(t(c(iter, log(lamx),Y)), file=paste(write,"-Y.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
     write.table(t(c(iter, log(lamx),WWx)), file=paste(write,"-WWx.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
     write.table(t(c(iter, log(lamx),uux)), file=paste(write,"-uux.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
     write.table(t(c(iter, log(lamx),res.fun$wt)), file=paste(write,"-wt.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
     }

     beta= solve(t(X)%*%diag(res.fun$wt)%*%X, t(X)%*%diag(res.fun$wt)%*%(Y-Z%*%b))
     vmat <- diag_solve_ZSZD(Z, WWx, lamx, res.fun$wt)
     #vmat <- solve1(Z, sinv=res.fun$wt, dinv=lamx*WWx)
     #b= vmat %*% t(Z*res.fun$wt)%*%(Y-X%*%beta)
     b= fast_wb_solve(Z, WWx, lamx, res.fun$wt, as.vector(Y-X%*%beta))
     b = thresh(c(b), penalty, zero.threshold)
     tole = sum((b-old.b)^2)
     iter = iter+1
     sigx <- sd(b)/sqrt(2)
     gc()
    } # end iteration for one given lambda

     b.matrix = cbind(b.matrix, c(b))  # random effects collection
     beta.matrix = cbind(beta.matrix, c(beta)) # fixed effects collection
     se.b.matrix = cbind(se.b.matrix, c(sqrt(vmat))) # se for b
     #Vtemp <- t(t(Z)*c(1/(lamx*WWx)))%*%t(Z)+diag(1/c(res.fun$wt)) # ZD^{-1}Z'
     Vtemp <- Z_D_Zt_fast(Z, c(1/(lamx*WWx))) + diag(1/c(res.fun$wt))
     Vmat <- solve(t(X)%*%solve(Vtemp)%*%X) # (X'V^{-1}X)^{-1}
     se.beta.matrix = cbind(se.beta.matrix, c(sqrt(diag(Vmat)))) # se for beta

     # AIC and BIC calculation
     loglik = sum(y*log(res.fun$p)+(1-y)*log(1-res.fun$p))
     #df.mat = vmat%*%t(Z*res.fun$wt)%*%Z
     #df=sum(diag(df.mat))
     df = trace_Ainv_B_fast(Z, WWx, lamx, res.fun$wt)
     AIC.vector = c(AIC.vector, -2*loglik+2*df)
     BIC.vector = c(BIC.vector, -2*loglik+log(n)*df)
     df.vector = c(df.vector, df)
     loglik.vector = c(loglik.vector, loglik)
     
    }# end iteration for different lambdas
    
    if(plot.all & opt.crit=="AIC"){
    plot(log(lambda.seq), AIC.vector, type="l", ylab="AIC", 
       xlab=expression(paste("log(",lambda,")")))
    } # end of if(plot.all & opt.crit=="AIC")
    if(plot.all & opt.crit=="BIC"){
    plot(log(lambda.seq), BIC.vector, type="l", ylab="BIC", 
       xlab=expression(paste("log(",lambda,")")))
    } # end of if(plot.all & opt.crit=="BIC")

    # optimal lambda, df, loglik, etc.
    if(opt.crit == "BIC"){
    opt.lambda <- lambda.seq[which.min(BIC.vector)]
    opt.df <- df.vector[which.min(BIC.vector)]
    opt.loglik <- loglik.vector[which.min(BIC.vector)]
    bhat <- b.matrix[,which.min(BIC.vector)]
    se.bhat <- se.b.matrix[,which.min(BIC.vector)]
    betahat <- beta.matrix[,which.min(BIC.vector)]
    se.betahat <- se.beta.matrix[,which.min(BIC.vector)]
     } else {
    opt.lambda <- lambda.seq[which.min(AIC.vector)]
    opt.df <- df.vector[which.min(AIC.vector)]
    opt.loglik <- loglik.vector[which.min(AIC.vector)]
    bhat <- b.matrix[,which.min(AIC.vector)]
    se.bhat <- se.b.matrix[,which.min(AIC.vector)]
    betahat <- beta.matrix[,which.min(AIC.vector)]
    se.betahat <- se.beta.matrix[,which.min(AIC.vector)]
     }
    
     estimates = list(lambda.seq = lambda.seq, log.lambda.seq = log(lambda.seq),
     AIC.vector = AIC.vector, BIC.vector = BIC.vector,
     df.vector = df.vector, loglik.vector = loglik.vector,
     b.matrix = b.matrix, beta.matrix = beta.matrix,
     se.b.matrix = se.b.matrix, se.beta.matrix = se.beta.matrix)
     
    
    result <- list(lambda=opt.lambda, df=opt.df, loglik=opt.loglik,
     bhat=bhat, se.bhat=se.bhat,
     z.bhat = bhat/se.bhat, pval.bhat = pnorm(-abs(bhat/se.bhat))*2,
     lower.ci.bhat = bhat+qnorm(0.025)*se.bhat, upper.ci.bhat = bhat+qnorm(0.975)*se.bhat,
     betahat = betahat, se.betahat = se.betahat, z.betahat = betahat/se.betahat,
     pval.betahat = pnorm(-abs(betahat/se.betahat))*2,
     lower.ci.betahat = betahat+qnorm(0.025)*se.betahat,
     upper.ci.betahat = betahat+qnorm(0.975)*se.betahat,
     estimates = estimates)
     
    return(result)
}# end of function log.sparse



fast_wb_solve <- function(Z, dD, lambda, sS, y,
                          nthreads = NULL,
                          use_chol = TRUE,
                          check_inputs = TRUE) {
  if (check_inputs) {
    if (!is.matrix(Z)) Z <- as.matrix(Z)
    storage.mode(Z) <- "double"
    n <- nrow(Z); p <- ncol(Z)

    stopifnot(length(dD) == p, length(sS) == n, length(y) == n)
    stopifnot(all(is.finite(dD)), all(is.finite(sS)), all(is.finite(y)))
    stopifnot(all(dD > 0), all(sS > 0), is.finite(lambda), lambda > 0)
  } else {
    n <- nrow(Z); p <- ncol(Z)
  }

  # ---- Optional: set BLAS threads via RhpcBLASctl ----
  if (!is.null(nthreads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {

    # RhpcBLASctl exports blas_get_num_procs() and blas_set_num_threads()
    # (there is no exported blas_get_num_threads()).
    old_threads <- tryCatch(
      RhpcBLASctl::blas_get_num_procs(),
      error = function(e) NA_integer_
    )

    tryCatch(
      RhpcBLASctl::blas_set_num_threads(as.integer(nthreads)),
      error = function(e) warning("Could not set BLAS threads: ", conditionMessage(e))
    )

    on.exit({
      if (is.finite(old_threads)) {
        tryCatch(
          RhpcBLASctl::blas_set_num_threads(as.integer(old_threads)),
          error = function(e) NULL
        )
      }
    }, add = TRUE)
  }

  invA <- 1.0 / (lambda * dD)     # length p
  Sy   <- sS * y                  # length n
  b    <- crossprod(Z, Sy)        # length p
  t1   <- invA * as.vector(b)     # length p
  r    <- as.vector(Z %*% t1)     # length n

  # K = S^{-1} + Z diag(invA) Z^T
  Zs <- sweep(Z, 2, invA, `*`)    # n x p
  G  <- tcrossprod(Zs, Z)         # n x n
  K  <- G
  diag(K) <- diag(K) + 1.0 / sS

  # Solve K q = r  (n x n, here n=100)
  if (use_chol) {
    Rchol <- chol(K)
    q <- backsolve(Rchol, forwardsolve(t(Rchol), r))
  } else {
    q <- solve(K, r)
  }

  ztq <- crossprod(Z, q)          # length p
  x   <- t1 - invA * as.vector(ztq)

  as.vector(x)
}



diag_solve_ZSZD <- function(Z, dD, lambda, sS,
                            nthreads = NULL,
                            block_size = 4000L,
                            check_inputs = TRUE,
                            verbose = FALSE) {
  # Returns diag( solve( t(Z)%*%S%*%Z + lambda*D ) )
  # where S=diag(sS), D=diag(dD), both positive diagonal vectors.

  if (check_inputs) {
    if (!is.matrix(Z)) Z <- as.matrix(Z)
    storage.mode(Z) <- "double"

    n <- nrow(Z); p <- ncol(Z)
    stopifnot(length(dD) == p, length(sS) == n)
    stopifnot(is.finite(lambda), lambda > 0)
    stopifnot(all(is.finite(dD)), all(dD > 0))
    stopifnot(all(is.finite(sS)), all(sS > 0))
  } else {
    n <- nrow(Z); p <- ncol(Z)
  }

  # ---- Optional: set BLAS threads via RhpcBLASctl ----
  # RhpcBLASctl exports blas_set_num_threads() and blas_get_num_procs().
  if (!is.null(nthreads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                            error = function(e) NA_integer_)
    tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(nthreads)),
             error = function(e) warning("Could not set BLAS threads: ", conditionMessage(e)))

    on.exit({
      if (is.finite(old_threads)) {
        tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(old_threads)),
                 error = function(e) NULL)
      }
    }, add = TRUE)
  }

  block_size <- as.integer(block_size)
  if (block_size < 1L) block_size <- 1L

  # a_j = 1/(lambda*d_j)  (diag of (lambda D)^-1)
  invA <- 1.0 / (lambda * dD)   # length p

  # ---- Build K = S^{-1} + Z diag(invA) Z^T  in blocks (avoid forming Zs) ----
  K <- diag(1.0 / sS, n, n)

  for (start in seq.int(1L, p, by = block_size)) {
    end <- min(p, start + block_size - 1L)
    idx <- start:end

    Zb <- Z[, idx, drop = FALSE]                       # n x b
    Wb <- sweep(Zb, 2, invA[idx], `*`)                 # n x b (Zb * invA)

    # Add Wb %*% t(Zb) to K; use tcrossprod for BLAS GEMM
    K <- K + tcrossprod(Wb, Zb)                        # n x n

    if (verbose) message("Built K block: ", start, "-", end)
  }

  # ---- Factorize K (n is small: 100) ----
  Rchol <- tryCatch(chol(K),
                    error = function(e) stop("Cholesky failed; K may not be SPD: ", conditionMessage(e)))

  # helper: solve K X = B using chol
  solveK <- function(B) {
    backsolve(Rchol, forwardsolve(t(Rchol), B))
  }

  # ---- Compute diag(Z^T K^{-1} Z) in blocks ----
  diagZtKinvZ <- numeric(p)

  for (start in seq.int(1L, p, by = block_size)) {
    end <- min(p, start + block_size - 1L)
    idx <- start:end

    Zb <- Z[, idx, drop = FALSE]        # n x b
    Qb <- solveK(Zb)                    # n x b  (K^{-1} Zb)

    # For each column j: z_j^T K^{-1} z_j = sum_i Z_ij * Q_ij
    diagZtKinvZ[idx] <- colSums(Zb * Qb)

    if (verbose) message("Diag block: ", start, "-", end)
  }

  # diag(A^{-1}) = invA - invA^2 * diag(Z^T K^{-1} Z)
  inv_diag <- invA - (invA * invA) * diagZtKinvZ
  as.vector(inv_diag)
}


trace_Ainv_B_fast <- function(Z, dD, lambda, sS,
                              nthreads = NULL,
                              block_size = 4000L,
                              check_inputs = TRUE,
                              verbose = FALSE) {
  # Computes: tr( (t(Z) S Z + lambda D)^(-1) * (t(Z) S Z) )
  # where S=diag(sS), D=diag(dD), both positive diagonal vectors.
  #
  # Equivalent to: sum(diag(solve(t(Z)%*%S%*%Z + lambda*D) %*% t(Z)%*%S%*%Z))
  #
  # Returns a single numeric scalar.

  if (check_inputs) {
    if (!is.matrix(Z)) Z <- as.matrix(Z)
    storage.mode(Z) <- "double"
    n <- nrow(Z); p <- ncol(Z)
    stopifnot(length(dD) == p, length(sS) == n)
    stopifnot(is.finite(lambda), lambda > 0)
    stopifnot(all(is.finite(dD)), all(dD > 0))
    stopifnot(all(is.finite(sS)), all(sS > 0))
  } else {
    n <- nrow(Z); p <- ncol(Z)
  }

  block_size <- as.integer(block_size)
  if (block_size < 1L) block_size <- 1L

  # ---- Optional: set BLAS threads via RhpcBLASctl ----
  # RhpcBLASctl exports blas_set_num_threads() and blas_get_num_procs(). 
  if (!is.null(nthreads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                            error = function(e) NA_integer_)
    tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(nthreads)),
             error = function(e) warning("Could not set BLAS threads: ", conditionMessage(e)))
    on.exit({
      if (is.finite(old_threads)) {
        tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(old_threads)),
                 error = function(e) NULL)
      }
    }, add = TRUE)
  }

  # a_j = 1/(lambda*d_j)
  invA <- 1.0 / (lambda * dD)   # length p

  # ---- Build K = S^{-1} + Z diag(invA) Z^T in blocks ----
  K <- diag(1.0 / sS, n, n)

  for (start in seq.int(1L, p, by = block_size)) {
    end <- min(p, start + block_size - 1L)
    idx <- start:end

    Zb <- Z[, idx, drop = FALSE]                 # n x b
    Wb <- sweep(Zb, 2, invA[idx], `*`)           # n x b

    # K += Wb %*% t(Zb)  (BLAS GEMM via tcrossprod)
    K <- K + tcrossprod(Wb, Zb)

    if (verbose) message("K build block: ", start, "-", end)
  }

  # ---- Cholesky of K (n=100, small) ----
  Rchol <- tryCatch(chol(K),
                    error = function(e) stop("Cholesky failed; K may not be SPD: ", conditionMessage(e)))

  solveK <- function(B) {
    backsolve(Rchol, forwardsolve(t(Rchol), B))
  }

  # ---- Accumulate sum_j invA_j * (z_j^T K^{-1} z_j) in blocks ----
  acc <- 0.0

  for (start in seq.int(1L, p, by = block_size)) {
    end <- min(p, start + block_size - 1L)
    idx <- start:end

    Zb <- Z[, idx, drop = FALSE]     # n x b
    Qb <- solveK(Zb)                 # n x b  (K^{-1} Zb)

    # diag(Zb^T K^{-1} Zb): length b
    zKz <- colSums(Zb * Qb)

    acc <- acc + sum(invA[idx] * zKz)

    if (verbose) message("Trace block: ", start, "-", end)
  }

  as.numeric(acc)
}




Z_D_Zt_fast <- function(Z, dD,
                        nthreads = NULL,
                        block_size = NULL,
                        check_inputs = TRUE,
                        verbose = FALSE) {
  # Compute: Z %*% D %*% t(Z) where D = diag(dD)
  # Z: n x p (here n=100, p=60000)
  # dD: length p, diagonal of D (numeric; can be positive/zero/negative, but typical is >=0)
  #
  # If block_size is NULL: one-shot BLAS (allocates a scaled copy of Z).
  # If block_size is integer: blocked computation (less temporary memory).

  if (check_inputs) {
    if (!is.matrix(Z)) Z <- as.matrix(Z)
    storage.mode(Z) <- "double"
    n <- nrow(Z); p <- ncol(Z)
    stopifnot(length(dD) == p)
    stopifnot(all(is.finite(dD)))
  } else {
    n <- nrow(Z); p <- ncol(Z)
  }

  # ---- Optional: set BLAS threads via RhpcBLASctl ----
  # Exported API includes blas_set_num_threads() and blas_get_num_procs(). 
  if (!is.null(nthreads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                            error = function(e) NA_integer_)
    tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(nthreads)),
             error = function(e) warning("Could not set BLAS threads: ", conditionMessage(e)))
    on.exit({
      if (is.finite(old_threads)) {
        tryCatch(RhpcBLASctl::blas_set_num_threads(as.integer(old_threads)),
                 error = function(e) NULL)
      }
    }, add = TRUE)
  }

  # Use sqrt weights: Z D Z^T = (Z * sqrt(d)) (Z * sqrt(d))^T
  # Works for dD >= 0. If dD may be negative, see note below.
  if (any(dD < 0)) {
    stop("dD has negative entries; sqrt(dD) is not real. If D can be negative, use the alternative formula shown in the notes.")
  }
  w <- sqrt(dD)

  # --- Option A: one-shot BLAS (usually fastest) ---
  if (is.null(block_size)) {
    Zs <- sweep(Z, 2, w, `*`)     # n x p (scaled columns)
    return(tcrossprod(Zs))        # n x n via BLAS GEMM
  }

  # --- Option B: blocked accumulation (memory-efficient) ---
  block_size <- as.integer(block_size)
  if (block_size < 1L) block_size <- 1L

  out <- matrix(0.0, n, n)
  for (start in seq.int(1L, p, by = block_size)) {
    end <- min(p, start + block_size - 1L)
    idx <- start:end

    Zb  <- Z[, idx, drop = FALSE]
    Zbs <- sweep(Zb, 2, w[idx], `*`)
    out <- out + tcrossprod(Zbs)

    if (verbose) message("Block ", start, "-", end, " done")
  }
  out
}




###########################################################################################################




cv.log.sparse <- function(y,X=1,Z, fold=5, penalty="HL", lambda.seq=exp(seq(-2,6,0.5)), tolerance=1e-4, write=NULL, plot.all=FALSE, seed=NULL,  opt.crit = "AIC", alpha.mix=0.5, alpha=1.00001, zero.threshold=NULL, fixed.b0=TRUE, epsilon=1e-3){
    # CROSS VALIDATION Logistic regression with HL random effects
    # for sparse estimation based on Lee and Oh (2006)
    # input: X, matrix of fixed predictors, default: 1 (fixed intercept)
    # Z, matrix of genomic data
    # penalty, either "L1", "SCAD", "ridge",  "HL", "HLnet", "HLIG", or "HLIGnet" (strict)
    # lambda.seq, sequence of lambda to be run
    # tolerance, convergence criterion
    # write, if not NULL, the path to store the value of some key quantities
    #      in the iteration -- leave it NULL, unless you know what they mean
    # plot.all, logical, shall a plot of criterion be done?
    # seed, if NULL, it will not be set (let the computer decides)
    # cv.log.sparse uses log.sparse for cross validation
    
    # Contact: A.Gusnanto@leeds.ac.uk
    
    y.pred.fun <- function (X = NULL, Z = NULL, beta = NULL, b = NULL) {
        Z = as.matrix(Z)
        beta = as.matrix(beta)
        eta = X %*% beta + (Z %*% b)
        p = abs(exp(eta)/(1 + exp(eta))-0.000001)
        return(round(p))
    }# End of Yfun function
    p.pred.fun <- function (X = NULL, Z = NULL, beta = NULL, b = NULL) {
        Z = as.matrix(Z)
        beta = as.matrix(beta)
        eta = X %*% beta + (Z %*% b)
        p = abs(exp(eta)/(1 + exp(eta))-0.000001)
        return(p)
    }# End of Yfun function
    
    n=dim(Z)[1]  # number of samples
    q=dim(Z)[2]  # number of random effects
    if(length(c(X))==1) X=rep(X,n)
    X = as.matrix(X)
    nc=dim(X)[2]
    
    # folding groups
    if(!is.null(seed)) set.seed(seed)
    n.y1 <- sum(y==1)
    n.y0 <- sum(y==0)
    fold1 <- sample(c(1:n.y1)%%fold)+1
    fold0 <- sample(c(1:n.y0)%%fold)+1
    fold.group <- vector()
    fold.group[y==1] <- fold1
    fold.group[y==0] <- fold0
    
    # vectors across different lambdas
    pred.error = NULL # matrix of prediction error across
    # different folds (row) and different lambdas (column)
    y.pred.matrix = NULL # matrix of y.pred across different folds (row)
    # and different lambdas (column)
    p.pred.matrix = NULL # matrix of p.pred across different folds (row)
    # and different lambdas (column)
    y.val.vector = NULL # vector of y.val across folds
    
    for(k in 1:fold){
        Z.val <- Z[fold.group==k,]
        Z.est <- Z[fold.group!=k,]
        y.val <- y[fold.group==k]
        y.est <- y[fold.group!=k]
        X.val <- X[fold.group==k,]
        X.est <- X[fold.group!=k,]
        
        y.val.vector <- c(y.val.vector, y.val) # to be compared to y.pred.matrix
        
        res.fold <- log.sparse(y=y.est,X=X.est,Z=Z.est, penalty=penalty, lambda.seq=lambda.seq, tolerance=tolerance, write=write, plot.all=plot.all,  opt.crit =opt.crit, alpha.mix=alpha.mix, alpha=alpha, zero.threshold=zero.threshold, fixed.b0=fixed.b0, epsilon=epsilon)
        
        pred.error.fold <- NULL
        y.pred.fold <- NULL
        p.pred.fold <- NULL
        
        for(kk in 1:length(lambda.seq)){
            y.pred <- y.pred.fun(X = X.val, Z = Z.val, beta = res.fold$estimates$beta.matrix[,kk],
                                 b = res.fold$estimates$b.matrix[,kk])
            p.pred <- p.pred.fun(X = X.val, Z = Z.val, beta = res.fold$estimates$beta.matrix[,kk],
                                 b = res.fold$estimates$b.matrix[,kk])
            pred.error.fold <- c(pred.error.fold, sum(y.val!=y.pred))
            y.pred.fold <- cbind(y.pred.fold, y.pred)
            p.pred.fold <- cbind(p.pred.fold, p.pred)
        } # End iteration across lambda
        pred.error <- rbind(pred.error, pred.error.fold)
        y.pred.matrix <- rbind(y.pred.matrix, y.pred.fold)
        p.pred.matrix <- rbind(p.pred.matrix, p.pred.fold)
    } # End of iteration across fold.
    
    mis.error <- apply(pred.error,2,sum)
    opt.error <- mis.error[which.min(mis.error)]
    
    #optimal lambda
    opt.lambda <- lambda.seq[which.min(mis.error)]
    
    estimates <- list(pred.error=pred.error, y.pred.fold = y.pred.matrix,
                      y.val.fold = y.val.vector, p.pred.fold=p.pred.matrix)
    result <- list(lambda=opt.lambda, lambda.seq=lambda.seq,
                   log.lambda.seq=log(lambda.seq), error=opt.error,
                   error.seq=mis.error, estimates=estimates)
    return(result)
}# end of function




## functions for simulation studies

library(DNAcopy)
library(CNAseg)
library(wavethresh)


seg = function(test.noise, denoise){
    n = ncol(test.noise)
    n.sampl = nrow(test.noise)
    #segmentation
    test.seg = matrix(0,n.sampl,n)
    if(denoise=="TGUHm"){
        require(CNAseg)
        require(wavethresh)
        for (i in 1:n.sampl) {
            test.seg[i,] = tguhm(test.noise[i,], chr=rep(1,n))$segmented
        }
    }
    if(denoise=="CBS"){
        require(DNAcopy)
        for (i in 1:n.sampl) {
            test.seg[i,] = CBS(test.noise[i,], chr=rep(1,n))
        }
    }
    
    return(test.seg)
}


#function for CBS segmentation
CBS = function(obj, chr = rep(1, length(obj))){
    CNA.obj = CNA(obj, chrom = chr, maploc = 1:length(obj), data.type = "binary")
    CBS = segment(CNA.obj)
    num.seg = length(CBS$segRows$startRow)
    CBS.seg = vector()
    for(k in 1:num.seg){
        CBS.seg[CBS$segRows$startRow[k]:CBS$segRows$endRow[k]] = CBS$output$seg.mean[k]
    }
    return(CBS.seg)
}

ndwt <- function(data, type="detail"){
    n = ncol(data)
    n.sampl = nrow(data)
    
    coef.ndwt =vector("list", log2(n))
    for (i in 1:log2(n)) {
        coef.ndwt[[i]] = matrix(NA, nrow=n.sampl, ncol=n)
    }
    
    scale=1
    for(k in (log2(n)-1):0){
        for(i in 1:n.sampl){
            if(type=="detail"){coef.ndwt[[scale]][i,] = temp = accessD(wd(data[i,], filter.number=1, family="DaubExPhase",type="station"), level=k)}
            if(type=="scaling"){coef.ndwt[[scale]][i,] = temp = accessC(wd(data[i,], filter.number=1, family="DaubExPhase",type="station"), level=k)}
            
            coef.ndwt[[scale]][i,(1+(2^(scale-1)-1)):n] = temp[1:(n-(2^(scale-1)-1))]
            coef.ndwt[[scale]][i,1:(2^(scale-1))] = temp[(n-(2^(scale-1)-1)):n]
        }
        scale=scale+1
    }
    
    return(coef.ndwt)
} # end of ndwt()


gen.cna <- function(n.sim=1, n.obs=100, p=1000, n.block=200, true.mu=NULL, effect.diff=0, seed=NULL, true.rho=0.9, block.cor=0, segment="CBS"){
    # function to generate CNA dataset(s)
    # Contact A.Gusnanto@leeds.ac.uk
    require(MASS)
    
    size.block <-  p/n.block
    if(p%%n.block != 0) stop("n.block needs to be common divisor of p.")
    
    if(is.null(true.mu)){
        true.mu <- c()
        cna.level <- c(0.5, 1, 1.5, 1)
        for(j in 1:n.block){
            true.mu <- c(true.mu, rep(cna.level[c(j%%4)+1], size.block))
        }# end for j n.block
        #true.mu <- rep(c(rep(1, size.block), rep(1.5, size.block), rep(1, size.block), rep(0.5, size.block)), p/(size.block*4))
    } # end if is.null true.mu
    if(length(true.mu)==1) true.mu <- rep(true.mu, p)
    
    true.d1 <- c(rep(effect.diff,size.block), rep(-effect.diff, size.block), rep(0,p-2*size.block))
    
    true.Sigma <- matrix(0,p,p)
    temp <- matrix(true.rho, size.block, size.block)
    for(k in 1:n.block){
        true.Sigma[((k-1)*size.block+1):(k*size.block), ((k-1)*size.block+1):(k*size.block)] <- temp
    }# end over number of block
    temp <- matrix(block.cor, size.block, size.block)
    for(k in 1:(n.block-1)){
        true.Sigma[(k*size.block+1):((k+1)*size.block), ((k-1)*size.block+1):(k*size.block)] <- temp
        true.Sigma[((k-1)*size.block+1):(k*size.block), (k*size.block+1):((k+1)*size.block)] <- temp
    }# end over number of block
    diag(true.Sigma) <- 1
    
    sim.data=list()
    for(j in 1:n.sim){
        if(!is.null(seed)) set.seed(seed+j-1)
        Z.sim <- mvrnorm(n.obs, mu=true.mu, Sigma=true.Sigma)
        Z.sim[1:(n.obs/2),] = t(t(Z.sim[1:(n.obs/2),])+true.d1)
        sample.CBS = seg(Z.sim, denoise = segment)
        attr(sample.CBS,"raw") <- Z.sim
        sim.data <- c(sim.data, list(sample.CBS))
    } # end of for j.
    
    sim.setting <- list(n.sim=n.sim, n.obs=n.obs, p=p, n.block=n.block,
                        true.mu=true.mu, effect.diff=effect.diff,
                        seed=seed, true.rho=true.rho, block.cor=block.cor,
                        segment=segment, size.block=size.block)
    attr(sim.data, "setting") <- sim.setting
    return(sim.data)
}# end of gen.cna()



require(pROC)
sim.cv.log.sparse <- function(sim.data, effect.size=1, type.b=1, seed=NULL, fold=5, penalty="HL", 
                              lambda.seq=exp(seq(-3,4,0.5)), tolerance=1e-4, 
                              write=NULL, plot.all=FALSE, opt.crit = "AIC",
                              alpha.mix=0.5, alpha=1.00001, zero.threshold=NULL,
                              fixed.b0=TRUE, epsilon=1e-3, starting=1, ending=NULL, saving=TRUE){
    # function to perform cross validation on the simulated data sim.data from gen.cna() 
    # A.Gusnanto@leeds.ac.uk
    
    get.res <- function(a,b,d){
        a$error <- c(a$error, b$error)
        a$lambda <- c(a$lambda, b$lambda)
        a$bhat <- cbind(a$bhat, d$bhat)
        return(a)
    }#end function
    
    res = list(error=NULL, lambda=NULL, bhat=NULL)
    nm <-deparse(substitute(sim.data))  
    sim.setting <- attr(sim.data,"setting")
    size.block <- sim.setting$size.block
    p <- sim.setting$p
    
    if(length(effect.size)==1 | length(effect.size)==p){
        if(length(effect.size)==1){
            if(type.b==1){
                true.b <- c(c(effect.size, rep(0,size.block-1)), c(-effect.size, rep(0,size.block-1)),
                            rep(0,p-2*size.block))
            } else {
                true.b <- c(rep(effect.size,size.block), rep(-effect.size, size.block),
                            rep(0,p-2*size.block))
            } # end if else type.b==1
        } # end if length true.b==1
    } else {
        stop("Length of true.b is not 1 or p.")
    }
    
    if(is.null(ending)) ending <- sim.setting$n.sim
    
    for(j in starting:ending){
        cat(nm, "j=", j, format(Sys.time()), "\n")
        eta.sim <- sim.data[[j]]%*%true.b
        y.sim <- round(exp(eta.sim)/(1+exp(eta.sim)))
        
        #temp <- log.sparse(y=y.sim,X=1,Z=sim.data[[j]], 
        #                   penalty=penalty, lambda.seq=lambda.seq, 
        #                   tolerance=tolerance, write=write, 
        #                   plot.all=plot.all,   
        #                   opt.crit = opt.crit, alpha.mix=alpha.mix, 
        #                   alpha=alpha, zero.threshold=zero.threshold,
        #                   fixed.b0=fixed.b0, epsilon=epsilon)
        temp <- log.sparse(y.sim,1,sim.data[[j]], penalty, lambda.seq, 
                           tolerance, write, plot.all, opt.crit, alpha.mix, alpha, 
                           zero.threshold, fixed.b0, epsilon)
        
        #temp.cv <-   cv.log.sparse(y=y.sim,X=1,Z=sim.data[[j]], fold=fold,
        #                           penalty=penalty, lambda.seq=temp$lambda, 
        #                           tolerance=tolerance, write=write, 
        #                           plot.all=plot.all, seed=seed,  
        #                           opt.crit = opt.crit, alpha.mix=alpha.mix, 
        #                           alpha=alpha, zero.threshold=zero.threshold,
        #                           fixed.b0=fixed.b0, epsilon=epsilon)
        temp.cv <- cv.log.sparse(y.sim,1,sim.data[[j]], fold, 
                                 penalty, temp$lambda, tolerance, write, 
                                 plot.all, seed,  opt.crit, 
                                 alpha.mix, alpha, zero.threshold, fixed.b0, epsilon)
            
        res <- get.res(res, temp.cv, temp)
        temp.auc <- auc(roc(c(temp.cv$estimates$y.val.fold), c(temp.cv$estimates$p.pred.fold)))
        gc()
            if(saving){
                write.table(temp.auc, file=paste(nm,"-", penalty,"-auc.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
                write.table(t(temp$bhat), file=paste(nm,"-", penalty,"-b.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
                write.table(temp.cv$error, file=paste(nm,"-", penalty,"-error.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
                write.table(temp.cv$lambda, file=paste(nm,"-", penalty,"-lambda.txt", sep=""), row.names=F, col.names=F, quote=F, append=T)
            }
            
    }# end for each simulated data
    
    return(res)
    
} #end of sim.cv.log.sparse()



