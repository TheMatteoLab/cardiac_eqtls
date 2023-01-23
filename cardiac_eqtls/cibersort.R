library(e1071)
library(parallel)
library(preprocessCore)
CoreAlg <- function(X, y, absolute, abs_method)
{
    #try different values of nu
    svn_itor <- 3
    res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
    }
    if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)
    t <- 1
    while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
    }
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
    if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]
    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}
doPerm <- function(perm, X, Y, absolute, abs_method)
{
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()
    while(itor <= perm){
        yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
        yr <- (yr - mean(yr)) / sd(yr)
        result <- CoreAlg(X, yr, absolute, abs_method)
        mix_r <- result$mix_r
        if(itor == 1) {dist <- mix_r}
        else {dist <- rbind(dist, mix_r)}
        itor <- itor + 1
    }
    newList <- list("dist" = dist)
}
CIBERSORT <- function(signature_matrix, mixture_matrix, perm=0, QN=TRUE, absolute=FALSE, abs_method='sig.score')
{
    if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score') stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
    X <- data.matrix(signature_matrix)
    Y <- data.matrix(mixture_matrix)
    X <- X[order(rownames(X)),]
    Y <- Y[order(rownames(Y)),]
    P <- perm #number of permutations
    if(QN == TRUE){
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }
    Yorig <- Y
    Ymedian <- max(median(Yorig),1)
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    YintX <- Ygns %in% Xgns
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY,]
    X <- (X - mean(X)) / sd(as.vector(X))
    if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
    header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
    if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999
    while(itor <= mixtures)
	{
        y <- Y[,itor]
        y <- (y - mean(y)) / sd(y)
        result <- CoreAlg(X, y, absolute, abs_method)
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse
        if(absolute && abs_method == 'sig.score') {
            w <- w * median(Y[,itor]) / Ymedian
        }
        if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
        out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
        if(absolute) out <- c(out, sum(w))
        if(itor == 1) {output <- out}
        else {output <- rbind(output, out)}
        itor <- itor + 1
    }
    #write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
    obj <- rbind(header,output)
    obj <- obj[,-1]
    obj <- obj[-1,]
    obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
    rownames(obj) <- colnames(Y)
    if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
    else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
    obj
}
