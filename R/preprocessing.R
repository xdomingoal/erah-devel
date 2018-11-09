
pre.process <- function(sampl, sigma.scans)
{
	sigma.model.filter <- round(sigma.scans)
	if(is.even(sigma.model.filter)) sigma.model.filter <- sigma.model.filter + 1 
	filt.order <- 3
	if(filt.order>=sigma.model.filter) filt.order <- sigma.model.filter - 1
	if(filt.order>=1) sampl <- soft.filter(sampl,filt.order,sigma.model.filter)
	sampl[sampl<0] <- 0
	
	k.filt <- trunc(sigma.scans*10)
	moving.maximas <- apply(sampl,2,function(x){max(x)})
	sampl[,moving.maximas<100] <- 0
	sampl <- apply(sampl,2,function(x) removeBaseline(x,k.filt))

	sampl[sampl<0] <- 0
	sampl
}

#' @importFrom stats sd

removeBaseline <- function (x, k) 
{
    if (max(x) == 0) return(x)
    x.min <- runningmin(x, k)
    k.m <- k
    if(k.m>=length(x)) k.m <- length(x) - 1
    if (is.even(k.m) == T) k.m <- k.m - 1
    x.dupCopy <- c(x[length(x):1],x,x[length(x):1])
    x.med <- stats::runmed(x.dupCopy, k.m, endrule='keep')[(length(x)+1):(2*length(x))]

    base.sd <- sd(x.min)
    x.base <- x.med
    x.base[x.base > (x.min + base.sd)] <- x.min[x.base > (x.min + base.sd)] + base.sd
    x.base <- runningmean(x.base, k)
    x.clean <- x - x.base
    x.clean[x.clean < 0] <- 0
    x.clean
}

#' @importFrom signal sgolayfilt

soft.filter <- function (mat, p, n = NULL) 
{
    if (is.null(n)) 
        n <- p + 3
        
    mat.f <- apply(mat, 2, function(x) sgolayfilt(x, p = p , n = n))
    mat.D <- normalize(mat.f)
    mat.D <- sweep(mat.D,2,apply(mat,2, max),"*")
    mat.D
}