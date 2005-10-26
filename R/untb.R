"display.untb" <- function(start, gens=100, prob.of.mutate=0.001, cex=3,
                      individually=TRUE, 
                      flash=FALSE,t1=0, flashsleep=0.1, ...){

  spp <- as.integer(as.census(start))
  n <- length(spp)
  cc <- rgb(red=runif(n),blue=runif(n),green=runif(n))
  n.sqrt <- floor(sqrt(n))
  
  e <- expand.grid(1:n.sqrt , 1:n.sqrt)

  plot(e,col=cc[(spp %% n)+1
  ],pch=16,cex=cex,axes=FALSE,xlab="",ylab="", ...)

  if(individually){
    for(i in 1:gens){
      j <- sample(n,1)
      if(runif(1)<prob.of.mutate){
        spp[j] <- max(spp)+1
      } else {
        spp[j] <- sample(spp,1)
      }
      if(flash){
        for(i in 1:8){
          points(e[j,],pch=1,cex=4,col="black")
          Sys.sleep(flashsleep)
          points(e[j,],pch=1,cex=4,col="white")
          Sys.sleep(flashsleep)
        }
      }
      points(e[j,],col=cc[spp[j]],pch=16,cex=cex, ...)
      Sys.sleep(t1)
    }
  } else {
    for(i in 1:gens){
      plot(e,col=cc[(spp %%n)+1
    ],pch=16,cex=cex,axes=FALSE,xlab="",ylab="", ...)
      spp <- select(spp,prob.of.mutate=prob.of.mutate)
      Sys.sleep(t1)
    }
  }
  return(invisible(spp))
}

"select" <- function(a, D=length(a), prob.of.mutate=0, meta=NULL){
  n <- length(a)
  died <- sample(n,D,replace=TRUE)
  mutated <- runif(length(died))<prob.of.mutate
  n1 <- sum(mutated)
  n2 <- sum(!mutated)

  a[died[!mutated]] <- sample(a,n2,replace=TRUE)
  if(is.null(meta)){
    a[died[mutated]] <- (1:n1) + max(a)
  } else {
    a[died[mutated]] <- sample(meta,n1,replace=TRUE)
  }
  return(a)
}
  
"untb" <- function(start, prob.of.mutate=0, D=1,
  gens=150, keep=FALSE, meta=NULL){
  if(is.census(start)|is.count(start)){
    a <- as.integer(as.census(start))
  } else {
    a <- start
  }
  n <- length(a)
  if(keep){
    aa <- matrix(NA,gens,n)
    for(i in 1:gens){
      aa[i,] <- a
      a <- select(a,D=D,prob.of.mutate=prob.of.mutate,meta=meta)
    }
    return(aa)
  } else {
    for(i in 1:gens){
      a <- select(a,D=D,prob.of.mutate=prob.of.mutate,meta=meta)
    }
    return(as.count(a))
  }
}

"species.count" <- function(x){
  if(is.vector(x)){x <- t(as.matrix(x))}
  apply(x,1,function(u){length(unique(u))})
}

"no.of.spp" <- function(x){
  length(as.count(x))
}

"no.of.ind" <- function(x){
  sum(as.count(x))
}

"species.table" <- function(x){
  if(is.vector(x)){x <- t(as.matrix(x))}
  drop(t(apply(x,1,tabulate,nbins=max(x))))
}

"abundance.curve" <- function(x, show.uncertainty=FALSE, n=10, ...){
  x <- as.count(x)
  plot(sort(x,decreasing=TRUE),log="y", col="red",pch=16,type="b")
  if(show.uncertainty){
    prob <- optimal.prob(x)
    for(i in 1:n){
      jj <- untb(start=x, prob.of.mutate=prob, ...)
      points(1:length(jj),jj,type="l",col="gray")
    }
  }
} 

"preston" <- function(x, n=8,original=FALSE){
  if(n<2){stop("n must be >= 2")}
  breaks <- c(0,2^(0:(n-2)),Inf)

  if(is.matrix(x)){
    out <- t(apply(x,1,function(u){preston(u,n=n)}))
    colnames(out) <- outnames
    return(out)
  }
  x <- as.count(x)  
  breaks[is.infinite(breaks)] <- max(x)
  r <- hist(x,plot=FALSE,breaks=breaks,right=TRUE)
  out <- r$counts
  if(original){
    transfer <- phi(x)[2^(0:(n-1))]/2
    transfer[is.na(transfer)] <- 0
    out <- out-transfer
    out[-1] <- out[-1]+transfer[-n]
    outnames <- 2^(0:(n-1))
  } else {
    breaks[length(breaks)] <- Inf
    if(n>2){
      outnames <- c("1 ", "2",paste(" ",breaks[-c(1:2,length(breaks))]+1,"-",breaks[-c(1:3)],sep=""))
    } else {
      outnames <- c("1 ","2-Inf")
    }
  }
  names(out) <- outnames
  class(out) <- "preston"
  return(out)
}




"print.preston" <- function(x, ...){
  x <- t(as.matrix(x))
  rownames(x) <- "number of species"
  class(x) <- "matrix"
  NextMethod("print")
}

"optimal.prob" <- function(x, interval=NULL, ...){
  if(!is.table(x)){
    x <- table(x)
  }
  J <- sum(x)
  if(is.null(interval)){
    interval <- c(0.001/J,J)
  }

  theta <-
optimize(f=theta.likelihood,interval=interval,maximum=TRUE,give.log=TRUE,x=x, ...)$maximum 
  prob <- theta/(2*J)
  return(prob)
}
 
#"species.abundance" <- function(x){
#  xx <- unique(x)
#  out <- rbind(species.id=xx, abundance=tabulate(match(x,xx)))
#  out <- out[,order(out[2,],decreasing=TRUE),drop=FALSE]
#  colnames(out) <- rep(" ",ncol(out))
#  out}

"theta.prob" <-
  function(theta, x=NULL, S=no.of.spp(x), J=no.of.ind(x), give.log=FALSE){
    if(!missing(x)){
      J <- no.of.ind(x)
      S <- no.of.spp(x)
    }
    if(give.log){
      mylog <- function(x){ifelse(x>0,log(x),0)}
      
      return(
             lgamma(J+1) + S*log(theta) - sum(jj*mylog(1:length(jj))) -
             sum(lgamma(jj+1)) - sum(log( (1:J)+theta-1))
             )
    }  else {
      return(
             (factorial(J)*theta^S)/
             (  prod((1:length(jj))^jj)*prod(factorial(jj))*prod((1:J)+theta-1))
             )
    }
  }

"theta.likelihood" <-
  function(theta, x=NULL, S=no.of.spp(x), J=no.of.ind(x), give.log=FALSE){
    if(!missing(x)){
      J <- no.of.ind(x)
      S <- no.of.spp(x)
    }
    if(give.log){
      return(
             (S-1)*log(theta) - sum(log((1:J) + theta - 1))
             )
    } else {
      return(
             theta^(S-1)/prod((2:J)+theta-1)
             )     
    }
  }

"phi" <- function(x,addnames=TRUE){
  x <- as.count(x)
  jj <- table(x)
  out <- tabulate(x)
  if(addnames){
    names(out)[out==1] <- names(x[x %in% as.numeric(names(jj[jj==1]))])
    names(out)[out != 1] <- ""
  }
  return(out) 
}

"is.census" <- function(a){
  inherits(a,"census")
}

"is.count" <- function(a){
  inherits(a,"count")
}

"census" <- function(a){
  out <- as.factor(rep(names(a),times=a))
  class(out) <- c("census","factor")
  return(out)
}

"as.census" <- function(a){
  census(as.count(a))
}

"count" <- function(a){
  out <- sort(as.table(a),decreasing=TRUE)
  class(out) <-  c("count","table")
  return(out)
}

"as.count" <- function(a,add=""){
  if(is.count(a)){
    out <- a
  } else if(is.data.frame(a)) {
    if(nrow(a) != 1){
      stop("data frame supplied: must have only one row")
    }
    out <- as.numeric(a)
    names(out) <- colnames(a)
    out <- count(out[out>0])
  } else if(is.table(a)) {
    out <- count(a)
  } else {
    out <- count(table(a))
  }    
  names(out) <- paste(add,names(out),sep="")
  return(out)
}

fishers.alpha <- function(N, S, give=FALSE){
#  print("gives Fisher's table 9, p55")
  f <- function(a){S+a*log((a/(a+N)))}
  a <- uniroot(f,interval=c(1/N,N))$root
  x <- N/(N+a)
  if(give){
    return(list(a=a,x=x))
  } else {
    return(a)
  }
}

"fisher.ecosystem" <- function(N, S, nmax, alpha=NULL, c=0){
  if(is.null(alpha)){alpha <- fishers.alpha(N,S)}
  x <- N/(N+alpha)
  j <- 1:nmax
  jj <- rpois(n=nmax,lambda=alpha*x^j/(j+c))
  return(rep(1:sum(jj) , rep(j,jj)))
}

"singletons" <- function(x){
  x <- as.count(x)
  names(x[x==1])
}

"no.of.singletons" <- function(x){
  x <- as.count(x)
  return(sum(x==1))
}

"summary.count" <- function(object, ...){
  cat("Number of individuals:", no.of.ind(object),"\n") 
  cat("Number of species:", no.of.spp(object),"\n") 
  cat("Number of singletons:", no.of.singletons(object),"\n") 
  cat("Most abundant species: ", names(object[1])," (",object[1],")","\n",sep="")
}

"summary.census" <- function(object, ...){
  summary(as.count(object))
}

"simpson" <- function(x){
  x <- as.count(x)
  J <- no.of.ind(x)
  return(1-sum(x*(x-1))/(J*(J-1)))
}
