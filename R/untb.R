"display.untb" <- function(n=400, gens=900000, prob.of.mutate=0.001, cex=3,
                      individually=TRUE, start.mono=TRUE,
                      flash=FALSE,t1=0, flashsleep=0.1, ...){

  cc <- rgb(red=runif(n),blue=runif(n),green=runif(n))
  n.sqrt <- floor(sqrt(n))
  
  e <- expand.grid(1:n.sqrt , 1:n.sqrt)
  if(start.mono){
    spp <- rep(1,n)
  } else {
    spp <- 1:n
  }
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

"select" <- function(a,prob.of.mutate=0){
  a <- sample(a,replace=TRUE)
  mutate <- runif(length(a))<prob.of.mutate
  if(any(mutate)){
    a[mutate] <- max(a) +(1:sum(mutate))
  }
  return(a)
}
  
"untb" <- function(n=100, prob.of.mutate=0.01, individually=FALSE,
  gens=150, keepall=FALSE, start.mono=FALSE){
  if(start.mono){
    a <- rep(1,n)
  } else {
    a <- 1:n
  }
  if(keepall){aa <- matrix(NA,gens,n)}
  if(individually){
    for(i in 1:gens){
      if(keepall){ aa[i,] <- a}
      j <- sample(n,1)
      if(runif(1) < prob.of.mutate){
        a[j] <- max(a)+1
      } else {
        a[j] <- sample(a,1)
      }
    }
  } else {
      if(keepall){
        for(i in 1:gens){
          aa[i,] <- a
          a <- select(a,prob.of.mutate)
        }
      } else {
        for(i in 1:gens){
          a <- select(a,prob.of.mutate)
        }
      }
    }
  if(keepall){
    return(aa)
  } else {
    return(a)
  }  
}

"normalize" <- function(x,keepnames=FALSE){
  x <- sort(table(x),decreasing=TRUE)
  out <- rep(1:length(x),x)
  if(keepnames){
    names(out) <- rep(names(x),x)
  } else {
    names(out) <- NULL
  }
  return(out)
}

"species.count" <- function(x){
  if(is.vector(x)){x <- t(as.matrix(x))}
  apply(x,1,function(u){length(unique(u))})
}

"species.table" <- function(x){
  if(is.vector(x)){x <- t(as.matrix(x))}
  drop(t(apply(x,1,tabulate,nbins=max(x))))
}

"species.curve" <- function(x, show.uncertainty=FALSE, n=10, ...){
  if(!is.table(x)){
    x <- table(x)
  }
  plot(sort(x,decreasing=TRUE),log="y", col="red",pch=16,type="b")
  if(show.uncertainty){
    prob <- optimal.prob(x)
    for(i in 1:n){
      jj <- untb(n=sum(x), prob.of.mutate=prob, ...)
      jj <- sort(table(jj),decreasing=TRUE)
      points(1:length(jj),jj,type="l",col="gray")
    }
  }
} 

"preston" <- function(x, n=8){
  if(n<2){stop("n must be > 2")}
  breaks <- c(0,2^(0:(n-2)),Inf)
  outnames <- paste("(",breaks[-n-1],",",breaks[-1],"]",sep="")
  if(is.matrix(x)){
    out <- t(apply(x,1,function(u){preston(u)}))
    colnames(out) <- outnames
    return(out)
  }

  if(!is.table(x)){
    x <- table(x)
  }
  breaks[is.infinite(breaks)] <- max(x)
  r <- hist(x,plot=FALSE,breaks=breaks,right=TRUE)
  out <- r$counts
  names(out) <- outnames
  out <- t(as.matrix(out))
  rownames(out) <- "number of species"
  return(out)
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
 
"species.abundance" <- function(x){
  xx <- unique(x)
  out <- rbind(species.id=xx, abundance=tabulate(match(x,xx)))
  out <- out[,order(out[2,],decreasing=TRUE),drop=FALSE]
  colnames(out) <- rep(" ",ncol(out))
  out}

"theta.prob" <- function(theta, x, give.log=FALSE){
  if(is.table(x)){
    J <- sum(x)
    S <- length(x)
    jj <- phi(x)
  } else {
    J <- length(x)
    S <- species.count(x)
    jj <-  tabulate(x)
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
  function(theta,x=NULL,S=species.count(x),J=length(x),give.log=FALSE){
    if(!missing(x)){
      if(is.table(x)){
        J <- sum(x)
        S <- length(x)
      } else {
        J <- length(x)
        S <- species.count(x)
      }
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

"phi" <- function(x){
  if(is.table(x)){
    return(tabulate(x))
  } else {
    return(tabulate(tabulate(normalize(x))))
  }
}

"count.to.census" <- function(a){
rep(names(a),a)
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

fisher.ecosystem <- function(N, S, nmax, alpha=NULL, c=0){
  if(is.null(alpha)){alpha <- fishers.alpha(N,S)}
  x <- N/(N+alpha)
  j <- 1:nmax
  jj <- rpois(n=nmax,lambda=alpha*x^j/(j+c))
  return(rep(1:sum(jj) , rep(j,jj)))
}
