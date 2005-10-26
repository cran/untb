"expected.abundance" <- function(J, theta){
  a <- parts(J)
  f <- function(x){theta.prob(theta=theta,x=extant(count(x)),give.log=FALSE)}
  probs <- apply(a,2,f)
  return(count(apply(sweep(a,2,probs,"*"),1,sum)))
}


"display.untb" <- function(start, gens=100, prob.of.mutate=0, cex=3,
                      individually=TRUE, ask=FALSE,
                      flash=FALSE, delay=0, cols=NULL, ...){

  spp <- as.integer(as.census(start))
  n <- length(spp)
  if(is.null(cols)){
    cols <- rgb(red=runif(n),blue=runif(n),green=runif(n))
  }
  n.sqrt <- floor(sqrt(n))
  
  e <- expand.grid(1:n.sqrt , 1:n.sqrt)
  colnames(e) <- c("" , "")
  plot(e,col=cols[(spp %% n)+1], pch=16, cex=cex,
       axes=FALSE, ...)

  if(individually){
    for(i in 1:gens){
      j <- sample(n,1)
      if(runif(1)<prob.of.mutate){
        spp[j] <- max(spp)+1
      } else {
        spp[j] <- sample(spp,1)
      }
      if(ask){
        readline("hit return to continue")
      }
      if(flash){
        for(i in 1:8){
          points(e[j,],pch=1,cex=4,col="black")
          Sys.sleep(0.15)
          points(e[j,],pch=1,cex=4,col="white")
          Sys.sleep(0.15)
        }
      }
      points(e[j,],col=cols[spp[j]],pch=16,cex=cex, ...)
      Sys.sleep(delay)

    }
  } else {
    for(i in 1:gens){
      plot(e,col=cols[(spp %%n)+1],
           pch=16,cex=cex,axes=FALSE, ...)
      spp <- select(spp,prob.of.mutate=prob.of.mutate)
            if(ask){
        readline("hit return to continue")
      }

      Sys.sleep(delay)
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

"no.of.spp" <- function(x,include.extinct=FALSE){
  if(include.extinct){
    return(length(as.count(x)))
  } else {
    return(sum(as.count(x) > 0))
  }
}

"no.of.ind" <- function(x){
  sum(as.count(x))
}

"no.of.extinct" <- function(x){
  sum(as.count(x)==0)
}

"species.table" <- function(x){
  if(is.vector(x)){x <- t(as.matrix(x))}
  drop(t(apply(x,1,tabulate,nbins=max(x))))
}

"plot.census" <- function(x, uncertainty=FALSE, expectation=FALSE, theta=NULL, n=10, ...){
  x <- as.count(x)
  plot.default(x,log="y", col="red",pch=16,type="b", xlab="species rank in abundance", ylab="abundance", ...)
  if(uncertainty|expectation){
    J <- no.of.ind(x)
    if(is.null(theta)){
      theta <- optimal.theta(x)
      }
  }
  if(uncertainty){
    for(i in 1:n){
      jj <- rand.neutral(J=J, theta=theta)
      points(1:no.of.spp(jj),jj,type="l",col="gray")
    }
  }
  if(expectation){
    points(1:J,expected.abundance(J=J,theta=theta),type="b",pch=16)
  }
} 

"plot.count" <- plot.census

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

"optimal.theta" <- function(x, interval=NULL, N=NULL, like=NULL, ...){
  x <- as.count(x)
  J <- no.of.ind(x)
  S <- no.of.spp(x)
  if(is.null(interval)){
    interval <- c(0.001/J,J)
  }
  jj <- 
    optimize(f=theta.likelihood, interval=interval, maximum=TRUE,
             give.log=TRUE, x=NULL, S=S,J=J, ...)
  
  theta.mle <- jj$maximum
  max.like <- jj$objective

  if( !is.null(N) & !is.null(like)){
    stop("N and like non-null.   Specify one only")
  } 
  if(!is.null(N)){ #run parametric resampling
    theta.dist <- rep(NA,N)
    for(i in 1:N){
      jj <- rand.neutral(J=J, theta=theta.mle)
      theta.dist[i] <- Recall(x=jj, ...)
    }
    return(theta.dist)
  }
  if(!is.null(like)){ #run likelihood estimate for credible interval
    g <- function(theta){
      theta.likelihood(theta=theta, S=S,J=J, give.log=TRUE)-max.like+like
    }
    return(c(
             lower=uniroot(f=g,lower=1/J,upper=theta.mle)$root,
             mle=theta.mle,
             upper=uniroot(f=g,lower=theta.mle,upper=J)$root
             )
           )
  }
  return(theta.mle)  # return MLE.
}

"optimal.prob" <- function(x, interval=NULL, N=NULL, like=NULL, ...){
  optimal.theta(x=x,interval=interval, N=N, like=like, ...)/(2*no.of.ind(x))
}


#"species.abundance" <- function(x){
#  xx <- unique(x)
#  out <- rbind(species.id=xx, abundance=tabulate(match(x,xx)))
#  out <- out[,order(out[2,],decreasing=TRUE),drop=FALSE]
#  colnames(out) <- rep(" ",ncol(out))
#  out}


"theta.prob" <-
  function(theta, x=NULL, give.log=TRUE){
    J <- no.of.ind(x)
    S <- no.of.spp(x)
    jj <- phi(x)
    out <-  (
             + theta.likelihood(theta=theta,S=S,J=J,give.log=TRUE)
             + lgamma(J+1)
             - sum(jj*log(1:length(jj)))
             - sum(lgamma(jj+1))
             )
    
    if(give.log){
      return(out)
    }  else {
      return(exp(out))
    }
  }


"theta.likelihood" <- function(theta, x=NULL, S=NULL, J=NULL, give.log=TRUE){
   if(!is.null(x)){
    J <- no.of.ind(x)
    S <- no.of.spp(x)
  }
  if(give.log){
    return(S*log(theta) + lgamma(theta) - lgamma(theta+J))
  } else {
    return(theta^S*exp(lgamma(theta)-lgamma(theta+J)))
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
  if(is.census(a)){
    return(a)
  } else {
    return(census(as.count(a)))
  }
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
  object <- as.count(object) 
  out <- list(
           no.of.ind         = no.of.ind(object),
           no.of.spp         = no.of.spp(object),
           no.of.singletons  = no.of.singletons(object),
           maximal.abundance = object[1],
           most.abundant.spp = names(object[1])
           )
  class(out) <- "summary.count"
  return(out)
}

"summary.census" <- summary.count

"print.summary.count" <- function(x, ...){
  cat("Number of individuals:", x[[1]],"\n") 
  cat("Number of species:", x[[2]],"\n") 
  cat("Number of singletons:", x[[3]],"\n") 
  cat("Most abundant species: ", x[[5]]," (",x[[4]]," individuals)","\n",sep="")
}

"simpson" <- function(x){
  x <- as.count(x)
  J <- no.of.ind(x)
  return(1-sum(x*(x-1))/(J*(J-1)))
}

"rand.neutral" <- function(J, theta=NULL, prob.of.mutate=NULL, string=NULL, pad=FALSE){
  if(!xor(is.null(theta) , is.null(prob.of.mutate))){
    stop("must supply exactly one of theta and prob.of.mutate")
  }
  if(is.null(theta)){
    theta <- 2*J*prob.of.mutate
  }
  out <- rep(NA,J)
  spp <- 1
  out[1] <- spp
  for(j in 2:J){
    sg <- theta/(theta+j-1)
    if(runif(1) < sg){
      #new species
      spp <- spp+1
      out[j] <- spp
    } else {
      #existing species
      out[j] <- sample(out[1:(j-1)],size=1)
    }
  }
  if(!is.null(string)){
    out <- paste(string,out,sep="")
  }
  out <- as.count(out)
  if(pad){
    extras <- J-no.of.spp(out)
    if(extras>0){
      jj <- rep(0,extras)
      names(jj) <- paste("extinct.",1:extras,sep="")
      out <- out + count(jj)
    }
  }
  return(out)
}

"extractor" <- function(x, index){
  jj <- names(x)
  x <- x[index,,drop=FALSE]
  if(isTRUE(all.equal(1,nrow(x)))){
        x <- as.numeric(x)
    names(x) <- jj
    return(count(x))
  } else {
    return(x)
  }
}

"extant" <- function(x){
  x <- as.count(x)
  count(x[x>0])
}

"extinct" <- function(x){
  x <- as.count(x)
  count(x[x==0])
}
  
"+.count" <- function(a,b){
  a <- as.count(a)
  b <- as.count(b)
  both <- c(a,b)
  as.count(as.table(tapply(both,names(both),sum)))
}

"logkda.a11" <- function(a){
  data(logS1)
  N <- no.of.ind(a)
  S <- no.of.spp(a)

  "f" <- function(x){
    total <- 0
    for(i in 1:length(x)){
      total <- total + logS1[a[i],x[i]]
    }
    total <- total + sum(lgamma(x)) - sum(lgamma(a))
    return(exp(total))
  }
  
  kda <- c(1,rep(NA,N-S))
  for(A in (S+1):N){
    jj <- 1+blockparts(A-S,a-1)
    kda[A-S+1] <- sum(apply(jj,2,f))
  }
  return(log(kda))
}
  
"logkda" <- function(a){
  a <- as.count(a)
  i <- 1
  maxabund <- max(a)
  jj <- rev(a)
  specabund <- rbind(unique(jj),table(jj))
  if(specabund[1,i]==1){
    i <- i+1
  }
  polyn <- 1

  Told <- 1:i
  for(n in 2:maxabund){
    Tnew <- (n > (1:n))*Told[pmin( (n-1),1:n)] + Told[pmax(1,(1:n)-1)]*( ((1:n)-1))/(n-1)
#    print(Tnew)
    if(specabund[1,i] == n){
      for(k0 in 1:specabund[2,i]){
        lenpolyn2 <- length(polyn) + length(Tnew)-1
        newpolyn <- rep(NA,lenpolyn2)
        for(k1 in 1:lenpolyn2){
          k2 <- max(1,k1+1-length(Tnew)):min(length(polyn),k1)
          newpolyn[k1] <- sum(polyn[k2]*Tnew[k1+1-k2])
        }
        polyn <- newpolyn
#        print("polyn starts")
#        print(polyn)
#        print("polyn ends")
      }
      i <- i+1
    }
    Told <- Tnew[1:n]
#    print("Told=",Told)
  }
  log(polyn)
}



"etienne" <- function(theta, m, D, log.kda=NULL, give.log=TRUE, give.like=TRUE){
  J <- no.of.ind(D)
  S <- no.of.spp(D)

  
  if(is.null(log.kda)){
    log.kda <- logkda(D)
  }
  A <- S:J
  if(m != 1){
    I <- m*(J-1)/(1-m)
    correction.factor <- 
      sum(  
          exp(
              log.kda
              +lgamma(theta+J)
              -lgamma(theta+A)
              +lgamma(I)
              -lgamma(I+J)
              +A*log(I)
              )
          )
  } else {
    correction.factor <- 1
  }
  if(give.like){
    if(give.log){
      return(theta.likelihood(theta=theta,S=S,J=J, give.log=TRUE) + log(correction.factor))
    } else {
      return(theta.likelihood(theta=theta,S=S,J=J, give.log=FALSE) * correction.factor)
    }
  } else {
    jj <- theta.prob(theta=theta,x=D) * correction.factor
    if(give.log){
      return(log(jj))
    } else {
      return(jj)
    }
  }
}

"optimal.params" <- function(D, start=NULL, give=FALSE, ...){
  log.kda <- logkda(D)
  if(is.null(start)){
    thetadash <- log(optimal.theta(D))
    mdash <- 0
  } else {
    thetadash <- log(start[1])
    mdash <- tan(  (pi/2)*(2*start[2]-1)  )
  }

  par <- c(thetadash=thetadash, mdash=mdash)  
  f <- function(p){
    thetadash=p[1]
    mdash=p[2]
    jj <- 
    -etienne(theta = exp(thetadash),
            m = 0.5*(1+(2/pi)*atan(mdash)),
            D = D, log.kda = log.kda, give.like = TRUE,give.log=TRUE)
   
    return(jj)
  }
  
  jj <- optim(par=par, fn=f, ...)
  if(give){
    return(jj)
  } else {
    jj <- jj$par
    names(jj) <- NULL
    return(c(
             theta = exp(jj[1]),
             m = 0.5*(1+(2/pi)*atan(jj[2]))
             ))
  }
}
