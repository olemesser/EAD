rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- ceiling(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }

  if (pos.only) while (any(vec <= 0)) {
    negs <- vec <= 0
    pos  <- vec >= 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

setupMC<-function(NUMB_CORES=1,logfile=""){
  library(parallel)
  library(doSNOW)
  cores_available<-parallel::detectCores(logical = T)
  if(cores_available<NUMB_CORES){
    NUMB_CORES<-cores_available
  }
  cl <- snow::makeSOCKcluster(NUMB_CORES,outfile=logfile)
  return(cl)
}


par_apply<-function(cl, X, FUN,export=NULL,packages=loadedNamespaces(),...){
  pb <- txtProgressBar(max=length(X),style=3)
  on.exit(close(pb))
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  out<-foreach(param=X,  .options.snow=opts,
               .packages = packages,
               .export = export,
               .errorhandling = 'remove') %dopar% {
                  temp<-tryCatch({withTimeout(FUN(param), timeout=time_limit)},
                           error=function(e){
                             message("Time limit reached!")
                             return(list())
                             })
            return(temp)
          }
  return(out)
}

makeMatrixsymmetric<-function(A){
  A_upper <- A
  A_lower <- A
  A_upper[lower.tri(A_upper)]<-t(A)[lower.tri(A)]
  A_lower[upper.tri(A_lower)]<-t(A)[upper.tri(A)]
  A<-A_upper+A_lower
  A[A>1]<-1
  return(A)
}


with_timeout <- function(expr, timeout=10){
  expr <- substitute(expr)
  envir <- parent.frame()
  setTimeLimit(cpu = timeout, elapsed = timeout, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  res<-tryCatch(eval(expr, envir = envir),
           error=function(e) NA)

  out<-list(res=res,
            message=ifelse(is.na(res),"error","success")[1])
  return(out)
}


try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf)
{
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE)
  if(inherits(y, "try-error")) NULL else y
}
