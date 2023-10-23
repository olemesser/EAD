run_MC<-function(cl,
                 X,
                 FUN,
                 packages=NULL,
                 errorhandlingNodes = 'remove',
                 export = NULL,
                 progress=T,
                 ...){
  suppressMessages(require(parallel))
  suppressMessages(require(doSNOW))
  suppressMessages(require(foreach))

  registerDoSNOW(cl)
  on.exit(stopCluster(cl))

  #### Install libs on remote workers ####
  sucess_lib_install <- install_libs_on_nodes(cl,libs = packages)
  print(sucess_lib_install)


  #### Create Progress Bar ###
  if(progress){
    pb <- txtProgressBar(max=length(X),style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    on.exit(close(pb))
  }else{
    opts <- NULL
  }


  #### Run Function in Multicore mode ####
  out<-foreach(input=X,
               .options.snow=opts,
               .packages = packages,
               .errorhandling = errorhandlingNodes,
               .export = export
  ) %dopar% {
    FUN(input,...)
  }

  return(out)
}


install_libs_on_nodes<-function(cl,libs=NULL){
  require(parallel)
  libs <- unique(c(libs, c("doSNOW", "parallel", "foreach")))
  localhost <- Sys.info()["nodename"]
  if (!is.null(libs)) {
    lib_path <- .libPaths(include.site = F)
    parallel::clusterExport(cl, envir = environment(), c("libs",
                                                         "lib_path", "localhost"))
    host_names<-sapply(cl, function(x) x$host)
    sucess <- clusterEvalQ(cl[match(unique(host_names),host_names)], {
      sucess <- sapply(libs, function(x) {
        if (Sys.info()["nodename"] != localhost) {
          if (!is.element(x, .packages(all.available = TRUE))) {
            .libPaths(c(lib_path))
            install.packages(x)
            x %in% .packages(all.available = TRUE)
            return(x %in% .packages(all.available = TRUE))
          }
          else {
            update.packages(oldPkgs = x)
            return(TRUE)
          }
        }else {
          return(TRUE)
        }
      })
      if (all(!is.null(unlist(sucess)))) {
        message(paste0("Installation on ", Sys.info()["nodename"],
                       " sucessfull for: ", sum(sucess), "/", length(libs)))
        if(any(unlist(sucess) == FALSE)){
          print(sucess)
          message("Not all packages were sucessfull installed on ",Sys.info()["nodename"])
        }
      }
      else {
        message("Nothing to do on ",Sys.info()["nodename"])
      }
      return(sucess)
    })
  }else {
    sucess <- NULL
  }
  return(sucess)
}
