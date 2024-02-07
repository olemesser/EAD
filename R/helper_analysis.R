data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      qt25 = as.numeric(quantile(x[[col]],0.25)),
      qt75 =  as.numeric(quantile(x[[col]],0.75)),
      IQR = IQR(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

cor_table<-function(res,filename,group=NULL){
  library(xlsx)
  require(dplyr)
  require(tidyr)
  require(data.table)
  require(tibble)

  filename<-paste0(Sys.Date(),"_",filename,".xlsx")

  CI<-res$ci %>%
    tibble::rownames_to_column("name") %>%
    mutate(src=sapply(strsplit(name,"-"),function(x) x[1]),
           tgt=sapply(strsplit(name,"-"),function(x) x[2]),
           CI=paste0("[",round(lower,2),",",round(upper,2),"]"),
           count=2,
           id=row_number(),
           method=as.character(res$Call)[3]) %>%
    tidyr::uncount(count) %>%
    group_by(id) %>%
    mutate(id_cnt=cumsum(id==first(id))) %>%
    ungroup() %>%
    mutate(src_temp=ifelse(id_cnt==1,src,tgt),
           tgt_temp=ifelse(id_cnt==1,tgt,src)) %>%
    select(src_temp,tgt_temp,CI,method) %>%
    data.table::as.data.table() %>%
    data.table::dcast("method + src_temp  ~ tgt_temp",value.var="CI") %>%
    as_tibble()

  if(!is.null(group)){
    CI<-CI %>%
      left_join(group,by=c("src_temp"="var")) %>%
      relocate(grp,.after = "src_temp") %>%
      arrange(method,grp)
    CI<-CI[,c(1,2,3,match(CI$src_temp,colnames(CI)))]
    order<-match(CI$src_temp,rownames(res$stars))
    res$stars <-res$stars[order,]
    res$stars <- res$stars[,order]
    order<-match(CI$src_temp,rownames(res$r))
    res$r <-res$r[order,]
    res$r <-res$r[,order]
    res$p <-res$p[order,]
    res$p <-res$p[,order]
  }else{
    CI<-CI[match(rownames(res$stars),CI$src_temp),]
    CI<-CI[,c(1,2,match(CI$src_temp,colnames(CI)))]
  }




  coef<-res$stars %>%
    as_tibble() %>%
    mutate(method=as.character(res$Call)[3],.before=1)

  coef_numb<-res$r %>%
    as_tibble() %>%
    mutate_all(round,digits=3) %>%
    mutate(method=as.character(res$Call)[3],.before=1)

  coef_pval<-res$p %>%
    as_tibble() %>%
    mutate_all(round,digits=3) %>%
    mutate(method=as.character(res$Call)[3],.before=1)


  write.xlsx(coef,
             file = filename,sheetName = "coef")
  write.xlsx(coef_numb,
              file = filename,sheetName = "coef_numb",append = T)
  write.xlsx(coef_pval,
             file = filename,sheetName = "pval",append = T)
  write.xlsx( CI,
              file = filename,sheetName = "CI",append = T)

  return(list(res=res,
              CI=CI))

}


# sem_effects <- function(fit){
#   require(lavaan)
#
#   beta <- lavInspect(fit, "est")$beta
#   y <- matrix(0,nrow = NROW(beta),ncol = NCOL(beta))
#   n <- 1
#   repeat{
#     effect <- eval(parse(text = paste(rep("beta",n),collapse = " %*% ")))
#     y <- y + effect
#     if(sum(effect)==0){
#       message(paste0("Effects of order ",n-1," were calculated."))
#       break
#     }else{
#       n <- n + 1
#     }
#   }
#   return(y)
# }


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}



matrix_plot_pre <- function(A){
  require(dplyr)
  require(data.table)
  colnames(A)<-1:NCOL(A)
  A<-A %>%
    as_tibble() %>%
    mutate(X=1:n()) %>%
    as.data.table() %>%
    melt.data.table(id.vars = "X",variable.name = "Y") %>%
    mutate_all(as.numeric) %>%
    mutate(color = case_when(
      value == 0 ~ 'white',
      value >0 ~ 'black'
    )) %>%
    mutate(color = ifelse(X == Y,'grey50',color))
  return(A)
}
