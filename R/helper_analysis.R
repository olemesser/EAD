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
