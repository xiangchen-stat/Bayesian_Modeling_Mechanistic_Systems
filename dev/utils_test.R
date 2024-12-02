compare_R_cpp <- function(res1, res2){
  out <- list()
  for (i in 1:length(res1)) {
    ele1 <- res1[[i]]
    ele2 <- res2[[i]]
    if(is.null(dim(ele1))){
      out[[i]] <- ele1 - ele2
    }else{
      out[[i]] <- norm(ele1 - ele2, "F")
    }
  }
  return(out)
}