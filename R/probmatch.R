require(reshape2)
require(dplyr)

prob.match <- function(X, N, Y=NULL, avail=1, match='10/10', title=NULL){
  i=1; blnk=T; while(blnk){L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[i,1], '~'))); blnk <- "blank" %in% L; i<-i+1}
  #L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[1,1], '~')))
  X <- cbind(colsplit(X$alleles, '~', L), freq=X[,2])
  X <- X[X$freq!=0,]
  match=unlist(strsplit(match, '/'))
  lmt <- as.integer(match[2]) - as.integer(match[1])
  sel <- switch(match[2], '8'=c('A','B','C','DRB1'), '10'=c('A','B','C','DRB1','DQB1'), '12'=c('A','B','C','DRB1','DQB1','DPB1'))
  X <- as.data.frame(summarise_at(group_by(X[,c(sel,'freq')], X[,sel]), vars(freq), ~ sum(.,na.rm=TRUE)))
  if(is.null(Y)){
    Y <- X
  } else {
    Y <- cbind(colsplit(Y$alleles, '~', L), freq=Y[,2])
    Y <- Y[Y$freq!=0,]
    Y <- as.data.frame(summarise_at(group_by(Y[,c(sel,'freq')], Y[,sel]), vars(freq), ~ sum(.,na.rm=TRUE)))
  }
  nx <- dim(X)[1]
  ny <- dim(Y)[1]
  if (dim(X)[2]!=dim(Y)[2]){
    stop()
  }
  m <- dim(X)[2]-1

  fact <- as.factor(unique(c(unlist(X[,sel]),unlist(Y[,sel]))))
  matx <- as.matrix(apply(X[,sel], c(1,2), function(x){as.numeric(factor(x,fact))}))
  maty <- as.matrix(apply(Y[,sel], c(1,2), function(y){as.numeric(factor(y,fact))}))

  P <- donorlk(matx, X$freq, maty, Y$freq, nx, ny, length(sel), lmt, avail, N)
  return(P)
}
