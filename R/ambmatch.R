###########################################
#'@export
###########################################
find.donors <- function(X, Y, dict, freq, match='8/8', verbose=F, title=NULL, covlkl=NULL){
  require(reshape2)
  require(dplyr)
  tm <- Sys.time()
  i=1; blnk=T; while(blnk){L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[i,1], '~')));
                            blnk <- "blank" %in% L; i<-i+1}
  #L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[1,1], '~')))
  X <- cbind(colsplit(X$alleles, '~', L), freq=X[,2])
  X <- X[X$freq!=0,]
  match=unlist(strsplit(match, '/'))
  lmt <- as.integer(match[2]) - as.integer(match[1])
  sel <- switch(match[2], '8'=c('A','B','C','DRB1'), '10'=c('A','B','C','DRB1','DQB1'),
                '12'=c('A','B','C','DRB1','DQB1','DPB1'))
  X <- as.data.frame(summarise_at(group_by(X[,c(sel,'freq')], X[,sel]), vars(freq), ~ sum(.,na.rm=TRUE)))

  A <- sapply(1:dim(Y)[1], function(a){
    as.character(unique(unlist(expandDonor(Y[a,], dict, freq)$tmat)))
  })

  n <- dim(X)[1]
  if(verbose){progressbar <- txtProgressBar(max=(n^2+n)/2, style=3)
  if(!is.null(title)){cat(title)}}

  P <- 0
  for(i in 1:n){
    hi <- X[i,'freq']
    for (j in i:n){
      k <- 1+(i!=j)
      hj <- X[j,'freq']
      hh <- unlist(X[c(i,j), sel])
      Yfilt <- which(sapply(A, function(a){sum(!hh %in% a)<=lmt}))
      if(length(Yfilt)==0){next}
      p <- match.likelihood(hh, Y[Yfilt,], lmt, dict, freq, covlkl)
      P <- P + k*hi*hj*p
      tt <- (i-1)*(n)-ifelse(i==1, 0, factorial(i-1))+j-i+1
      if(verbose){setTxtProgressBar(progressbar, value = tt)}
    }
  }
  if(verbose){close(progressbar); print(Sys.time()-tm)}
  return(P)
}


###########################################
#'@export
###########################################
find.donorsMostLk <- function(X, Y, dict, freq, match='8/8', verbose=F, title=NULL, covlkl=0.95){
  require(reshape2)
  require(dplyr)
  tm <- Sys.time()
  i=1; blnk=T; while(blnk){L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[i,1], '~')));
  blnk <- "blank" %in% L; i<-i+1}
  #L <- sub("^([^\\*]+)\\*[^\\*]+", "\\1", unlist(strsplit(X[1,1], '~')))
  X <- cbind(colsplit(X$alleles, '~', L), freq=X[,2])
  X <- X[X$freq!=0,]
  match=unlist(strsplit(match, '/'))
  lmt <- as.integer(match[2]) - as.integer(match[1])
  sel <- switch(match[2], '8'=c('A','B','C','DRB1'), '10'=c('A','B','C','DRB1','DQB1'),
                '12'=c('A','B','C','DRB1','DQB1','DPB1'))
  X <- as.data.frame(summarise_at(group_by(X[,c(sel,'freq')], X[,sel]), vars(freq), ~ sum(.,na.rm=TRUE)))

  y <- compileMostLk(Y, dict, freq, CRA=c("mean", "CRA"), thr=covlkl)

  n <- dim(X)[1]
  if(verbose){progressbar <- txtProgressBar(max=(n^2+n)/2, style=3)
  if(!is.null(title)){cat(title)}}

  P <- 0
  for(i in 1:n){
    hi <- X[i,'freq']
    for (j in i:n){
      k <- 1+(i!=j)
      hj <- X[j,'freq']
      hh <- unlist(X[c(i,j), sel])

      Yfilt <- which(apply(apply(y$tmat, c(1,2), as.character), 1, MMcomp, p=hh, type="both")<=lmt)
      if(length(Yfilt)==0){next}
      temp <- data.frame("IND"=y$IND[Yfilt], "p"=y$p[Yfilt])
      p <- 1-prod(1-as.data.frame(summarise_at(group_by(temp, temp$IND), vars(p), ~ sum(.,na.rm=TRUE)))$p)

      P <- P + k*hi*hj*p
      tt <- (i-1)*(n)-ifelse(i==1, 0, factorial(i-1))+j-i+1
      if(verbose){setTxtProgressBar(progressbar, value = tt)}
    }
  }
  if(verbose){close(progressbar); print(Sys.time()-tm)}
  return(P)
}
