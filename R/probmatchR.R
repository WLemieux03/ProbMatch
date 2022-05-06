require(reshape2)
require(dplyr)
require(foreach)

cross.find <- function(h1, N, h2=NULL, col1=NULL, col2=NULL,  avail=1, match='10/10', verbose=F){
  if(is.null(h2)){
    h2 <- h1
  }
  if (length(N)==1){
    N <- rep(N, length(colnames(h2)[-1]))
  }
  if (!is.null(col2) & length(col2)==1){
    K <- foreach(cra=colnames(h2)[-1], .combine = cbind) %do% {
      find.donors(h1[,c('alleles',col)], N[cra], h2[,c('alleles',cra)], avail=avail, match=match, verbose = verbose)
    }
    colnames(K) <- colnames(h2)[-1]
  } else if ((is.null(col2) & length(colnames(h2)[-1])>1) | length(col2)!=1 ) {
    if (is.null(col2)){
      col2 <- colnames(h2)[-1]
    }
    if (is.null(col1)){
      col1 <- colnames(h1)[-1]
    }
    K <- foreach(CRA=col1, .combine = rbind) %:% foreach(cra=col2, .combine = cbind) %do% {
      prob.matchR(h1[,c('alleles',CRA)], N[cra], h2[,c('alleles',cra)], avail=avail, match=match, verbose = verbose, title=paste0(CRA,"<-",cra))
    }
    rownames(K) <- col1
    colnames(K) <- col2
  }
  return(K)
}

#X=h6l_PN[,c('alleles',CRA)]; N=1000;Y= h6l[,c('alleles',cra)]; avail=0.7
prob.matchR <- function(X, N, Y=NULL, avail=1, match='10/10', verbose=F, title=NULL){
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
  n <- dim(X)[1]
  if(verbose){progressbar <- txtProgressBar(max=(n^2+n)/2, style=3)
              if(!is.null(title)){cat(title)}}

  P <- 0
  for(i in 1:n){
    hi <- X[i,'freq']
    for (j in i:n){
      if(verbose){setTxtProgressBar(progressbar, value = i*(2*n-1-i)/2+j)}
      k <- 1+(i!=j)
      hj <- X[j,'freq']
      hh <- unlist(X[c(i,j), sel])
      p <- donorlkR(hh, Y[apply(Y[,sel], 1, function(x){sum(is.na(match(x, hh)))})<=lmt,], sel = sel, lmt = lmt)
      P <- P + k*hi*hj*(1-(1-p)^(avail*N))
    }
  }
  if(verbose){close(progressbar)}
  return(P)
}

# A, B, C, DR, DQ. (DP > TCE)

donorlkR <- function(hh, Y, sel, lmt){
  p <- 0
  n <- dim(Y)[1]
  if (n==0){return(0)}
  for(i in 1:n){
    #if(){next}
    hi <- Y[i,'freq']
    for (j in i:n){
      k <- 1+(i!=j)
      hj <- Y[j,'freq']
      w <- match(unlist(Y[c(i),sel]), hh)
      w <- c(w,match(unlist(Y[c(j),sel]), hh[-w[!is.na(w)]]))
      w <- sum(is.na(w))<=lmt
      p <- p + k*hi*hj*w
    }
  }

  return(p)
}


#EG CRA=12 & 8/8; 40:2.6%, 200:11.6%, 1000:37.3%, 5000:67.2%
#EG CRA=12 & 10/10; 1000:34.5%

#EG CRA=01 ds 01, avail=1; 1223:42.78%
#EG CRA=01 ds 02, avail=1; 1724:20.48%
#EG CRA=01 ds 03, avail=1; 4678:28.96%
#EG CRA=01 ds 04, avail=1; 1693:
#EG CRA=01 ds 05, avail=1; 2062:
#EG CRA=01 ds 06, avail=1; 12789:35.68%
#EG CRA=01 ds 07, avail=1; 2480:
#EG CRA=01 ds 08, avail=1; 914:17.57%
#EG CRA=01 ds 12, avail=1; 2674:
#EG CRA=01 ds 13, avail=1; 2735:
#EG CRA=01 ds 14, avail=1; 3240:
#EG CRA=01 ds 15, avail=1; 3903:
#EG CRA=01 ds 16, avail=1; 9906:
#EG CRA=01 ds 17, avail=1; 1555:
