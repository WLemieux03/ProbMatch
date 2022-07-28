###########################################
#'@export
###########################################
prob.matchOpt <- function(X, N, Y=NULL, avail=1, match="8/8", type=c("GvH", "both", "HvG"), raw=F, ignore.null=T){
  tm <- Sys.time()
  require(reshape2)
  require(dplyr)

  type <- match.arg(type)

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
  nX <- dim(X)[1]
  nY <- dim(Y)[1]

  allk <- sort(unique(c(unlist(X[,sel]), unlist(Y[,sel]))))
  allkid <- 1:length(allk)
  names(allkid) <- allk
  if(ignore.null){
    allkid[grepl("\\d{2,4}:\\d{2,4}N", allk)] <- 0
  }

  hX <- cbind(apply(X[,sel], c(1,2), function(x){allkid[x]}), freq=X$freq)
  hY <- cbind(apply(Y[,sel], c(1,2), function(x){allkid[x]}), freq=Y$freq)

  genoX <- expand.grid(1:nX, 1:nX)
  genoX <- genoX[apply(genoX, 1, function(x){x[1]>=x[2]}),c(2,1)]
  colnames(genoX) <- c("h1", "h2")
  genoX$k <- apply(genoX, 1, function(x){(x[1]!=x[2])+1})
  genoX$hh <- apply(genoX, 1, function(x){paste(sort(unlist(hX[as.integer(c(x[["h1"]], x[["h2"]])), sel])), collapse = "~")})
  genoX$freq <- apply(genoX, 1, function(x){hX[as.integer(x[["h1"]]), "freq"] * hX[as.integer(x[["h2"]]), "freq"]})

  genoY <- expand.grid(1:nY, 1:nY)
  genoY <- genoY[apply(genoY, 1, function(x){x[1]>=x[2]}),c(2,1)]
  colnames(genoY) <- c("h1", "h2")
  genoY$k <- apply(genoY, 1, function(x){(x[1]!=x[2])+1})
  genoY$hh <- apply(genoY, 1, function(x){paste(sort(unlist(hY[as.integer(c(x[["h1"]], x[["h2"]])), sel])), collapse = "~")})
  genoY$freq <- apply(genoY, 1, function(x){hY[as.integer(x[["h1"]]), "freq"] * hY[as.integer(x[["h2"]]), "freq"]})

  P <- apply(genoX, 1, function(x){
    sum(apply(genoY, 1, function(y){
      (MMcomp(unlist(strsplit(x[["hh"]], "~")), unlist(strsplit(y[["hh"]], "~")), type)<=lmt)*as.integer(y[["k"]])*as.numeric(y[["freq"]])
    }))
  })

  genoX$p <- P

  if(raw){
    return(genoX)
  } else {
    P <- sum(apply(genoX, 1, function(x){
      as.integer(x[["k"]])*as.numeric(x[["freq"]])*(1-(1-as.numeric(x[["p"]]))^(avail*N))
    }))
    return(P)
  }
}

###########################################
#'@export
###########################################
recomp.matchOpt <- function(genoX, N, avail=1){
  P <- sum(apply(genoX, 1, function(x){
    as.integer(x[["k"]])*as.numeric(x[["freq"]])*(1-(1-as.numeric(x[["p"]]))^(avail*N))
  }))
  return(P)
}

###########################################
#'@export
###########################################
max.matchOpt <- function(genoX){
  P <- sum(apply(genoX, 1, function(x){
    as.integer(x[["k"]])*as.numeric(x[["freq"]])*(as.numeric(x[["p"]])!=0)
  }))
  return(P)
}
