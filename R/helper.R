match.likelihood <- function(hh, Y, lmt, dict, freq, covlkl=NULL, mmdir=c("both", "patient", "donor")){
  q <- 1
  n <- dim(Y)[1]

  mmdir <- match.arg(mmdir)
  if(!is.null(covlkl) && (!is.numeric(covlkl) | covlkl<=0 | covlkl>1)){
    stop("Covered likelihood threshold must be NULL or be a value between ]0, 1]")
  }

  if (n==0){return(0)}
  for(i in 1:n){
    pi <- 0
    if(is.null(covlkl)){
      y <- expandDonor(Y[i,], dict, freq)
    } else {
      y <- expandMostLk(Y[i,], dict, freq, thr=covlkl)
    }

    if(any(apply(y$tmat, 1, function(temp){as.character(temp)==hh}))){
      w <- apply(y$tmat, 1, function(temp){yj<-as.character(temp);MMcomp(hh,yj)[mmdir]<=lmt})*1
      pi <- sum(y$p*w)

    }
    q <- q*(1-pi)
  }
  return(1-q)
}


MMcomp <- function(p, d, type=c("both", "GvH", "HvG")){
  if(length(p)!=length(d)){stop("not same number of loci!")}
  type <- match.arg(type, several.ok=T)
  mp <- pairsum(!p %in% d & !is.na(p))
  md <- pairsum(!d %in% p & !is.na(d))
  mm <- pairmax(mp, md)
  res <- c(sum(mm), sum(mp), sum(md))
  names(res) <- c("both", "GvH", "HvG")
  return(res[type])
}


pairsum <- function(a){
  nloc <- length(a)/2
  b <- NULL
  for(i in 1:nloc){
    b <- c(b,sum(a[c(2*i-1,2*i)]))
  }
  return(b)
}


pairmax <- function(p, d){
  if(length(p)!=length(d)){stop("not same number of loci!")}
  n <- length(p)
  mm <- sapply(1:n, function(a){max(p[a], d[a])})
  return(mm)
}


MMpart <- function(p, d, type=c("both", "patient", "donor")){
  if(length(p)!=2*length(d)){stop("not same number of loci!")}
  type <- match.arg(type, several.ok=T)
  mp <- pairsum(!p %in% d)
  md <- !d %in% p
  mm <- pairmax(mp, md)
  res <- c(sum(mm), sum(mp), sum(md))
  names(res) <- c("both", "patient", "donor")
  return(res[type])
}


expandDonor <- function(X, dict, freq, CRA=c("mean", "CRA"), simplify=T){
  CRA=match.arg(CRA)

  if(CRA=="mean" | is.na(X$CRA)){
    P <- lapply(freq, function(x){apply(x, 1, mean)})
  } else if(CRA=="CRA"){
    P <- lapply(freq, function(x){tmp <- x[,X$CRA]; names(tmp) <- rownames(x); tmp})
  }
  temp <- P[[1]]
  for (i in 2:length(P)){
    temp <- c(temp, P[[i]])
  }
  P <- temp; rm(temp)

  tlist <- list()
  for (l in setdiff(colnames(X), c("IND", "CRA"))){
    if (grepl("\\d{2,4}:\\d{2,4}|\\d{2,3}:N$", X[,l]) | is.na(X[,l])){
      tlist[l] <- X[, l]
    } else {
      tlist[l] <- dict[X[, l]]
    }
  }
  tlist <- sapply(names(tlist), function(x){if(is.null(tlist[[x]])){
    grep(paste0("^", gsub("_\\d", "", x)), names(P), value=T)}else{tlist[[x]]
    }})
  tmat <- expand.grid(tlist)
  if(dim(tmat)[2]==1){
    tmat <- as.data.frame(t(tmat))
    colnames(tmat) <- names(tlist)
    rownames(tmat) <- "1"
  }

  corr.all <- NULL
  pmat <- apply(tmat, c(1,2), function(x){P[x]})
  pmat[is.na(pmat)] <- 0
  null.all <- as.character(unique(unlist(tmat[,apply(pmat, 2, function(k){all(k==0)})])))
  if(length(null.all)==0){
    null.all <- NULL
  }
  pmat[,apply(pmat, 2, function(k){all(k==0)})] <- 0.000001
  p <- apply(pmat, 1, prod)

  if(simplify){
    tmat <- tmat[p!=0,]
    p <- p[p!=0]
  }

  return(list(IND=X$IND, tmat=tmat, p=p/sum(p), lkl=sum(p), CRA=X$CRA, null.all=null.all))
}


expandMaxLk <- function(X, dict, freq, CRA=c("mean", "CRA")){
  CRA <- match.arg(CRA)
  y <- expandDonor(X, dict, freq, CRA)
  wm <- which.max(y$p)
  y$tmat <- y$tmat[wm,]
  y$lkl <- y$lkl*y$p[wm]
  y$p <- y$p[wm]
  return(y)
}

###########################################
#'@export
###########################################
expandMostLk <- function(X, dict, freq, CRA=c("mean", "CRA"), thr=0.95, grp=NULL){
  if(!is.numeric(thr) | thr<=0 | thr>1){stop("Likelihood threshold must be a value between ]0, 1]")}
  CRA <- match.arg(CRA)
  y <- expandDonor(X, dict, freq, CRA)

  if(!is.null(grp)){
    y$tmat <- apply(y$tmat, c(1,2), function(x){if(x %in% unlist(grp) & !grepl("\\d{2,4}:\\d{2,4}N", x)){names(grp)[which(sapply(grp, function(y){x %in% y}))]}else{x}})
    temp <- data.frame(y$tmat, "p"=y$p)
    temp <- as.data.frame(summarise_at(group_by(temp[,c(colnames(y$tmat),'p')], temp[,colnames(y$tmat)]), vars(p), ~ sum(.,na.rm=TRUE)))
    y$tmat <- temp[,colnames(y$tmat)]
    y$p <- temp[,"p"]

  }

  yp <- sort(y$p, T)
  ym <- yp[min(which(cumsum(yp)>=thr))]
  wm <- which(y$p>=ym)
  y$tmat <- y$tmat[wm,]
  y$lkl <- y$lkl*sum(y$p[wm])
  y$p <- y$p[wm]
  y$p_rel <- sum(y$p)
  y$p <- y$p/sum(y$p)
  return(y)
}

###########################################
#'@export
###########################################
compileMostLk <- function(Y, dict, freq, CRA=c("mean", "CRA"), thr=0.95, grp=NULL){
  require(foreach)
  temp <- foreach(x=1:dim(Y)[1], .combine=rbind) %do% {y <- expandMostLk(Y[x,], dict, freq, CRA=CRA, thr=thr, grp=grp); cbind("IND"=Y[x,"IND"], y$tmat,"p"=y$p)}
  y <- list()
  y$IND <- temp$IND
  y$tmat <- apply(temp[,setdiff(colnames(temp), c("IND","p"))], c(1,2), as.character)
  y$p <- temp$p
  return(y)
}
