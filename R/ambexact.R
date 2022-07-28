###########################################
#'@export
###########################################
exact.donorsMostLk <- function(X, Y, dict, freqX, freqY, match='8/8', type=c("GvH", "both", "HvG"), CRA=c("mean", "CRA"), thrX=0.95, thrY=0.95, grp=NULL){
  tm <- Sys.time()
  require(reshape2)
  require(dplyr)
  match.arg(type)

  match=unlist(strsplit(match, '/'))
  lmt <- as.integer(match[2]) - as.integer(match[1])
  sel <- switch(match[2], '8'=c('A','B','C','DRB1'), '10'=c('A','B','C','DRB1','DQB1'),
                '12'=c('A','B','C','DRB1','DQB1','DPB1'))

  x <- compileMostLk(X, dict, freqX, CRA=CRA, thr=thrX, grp=grp)

  y <- compileMostLk(Y, dict, freqY, CRA=CRA, thr=thrY, grp=grp)

  n <- dim(x$tmat)[1]

  P <- apply(x$tmat, 1, function(hh){
    isn <- grepl("\\d{2,4}:\\d{2,4}N", hh)
    hh[isn] <- NA

    Yfilt <- which(apply(apply(y$tmat, c(1,2), as.character), 1, MMcomp, p=hh, type="GvH")<=lmt)
    if(length(Yfilt)==0){return(0)}
    temp <- data.frame("IND"=y$IND[Yfilt], "p"=y$p[Yfilt])
    p <- 1-prod(1-as.data.frame(summarise_at(group_by(temp, temp$IND), vars(p), ~ sum(.,na.rm=TRUE)))$p)
    return(p)
  })

  res <- data.frame(IND=x$IND, p=P*x$p)
  res <- summarise_at(group_by(res, res$IND), vars(p), ~ sum(.,na.rm=TRUE))

  return(sum(res$p)/length(res$p))
}
