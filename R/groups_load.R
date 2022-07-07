###########################################
#'@export
###########################################
group_load <- function(grp=c("p", "g"), twofield=T, forceG2f=F, pth=NULL){
  grp <- match.arg(grp)

  if(is.null(pth)){
    pth <- system.file(paste0("extdata/hla_nom_", grp, ".txt"), package = "ProbMatch")
    if(pth==""){stop()}
  }
  print(pth)

  Graw <- read.table(pth, sep = ";", skip = 6)
  Graw <- Graw[Graw[,3]!="",]

  if(twofield){
    Ggrp <- list()
    for (i in 1:dim(Graw)[1]){
      Ggrp[[paste0(Graw[i,c(1,3)], collapse="")]] <-
        paste0(Graw[i,1], sort(unique(paste0(regmatches(unlist(strsplit(Graw[i,2], "/")),
           regexpr("\\d{2,4}:\\d{2,4}", unlist(strsplit(Graw[i,2], "/")))),
              paste0(c("","N")[(regexpr("N$", unlist(strsplit(Graw[7,2], "/")))!=-1)+1],
                     c("","Q")[(regexpr("Q$", unlist(strsplit(Graw[7,2], "/")))!=-1)+1],
                     c("","L")[(regexpr("L$", unlist(strsplit(Graw[7,2], "/")))!=-1)+1])))))
    }
  } else {
    Ggrp <- list()
    for (i in 1:dim(Graw)[1]){
      Ggrp[[paste0(Graw[i,c(1,3)], collapse="")]] <-
        paste0(Graw[i,1], grep("\\d{2,4}:\\d{2,4}", unlist(strsplit(Graw[i,2], "/")), value=T))
    }
  }
  if(grp=="g" & twofield & forceG2f){
    Ggrp <- merge_g_2field(Ggrp)
  }
  return(Ggrp)
}


###########################################
#'@export
###########################################
merge_g_2field <- function(Ggrp){
  mname <- paste0(regmatches(names(Ggrp), regexpr(".{1,4}\\*\\d{2,4}:\\d{2,4}", names(Ggrp))), "G")
  uname <- unique(mname)
  ggrp <- list()
  for (u in uname){
    ggrp[[u]] <- sort(unique(unlist(Ggrp[which(mname==u)])))
  }
  return(ggrp)
}

###########################################
#'@export
###########################################
grab_groups <- function(dest=getwd()){
  download.file("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt", paste0(dest, "/hla_nom_p.txt"))
  download.file("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt", paste0(dest, "/hla_nom_g.txt"))
}
