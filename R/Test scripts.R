require(tidyr)
simA <- crossing(c("A*01", "A*02", "A*03", "A*11"),
         c("B*07", "B*13", "B*51"),
         c("C*01", "C*02", "C*03"),
         c("DRB1*01", "DRB1*02", "DRB1*04"))
simH <- apply(simA, 1, paste, collapse="~")
set.seed(20220503); simF <- runif(length(simH))
simF <- simF/sum(simF)
simX <- data.frame(alleles=simH, freq=simF)

microbenchmark(
  prob.match(simX, 1000, match="8/8", twoSides = F),
  prob.matchR(simX, 1000, match="8/8"),
  times = 5L
)
#Unit: seconds
#nx: 108                                   expr      min      lq       mean     median   uq       max       neval cld
#prob.match(simX, 1000, match = "8/8", twoSides = F) 3.045275 3.058261 3.079038 3.062833 3.074749 3.154069     5  a
#prob.matchR(simX, 1000, match = "8/8")              49.81076 49.81186 49.96077 49.89218 49.98675 50.30231     5   b



simA <- crossing(c("A*01", "A*02", "A*03", "A*11", "A*24", "A*30"),
                 c("B*07", "B*13", "B*51", "B*57"),
                 c("C*01", "C*02", "C*03", "C*16"),
                 c("DRB1*01", "DRB1*02", "DRB1*04", "DRB1*11"))
simH <- apply(simA, 1, paste, collapse="~")
set.seed(20220503); simF <- runif(length(simH))
simF <- simF/sum(simF)
simX <- data.frame(alleles=simH, freq=simF)

microbenchmark(
  prob.match(simX, 1000, match="8/8", twoSides = F),
  prob.matchR(simX, 1000, match="8/8"),
  times = 5L
)
#Unit: seconds
#nx: 384                                   expr      min       lq        mean      median    uq         max       neval cld
#prob.match(simX, 1000, match = "8/8", twoSides = F) 159.6931  160.4584  161.9433  162.5754  162.7357  164.254      5  a
#prob.matchR(simX, 1000, match = "8/8")              855.2778  856.0939  861.8312  856.9627  857.8357  882.9861     5  b



simA <- crossing(c("A*01", "A*02", "A*03", "A*11", "A*24", "A*30"),
                 c("B*07", "B*13", "B*35", "B*40", "B*51", "B*57"),
                 c("C*01", "C*02", "C*03", "C*16", "C*17", "C*18"),
                 c("DRB1*01", "DRB1*02", "DRB1*04", "DRB1*11", "DRB1*14", "DRB1*16"))
simH <- apply(simA, 1, paste, collapse="~")
set.seed(20220503); simF <- runif(length(simH))
simF <- simF/sum(simF)
simX <- data.frame(alleles=simH, freq=simF)

microbenchmark(
  prob.match(simX, 1000, match="8/8", twoSides = F),
  prob.matchR(simX, 1000, match="8/8"),
  times = 5L
)
#Unit: seconds
#nx: 1296                                   expr      min       lq        mean      median    uq         max       neval cld
#prob.match(simX, 1000, match = "8/8", twoSides = F)  6937.772  6948.248  6953.98  6955.567  6960.043  6968.273     5  a
#prob.matchR(simX, 1000, match = "8/8")               15234.740 15251.060 15262.84 15256.191 15257.944 15314.248     5   b
