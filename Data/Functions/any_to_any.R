#' any2any
#'
#' Converts most common test statistics into most common (signed) effect sizes.
#'
#' @param st     Value(s) of a test statistic.
#' @param df1     Degrees of freedom
#' @param df2     NULL or degrees of freedom of the denominator for the f-distribution.
#' @param N     Number of data points used in calculation of test-statistic.
#' @param n1    Number of data points in sample 1.
#' @param n2    Number of data points in sample 2.
#' @param esType     Type of test statistic. One of: "t", "lm.t", "f", "lm.f", "r", "X2", "Z", "lm.Z"
#' @param CIcalc     If \code{TRUE} (default) the Confidence Interval for the test statistic in \code{x} will be calculated using the "Confidence limits for noncentral parameters" functions in package (e.g., for type - "t": \link[MBESS]{conf.limits.nct}).
#' @param CL    Confidence Limit (default: .95).
#' @param rID    Correlation among predictor values in a linear model.
#' @param q      Number of predictors in the model.
#' @param alternative     Alternative hypothesis (defult = "two").
#' @param keepSign     Return effect size with sign of test statistic? (default = TRUE).
#' @param keepSignNames     Which effect sizes should keep the sign if \code{keepSign = TRUE}? Default is to keep the sign for: "r","l.r","u.r","fisher.z","l.z","u.z".
#'
#' @details The procedure to calculate a variety of effect sizes is as follows:
#'
#' \itemize{
#' \item If \code{CIcalc == FALSE}, \code{package::compute.es} will be used to convert the test statistic to a large number of effect size estimates. The confidence intervals around the effect size estimates will be based meta-analytic estimates of effect size variance (e.g., for type - "t": \link[compute.es]{tes}).
#' \item If \code{CIcalc == TRUE}, \code{package::MBESS} will be used to calculate the confidence interval for the test statistic based on its noncentral distribution (e.g., for type - "t": \link[MBESS]{conf.limits.nct}). Subsequently the test statistic, as well as its lower and upper confidence limit will each be passed to \code{compute.es} seperately.
#' \item If \code{keepSign == TRUE} the sign of the test statistic will be copied to all the effect sizes in \code{keepSignNames}.
#' }
#'
#' @note The prefix "lm" is currently disregarded, but will be implemented in future versions to indicate the test statistic is in fact a fixed factor in a linear model.
#'
#' @author
#' Fred Hasselman (inspired by RP:P function \code{any2r} by CHJ Hartgerink)
#'
#' @return The effect sizes calculated by \code{compute.es} corresponding to the test statistic(s), with either meta-analytic, or, exact CI.
#'
any2any <- function(testInfo,
                    df1 = NULL,
                    df2 = NULL,
                    N   = NULL, n1 = NULL, n2 = NULL,
                    esType  = NA,
                    var.lor = NA,
                    CIcalc  = TRUE, CL = .95, rID = 0, q = 1,
                    alternative   = "two", keepDirection = TRUE,
                    keepSign      = TRUE,
                    keepSignNames = c("r","l.r","u.r","fisher.z","l.z","u.z")){
   require(MBESS)
   require(compute.es)
  esType.cl <- NA

  ifelse(grepl("two",alternative),{alternative <<- "two"},{alternative <<- "one"})

  #  ifelse(any(grepl("estimate",colnames(testInfo))), {st <- testInfo$estimate},{st <- testInfo$statistic})
  if(any(grepl("(.r)+", esType),(esType%in%c("OR")))){
    st <- testInfo$estimate
  } else {
    st <- testInfo$statistic
  }

  if(grepl("Z.f", esType, fixed = TRUE)){
    if(length(testInfo$estimate)>1){
      esType <- esType.cl <- "Z"
      st <- testInfo$statistic
    }
  }


  if(grepl("Z", esType, fixed = TRUE)){
    if(alternative=="one"){
      n1<-N/2
      n2<-N/2
      alternative<-"two"
    }
  }


  if(is.null(st)){stop("No test statistic to caclulate ES-CI.")}

  # Use Cohen's dz for paired t-test
  if(esType%in%"t.p"){
    n1 <- n2 <- N/2
    st <- testInfo$statistic / sqrt(N)
  }


  # Check for model-based t-statistics
  if(grepl("(lm.t)+",esType)){
    # Make df of the predictor
    n1 <- n2 <- (df1+1)/2
  }


  # Treat model based stats as 'regular'
  esType.cl <- gsub("(lm.)","",esType)


  # Settings for model-based OR
  if(grepl("(lm.OR)+",esType)){
    st        <- exp(st)
    esType.cl <- "Asym"
    CIcalc    <- FALSE
  }

  if(is.na(esType.cl)){esType.cl<-esType}

  if(is.na(st)|is.null(st)|is.nan(st)){

    ES <- compute.es::tes(t=2, level= 95, n.1 = 100, n.2 = 100, verbose = FALSE, dig = 5)
    ES[seq_along(ES)] <- NA
    ES <- c(st,st,st,ES)

    CIcalc <- FALSE
  }


  if(CIcalc){

    getCI <- TRUE

    if(esType=="OR"){
      # Fisher exact test gives exact noncentral hypergeometric CI
      sCI <- cbind(ncp    = testInfo$estimate,
                   ncp.lo = testInfo$conf.low,
                   ncp.hi = testInfo$conf.high)
      getCI <- FALSE
    }

    if((grepl("Z.f", esType, fixed = TRUE))&(length(testInfo$estimate)==1)){
      sCI <- cbind(ncp    = testInfo$estimate,
                   ncp.lo = testInfo$conf.low,
                   ncp.hi = testInfo$conf.high)
      esType <- esType.cl <- "r"
      getCI <- FALSE
    }

    if(getCI){
      sCI <- get.ncpCI(st, df1, df2, N, esType.cl, CL, keepSign, alternative)
      if(esType=="f"){sCI[1,is.na(sCI)]<-1}
      if(esType%in%c("t","t.p","t.r","Z")){sCI[1,is.na(sCI)]<-0}
    }
    # no CI
  } else {
    sCI <- cbind(ncp  = st)
  }

  esComp <- list()
  cnt    <- 0

  for(cnt in seq_along(sCI)){

    x <- sCI[cnt]
    if(x==0|is.na(x)|is.nan(x)){
      x <- rnorm(1)*1e-12
      disp("ES converison: A test statistic of 0 (or NA, or NaN) was changed to rnorm(1) * 1e-12 in order to enable ES conversion.", header= FALSE, footer = FALSE)}

    # This effectively ignores model based stats
    esType <- gsub("lm.","",esType,fixed=TRUE)
    #esType <- gsub("Z.f","r",esType,fixed=TRUE)

    ncCI <- list()
    switch(esType,
           t.p  = esComp[[cnt]] <- compute.es::des(d   = x, n.1 = (N/2), n.2 = (N/2), level=CL*100, verbose = FALSE, dig = 5),
           t    = esComp[[cnt]] <- compute.es::tes(t   = x, level=CL*100,
                                                   n.1 = n1, n.2 = n2, verbose = FALSE, dig = 5),
           lm.t = esComp[[cnt]] <- compute.es::a.tes(t=x, level=CL*100,
                                                     n.1 = n1, n.2 = n2, R = rID, q = q,
                                                     verbose = FALSE, dig = 5),
           t.r  = esComp[[cnt]] <- compute.es::res(r = x, level=CL*100, var.r = ((1-x^2)^2)/(N-1),
                                                   n = N, verbose = FALSE, dig = 5),
           r  = esComp[[cnt]] <- compute.es::res(r = x, level=CL*100, var.r = NULL,
                                                 n = N, verbose = FALSE, dig = 5),
           #compute.es::res(r=x, level=CL, n=N, verbose = FALSE, dig = 5),
           f    = esComp[[cnt]] <- compute.es::fes(f=x, level=CL*100,
                                                   n.1 = n1, n.2 = n2, verbose = FALSE, dig = 5),
           lm.f = esComp[[cnt]] <- compute.es::a.fes(f=x, level=CL*100,
                                                     n.1 = n1, n.2 = n2, R = rID, q = q,
                                                     verbose = FALSE, dig = 5),
           X2   = esComp[[cnt]] <- compute.es::chies(chi.sq = x, level = CL*100,
                                                     n = N, verbose = FALSE, dig = 5),
           Z    = esComp[[cnt]] <- compute.es::pes(p = pnorm(abs(x), lower.tail= FALSE)*2, level = CL*100,
                                                   n.1 = n1, n.2 = n2, tail = "two", verbose = TRUE, dig = 5),
           lm.Z  = esComp[[cnt]] <- compute.es::a.pes(p = pnorm(abs(x), lower.tail= FALSE)*2, level = CL*100,
                                                      n.1 = n1, n.2 = n2, R = rID, q = q,
                                                      tail = alternative, verbose = FALSE, dig = 5),
           OR    = esComp[[cnt]] <- compute.es::lores(lor=log(x), n.1 = n1, n.2 = n2,
                                                      var.lor = var.lor, verbose = FALSE, dig = 5, level = CL*100)
    )
  }

  # This section re-calculates CI based on the exact CI for the test statistic obtained from MBESS in function get.ncpCI
  if(cnt>1){

    if(esType%in%c("r","t.r")){

      ncp <- compute.es::tes(t=2, level= 95, n.1 = 100, n.2 = 100, verbose = FALSE, dig = 5)
      ncp[seq_along(ncp)] <- NA
      id.l <- c("l.d","l.r", "l.z", "l.or", "l.lor")
      id.u <- c("u.d","u.r", "u.z", "u.or", "u.lor")
      id.e <- c("d", "r", "fisher.z", "OR", "lOR")
      rNames <- names(compute.es::res(r=1,var.r=.5, n=100, level=95,dig=5,verbose = FALSE))
      ncp[,rNames] <- esComp[[1]][,rNames]
      ncp$N.total <- N
      ncp$n.1 <- n1
      ncp$n.2 <- n2

      if(esType=="t.r"){
        sCI <- get.ncpCI(testInfo$statistic, df1, df2, N, "t", CL, keepSign)
      }
    } else {

      ncp  <- esComp[[1]]
      id.l <- c("l.d", "l.g", "l.r", "l.z", "l.or", "l.lor")
      id.u <- c("u.d", "u.g", "u.r", "u.z", "u.or", "u.lor")
      id.e <- c("d", "g", "r", "fisher.z", "OR", "lOR")
    }

    ncp[,id.l] <- esComp[[2]][,id.e]
    ncp[,id.u] <- esComp[[3]][,id.e]

    ES <- cbind(sCI, ncp[,colnames(ncp)!="NNT"])

  } else {

    if(esType.cl%in%"Asym"){
      sCI <- cbind(ncp    = esComp[[1]]$lOR,
                   ncp.lo = esComp[[1]]$l.lor,
                   ncp.hi = esComp[[1]]$u.lor)
    }
    ncp  <- esComp[[1]]
    ES   <- cbind(sCI, ncp[,colnames(ncp)!="NNT"])

  }

  colnames(ES)[1:cnt] <- c("ncp","ncp.lo","ncp.hi")[1:cnt]

  # compute.es keeps the sign for d and related ES, but not for r if tes and des are used,.
  # If keepSign = TRUE the sign from d will be copied to r and related es.

  # unique(ML2.key$stat.type)
  #  "t"    "t.r"  "OR"   "lm.t" "Z"    "f"    "lm.Z"
  if(!all((sign(ES$ncp)==sign(ES[ ,c("d","r")])),(sign(ES$ncp.lo)==sign(ES[ ,c("l.d","l.r")])),(sign(ES$ncp.hi)==sign(ES[ ,c("u.d","u.r")])), na.rm = TRUE) & !esType%in%c("OR","t.r","r")){
    if(keepSign){
      if(esType%in%c("X2","f")){
        id.l <- which(colnames(ES) %in% c("l.d", "l.g", "l.r", "l.z"))
        id.u <- which(colnames(ES) %in% c("u.d", "u.g", "u.r", "u.z"))
        id.e <- which(colnames(ES) %in% c("d", "cliffs.d", "g", "r", "fisher.z"))
        col.id <- c(id.e,id.l,id.u)
      }
      if(any(esType%in%c("lm.Z","Z"))){ # esType=="Z"|esType=="lm.Z"
        col.id <-which(colnames(ES)%in%c("d","l.d","u.d",keepSignNames))
      }
      if(esType%in%c("t","lm.t")){
        col.id <-which(colnames(ES)%in%keepSignNames)
      }
      col.id <- sort(col.id)
      ES[ ,col.id] <- ES[ ,col.id] * sign(sCI)[1:cnt]
    }
  }

  return(ES)
}