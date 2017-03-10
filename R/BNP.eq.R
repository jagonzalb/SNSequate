#' @name BNP.eq.predict
#' 
#' @title Prediction step for Bayesian non-parametric model for test equating
#' 
#' @description The 
#' 
#' @param model A 'BNP.eq' object.
#' @param from Numeric. A vector of indices indicating from which patterns equating should be performed.
#'        The covariates involved are integrated out.
#' @param into Numeric. A vector of indices indicating into which patterns equating should be performed.
#'        The covariates involved are integrated out.
#' @param alpha Numeric. Significance for credible bands.
#'
#' @details Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
#' Vivamus finibus vitae eros quis dictum. Donec lacus risus, facilisis quis 
#' tincidunt et, tincidunt sed mi. Nullam ullamcorper eros est, sed 
#' fringilla metus volutpat eu. Etiam ornare nulla id lorem posuere, eu vehicula 
#' urna vestibulum. Quisque luctus, diam ac mattis faucibus, leo felis tincidunt urna, eu 
#' tempus massa neque nec nibh. Aliquam erat volutpat. Fusce tempor mattis enim quis pretium. 
#' Aliquam volutpat luctus felis, nec fringilla enim tincidunt sed. 
#' Nam nec leo quis erat lobortis vulputate ac at neque
#' 
#' @return  
#'  A 'BNP.eq.predict' object, which is list containing the following items:
#' @return pdf A list of PDF's.
#' @return cdf A list of CDF's.
#' @return equ Numeric. Equated values.
#' @return grid Numeric. Grid used to evaluate pdf's and cdf's.
#' 
#' @references asdasd
#' 
#' @author Daniel Leon A. \email{dnacuna@mat.uc.cl}, Felipe Barrientos \email{afb26@stat.duke.edu}.
#' 
#' @keywords BNP equating, Bayesian non-parametrics, equating
BNP.eq.predict <- function(model, from=NULL, into=NULL, alpha=0.05){
  if(class(model) != "BNP.eq")
    stop("Fitted object must be BNP.eq.")

  fit <- model$fit
  patterns <- model$patterns
  patterns_freq <- model$patterns_freq
  prior <- model$fit$prior
  mcmc <- model$fit$mcmc
  Y <- model$Y
  X <- model$X
  max_score <- model$max_score

  patterns_cols <- ncol(patterns)
  patterns_rows <- nrow(patterns)
  
  grid <- seq(0, max_score, by=0.25) / max_score  ## This should be an input
  # grid <- seq(0, max_score, length.out=100)
  
  requireNamespace("progress")
  
  ## Only compute PDF and CDF if has not been done before
  if(!exists("sample_cdf", where=model) || !exists("sample_pdf", where=model)){
    ##################################
    # Prediction of densites and CDFs#
    ##################################
    cat("Computing PDF/CDF's and weights for BNP model. \n")
    
    Density <- list()
    CDF <- list()

    k <- fit$save.state$thetasave[,1]

    # Atoms
    #                         Max N
    #                  ------------------->
    #                 |
    #    MCMC Sample  |   Sample Matrix
    #                 |
    #                 v  
    tmpb <- fit$save.state$randsave[, 1:(patterns_cols*(prior$maxn-1))]
    theta <- fit$save.state$randsave[, (patterns_cols*(prior$maxn-1)+1):(patterns_cols*(prior$maxn-1)+prior$maxn)]
    theta <- ceiling(matrix(k, ncol=prior$maxn, nrow=mcmc$nsave)*theta)

    b <- lapply(split(tmpb, ceiling(col(tmpb)/patterns_cols)),
                matrix, ncol=patterns_cols)
    
    w_samples <- list()
    
    
    for(j in 1:patterns_rows){
      x <- as.numeric(patterns[j,])

      # Weights
      v <- do.call(cbind,
                   lapply(1:length(b), function(i) 1/(1+exp(-as.vector(b[[i]]%*%x)))))
      v <- cbind(v, 1) # Last column have only 1's

      w <- v;
      w[,2] <- v[,2]*(1-v[,1])

      for(i in 3:prior$maxn)
        w[, i] <- v[, i] * apply(1-v[,1:(i-1)], 1, prod)
      
      w_samples[[j]] <- w  # Save samples for equating process
      
      # Computing pdf or cdf
      f1 <- function(i, y, type="pdf") {
        tmpt <- theta[i, ]
        tmpw <- w[i, ]
        tmpk <- k[i]
        if(type == "pdf"){
          sum(tmpw * dbeta(y, tmpt, tmpk-tmpt+1))
        } else if(type == "cdf"){
          sum(tmpw * pbeta(y, tmpt, tmpk-tmpt+1))
        }

      }
      
      f11 <- function(y, nsave, type) {
        pb$tick()
        apply(cbind(1:nsave), 1, FUN=f1, y=y, type=type)
      }

      ## PDF
      pb <- progress_bar$new(
        format = paste0("Computing PDF ", j, "/", patterns_rows, ": [:bar] :percent eta: :eta"),
        total = length(grid)
      )
      pb$tick(0)
      tmp.pdf <- lapply(grid, FUN=f11, mcmc$nsave,type="pdf")
      dens <- do.call(cbind, tmp.pdf)

      ## CDF
      pb <- progress_bar$new(
        format = paste0("Computing CDF ", j, "/", patterns_rows, ": [:bar] :percent eta: :eta"),
        total = length(grid)
      )
      pb$tick(0)
      tmp.cdf <- lapply(grid, FUN=f11, mcmc$nsave, type="cdf")
      cdf <- do.call(cbind, tmp.cdf)

      Density[[j]]<-dens
      CDF[[j]]<-cdf
    }

    eval.parent(substitute(model[["sample_pdf"]] <- Density))
    eval.parent(substitute(model[["sample_cdf"]] <- CDF))
    eval.parent(substitute(model[["k"]] <- k))
    eval.parent(substitute(model[["b"]] <- b))
    eval.parent(substitute(model[["theta"]] <- theta))
    eval.parent(substitute(model[["w"]] <- w_samples))

  } else {
    cat("Retrieving cached results. \n")
    Density <- model$sample_pdf
    CDF <- model$sample_cdf
    k <- model$k
    b <- model$b
    theta <- model$theta
    w_samples <- model$w
  }


  ##############################################
  # Prediction of cdfs - Intregating F1 and F2 #
  ##############################################
  
  # Function to integrate samples
  integ_sample <- function(from, sample, fix_bounds=TRUE){
    CDFI <- Reduce('+', lapply(from, function(i) patterns_freq[i] * sample[[i]]) ) / sum(patterns_freq[from])
    if(fix_bounds){
      CDFI <- ifelse(CDFI > 1, 1, CDFI)
      CDFI <- ifelse(CDFI < 0, 0, CDFI)
    }
    CDFI
  }
  
  # Compute credible bands
  credible_bands <- function(sample, alpha){
    list(sup=apply(sample, 2, quantile, probs=1-(alpha/2)), 
         inf=apply(sample, 2, quantile, probs=alpha/2))
  }
  
  # If 'from' is not null, integrate that
  if(!is.null(from)){
    if(!all(from %in% 1:patterns_rows))
      stop("Indexes where integration occurs must exist.")
    
    CDFI <- integ_sample(from, CDF)
    
  } else {  # Integrate everything
    CDFI <- integ_sample(1:patterns_rows, CDF)
  }
  
  if(!is.null(into)){  ## If there is 'into' parameter, do equating
    if(!all(into %in% 1:patterns_rows))
      stop("Indexes where integration occurs must exist.")
    
    # Function to get an averaged point sample from 'into' cdf
    integ_sample_point <- function(into, Z, row) {
      if(Z == 0)
        return(0)
      if(Z == 1)
        return(1)
      
      t <- as.numeric(theta[row, ])
      k <- as.numeric(k[row])
      
      tmp <- sum(sapply(into, function(i) {
        patterns_freq[i] * sum(w_samples[[i]][row, ] * pbeta(Z, t, k-t+1))
      }))
      tmp / sum(patterns_freq[into])
    }
    
    # Find root to equate
    f2 <- function(i1, i2, cdf) {
      g <- function(Z) {
        integ_sample_point(into, Z, i1) - CDFI[i1, i2]
      }
      
      uniroot(g, c(0, 1))$root
    }
    
    # Basically the 'i' in a double 'for'
    f21 <- function(i2, nsave) {
      pb$tick()
      
      apply(cbind(1:nsave), 1, f2, i2=i2)
    }
    ## TODO: Add parallel comp if available
    # The 'j' in a double 'for'
    #                            Grid Points (i2)
    #                     --------------------------->
    #                     |
    #   MCMC sample (i1)  |     CDF Sample Matrix
    #                     |
    #                     v
    pb <- progress_bar$new(
      format = "Computing equated values: [:bar] :percent eta: :eta",
      total = length(grid)
    )
    pb$tick(0)
    equ <- do.call(cbind, lapply(1:length(grid), FUN=f21, mcmc$nsave))
    
    pdf_X <- integ_sample(from, Density, fix_bounds=FALSE)
    pdf_Y <- integ_sample(into, Density, fix_bounds=FALSE)
    
    cdf_X <- CDFI
    cdf_Y <- integ_sample(into, CDF)
    
    tmp <- apply(equ,2,mean)
    eYx <- tmp[c(1, (1:max_score)*4+1)]*max_score  # Skinning the grid into the actual scores
    
    ret <- list(pdf_X=list(pdf=(apply(pdf_X, 2, mean) ), cb=credible_bands(pdf_X, alpha)),
                pdf_Y=list(pdf=(apply(pdf_Y, 2, mean) ), cb=credible_bands(pdf_Y, alpha)),
                cdf_X=list(cdf=(apply(cdf_X, 2, mean) ), cb=credible_bands(cdf_X, alpha)),
                cdf_Y=list(cdf=(apply(cdf_Y, 2, mean) ), cb=credible_bands(cdf_Y, alpha)),
                equ=list(equ=tmp, cb=credible_bands(equ, alpha)), 
                eYx=eYx, grid=grid, max_score=max_score)
    class(ret) <- 'BNP.eq.predict'
    ret
    
  } else {  ## If there is no 'into', only integrate and return PDF/CDF
    pdf_X <- integ_sample(from, Density)
    cdf_X <- CDFI
    
    ret <- list(pdf_X=list(pdf=apply(pdf_X, 2, mean), cb=credible_bands(pdf_X, alpha)),
                cdf_X=list(cdf=apply(cdf_X, 2, mean), cb=credible_bands(cdf_X, alpha)),
                grid=grid, max_score=max_score)
    class(ret) <- 'BNP.eq.predict'
    ret
  }
  
}

print.BNP.eq.predict <- function(x, ...) {
  if(exists("eYx", where=x)){
    cat("Equated values: \n\n")
    df <- data.frame(Score=0:x$max_score, eYx=x$eYx)
    print(df, row.names=FALSE, digits=3)
  }
    
}

plot.BNP.eq.predict <- function(x, which="Equating", what=NULL, add=FALSE, ...){

  max_score <- x$max_score
  yname <- ""
  xname <- "X"
  
  if(which == "X"){
    if(what == "PDF"){
      data <- x$pdf_X$pdf * 12
      cb <- x$pdf_X$cb
      yname <- "f(X)"
    } else if(what == "CDF"){
      data <- x$cdf_X$cdf
      cb <- x$cdf_X$cb
      yname <- ""
    }
    
  } else if(which == "Y") {
    if(what == "PDF"){
      data <- x$pdf_Y$pdf
      cb <- x$pdf_Y$cb
      yname <- "f(Y)"
    } else if(what == "CDF"){
      data <- x$cdf_Y$cdf
      cb <- x$cdf_Y$cb
      yname <- ""
    }
    xname <- "Y"
  } else if(which == "Equating") {  ## Only equating
    if(!is.null(what)) {  ## Integrated values only
      if(what == "PDF") {
        data <- x$pdf_X$pdf
        cb <- x$pdf_X$cb
        yname <- "f(X)"
      } else if(what == "CDF"){
        data <- x$cdf_X$cdf
        cb <- x$cdf_X$cb
        yname <- ""
      }
    } else {
      data <- x$equ$equ * max_score
      cb <- x$equ$cb
      cb$sup <- cb$sup * max_score
      cb$inf <- cb$inf * max_score
      yname <- expression(varphi[X])
    }
  }
  
  tmpx <- c(x$grid, rev(x$grid), x$grid[1]) * max_score
  tmpy <- c(cb$sup, rev(cb$inf), cb$sup[1])
  
  if(add) {
    polygon(tmpx, tmpy, col=rgb(0, 0, 1, 0.2), border=NA)
    lines(x$grid*max_score, data, lwd=2, type='l', ...)
  } else {
    plot(x$grid*max_score, data, lwd=2, type='l', ylim=c(0, max(tmpy, data)),
         xlab=xname, ylab=yname)
    polygon(tmpx, tmpy, col='gray', border=NA)
    lines(x$grid*max_score, data, lwd=2, type='l', ...)
  }
  
}

#' @name BNP.eq
#' 
#' @title Bayesian non-parametric model for test equating
#' 
#' @description The Bayesian nonparametric (BNP) approach (Ghoshand Ramamoorthi, 2003; Hjort et al., 2010) 
#' starts by focusing on spaces of distribution functions, so that uncertainty is expressed on F itself. 
#' The prior distribution p(F) is defined on the space F of all distribution functions defined on X . If X 
#' is an infinite set then F is infinite-dimensional, and the corresponding prior model 
#' p(F) on F is termed nonparametric. The prior probability model is also referred to
#' as a random probability measure (RPM), and it essentially corresponds to a distribution 
#' on the space of all distributions on the set X . Thus Bayesian nonparametric models 
#' are probability models defined on a function space (Muller and Quintana, 2004).
#' 
#' Gonzalez et al. (2015) proposed a Bayesian non-parametric approach for equating. The main
#' idea consists of introducing covariate dependent BNP models for a collection of 
#' covariate-dependent equating transformations
#' 
#'  \eqn{ \left\{ \boldsymbol{\varphi}_{\boldsymbol{z}_f, \boldsymbol{z}_t} (\cdot): 
#'            \boldsymbol{z}_f, \boldsymbol{z}_t \in \mathcal{L}
#'        \right\} 
#'   }
#' 
#' 
#' 
#' 
#' @param scores_x Vector.  Scores of form X.
#' @param scores_y Vector.  Scores of form Y.
#' @param range_scores Vector of length 2.  Represent the minimum and maximum scores in the test.
#' @param design Character.  Only supports 'EG' design now.
#' @param covariates Data.frame.  A data frame with factors, containing covariates 
#'        for test X and Y, stacked in that order.
#' @param prior List.  Prior information for BNP model. 
#' @param mcmc List.  MCMC information for BNP model. 
#' @param normalize Logical.  Whether normalize or not the 
#'        response variable. This is due to Berstein's polynomials. Default is TRUE.
#'
#' @details Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
#' Vivamus finibus vitae eros quis dictum. Donec lacus risus, facilisis quis 
#' tincidunt et, tincidunt sed mi. Nullam ullamcorper eros est, sed 
#' fringilla metus volutpat eu. Etiam ornare nulla id lorem posuere, eu vehicula 
#' urna vestibulum. Quisque luctus, diam ac mattis faucibus, leo felis tincidunt urna, eu 
#' tempus massa neque nec nibh. Aliquam erat volutpat. Fusce tempor mattis enim quis pretium. 
#' Aliquam volutpat luctus felis, nec fringilla enim tincidunt sed. 
#' Nam nec leo quis erat lobortis vulputate ac at neque
#' 
#' @return  
#'   A 'BNP.eq' object, which is list containing the following items:
#' @return Y Response variable.
#' @return X Design Matrix.
#' @return fit DPpackage object. Fitted model with raw samples.
#' @return max_score Maximum score of test.
#' @return patterns A matrix describing the different patterns formed
#'         from the factors in the covariables.
#' @return patterns_freq The normalized frequency of each pattern.
#' 
#' @references asdasd
#' 
#' @author Daniel Leon A. \email{dnacuna@mat.uc.cl}, Felipe Barrientos \email{afb26@stat.duke.edu}.
#' 
#' @keywords BNP equating, Bayesian non-parametrics, equating
BNP.eq <- function(scores_x, scores_y, range_scores=NULL, design="EG", covariates=NULL, 
                   prior=NULL, mcmc=NULL, normalize=TRUE){

  ## Response
  Y <- c(scores_x, scores_y)

  if(is.null(range_scores)){
    max_score=max(Y)
  } else{
    max_score=range_scores[2]
  }
  ## Bernstein poly assume response on [0, 1]
  if(normalize)
    Y <- Y / max_score

  ## Covariates.
  form_cov <- factor( c( rep("X", length(scores_x)),
                         rep("Y", length(scores_y))) )
  X <- model.matrix(~ ., data=cbind(Form=form_cov, covariates))

  ## Count patterns in the design matrix. Basically it
  ## extract each unique pattern from the discrete covariates
  Grid <- count(X)
  patterns <- Grid[, 1:(ncol(Grid)-1)]
  patterns_freq <- as.numeric(Grid$freq)
  patterns_freq <- patterns_freq/sum(patterns_freq)

  # Fitting single-atoms LDBPP model
  # mcmc parameters
  if(is.null(mcmc))
    mcmc <-  list(nburn = 5000, nskip = 20, ndisplay = 10, nsave = 1000)

  ## Default params, not needed for our purposes
  # Predictions will be made later
  grid <- 0.5
  npred <- 1
  xpred <-  patterns[1, ]  ## Just to get results. Prediction isn't made by DPpackge

  # List of prior information for DPpackage
  if(is.null(prior))
    prior <- list(maxn = 25, a1=1, a2=1,
                  lambda = 25, nu = 6,
                  psiinv = diag(1000, ncol(patterns)),
                  m0 = rep(0, ncol(patterns)),
                  S0 = diag(1000, ncol(patterns)))

  # State
  state <- NULL

  # fitting the model
  cat("Calling 'DPpackage' to sample model: \n\n")
  fit <- tLDBPPdensity(formula=Y~X[, -1],xpred=xpred,  # We take out the intercept because DPpackage adds it
                       prior=prior,
                       mcmc=mcmc,
                       state=state,status=TRUE,
                       grid=grid,
                       compute.band=FALSE,type.band="PD")

  ret <- list(scores_x=scores_x, scores_y=scores_y,
              Y=Y, X=X,fit=fit, max_score=max_score,
              patterns=patterns, patterns_freq=patterns_freq)
  class(ret) <- "BNP.eq"
  ret
}
