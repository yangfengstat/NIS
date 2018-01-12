penGAM <- function(x, y, lambda.pen, lambda.curv, knots,
                   model = LinReg(), save.x = TRUE,
                   control = grpl.control()){
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 13 Mar 2008, 11:30

  if(any(lambda.pen > 1))
    stop("lambda.pen is on a *relative* scale with respect to lambda.max")

  n <- NROW(x)
  p <- NCOL(x)

  ## Construct the design matrix and the penalty matrices
  d <- getDesign(x, knots = knots)

  ## Extract the necessary information
  m.orig <- d$m     ## list of p design matrices
  pen    <- d$pen   ## list of p penalty matrices
  index  <- d$index ## index vector (without intercept)

  ## Coefficients for all combinations of lambda.pen * lambda.curv * coef
  coefficients <- array(0,
                        dim = c(length(lambda.pen),
                                length(lambda.curv),
                                length(index) + 1))

  lambda.max.all <- numeric(length(lambda.curv))

  for(u in 1:length(lambda.curv)){
    if(control@trace >= 1)
      cat("\n*** lambda.curv", lambda.curv[u], "***\n")

    lcurv <- lambda.curv[u]

    ## m contains the whole design matrix but *without* intercept column
    m           <- matrix(0, nrow = n, ncol = length(index))
    colnames(m) <- unlist(lapply(m.orig, colnames))

    ## We have to calculate the penalty matrix for each value of lambda.curv
    ## Transform each predictor group such that we can work with an
    ## "ordinary" Group Lasso model

    chol.decomp <- list(); length(chol.decomp) <- p

    for(j in 1:p){
      ## This will be the penalty matrix for each block
      penmat  <- lcurv * pen[[j]] + crossprod(m.orig[[j]])

      ## Now calculate the cholesky decomposition
      chol.decomp[[j]] <- chol(penmat)
      m[,index == j] <-
        t(forwardsolve(t(chol.decomp[[j]]), t(m.orig[[j]])))
    }

    ## Now we are ready to get the estimates using the grplasso function
    ## for *all* values of lambda.pen

    ## We do *not* have to rescale the penalty
    mypenscale <- function(x) rep(1, length(x))

    ## Calculate lambda.max
    lambda.max <- lambdamax(x = cbind(1, m), y = y, index = c(NA, index),
                            penscale = mypenscale, model = model,
                            standardize = FALSE)

    lambda.max.all[u] <- lambda.max

    ## Actual lambda values
    lambda.use <- lambda.max * lambda.pen

    coef.init    <- rep(0, NCOL(m) + 1)
    coef.init[1] <- model@link(mean(y))

    ## Call now group lasso algorithm on the *transformed* data
    fit.grp <- grplasso(x = cbind(1, m), y = y, index = c(NA, index),
                        lambda = lambda.use, model = model,
                        coef.init = coef.init,
                        standardize = FALSE, penscale = mypenscale,
                        control = control)

    ## Transform the estimates back to the original scale using backsolve
    coef.pre <- coef(fit.grp)[-1,,drop = FALSE] ## without intercept
    coef.out <- matrix(0, nrow = NROW(coef.pre), ncol = NCOL(coef.pre))
    for(j in 1:p){
      sel            <- index == j
      coef.out[sel,] <- backsolve(chol.decomp[[j]], coef.pre[sel,])
    }
    ## Now transform the values back to the original spline-basis and
    ## add intercept.
    ## Put everything into the "coefficients" array.
    coefficients[,u,] <- t(rbind(coef(fit.grp)[1,], coef.out))
  } ## end for(u in 1:length(lambda.curv))

  dimnames(coefficients) <- list(lambda.pen,
                                 lambda.curv,
                                 c("(Intercept)", colnames(m)))

  if(is.null(colnames(x))){
    colnam <- 1:NCOL(x)
  }else{
    colnam <- colnames(x)
  }

  rownam <- rownames(x)

  if(!save.x)
    x <- NULL

  out <- list(coefficients = coefficients,
              x            = x,
              colnames     = colnam,
              rownames     = rownam,
              splineobj    = list(basis = d$basis, index = d$index),
              lambda.max   = lambda.max.all,
              lambda.pen   = lambda.pen,
              lambda.curv  = lambda.curv,
              model        = model)
  structure(out, class = "penGAM")
}





predict.penGAM <- function(object, newdata, type = c("link", "response"), ...)
{
  ## Purpose: Predict fitted object on a new data set
  ## ----------------------------------------------------------------------
  ## Arguments: fit object
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 13 Mar 2008, 12:09

  type <- match.arg(type)

  nr.lambda.pen  <- length(object$lambda.pen)
  nr.lambda.curv <- length(object$lambda.curv)

  if(missing(newdata) || is.null(newdata)){ ## if we don't have new data object
    if(is.null(object$x))
      stop("No x matrix found in 'object'. Set save.x = TRUE or provide newdata!")

    ## Create the whole design matrix (*all* components)
    ## Re-create the design matrix using the given spline basis
    m.list <- list()
    for(j in 1:NCOL(object$x))
      m.list[[j]] <- eval.basis(x[,j], basisobj = object$splineobj$basis[[j]])

    n    <- NROW(m.list[[1]])
    pred <- array(0, dim = c(nr.lambda.pen, nr.lambda.curv, n))

    ## Create predictions
    m <- matrix(unlist(m.list), nrow = n)
    for(l in 1:nr.lambda.pen){
      for(u in 1:nr.lambda.curv){
        pred[l,u,] <- cbind(1, m) %*% coef(object)[l,u,]
      }
    }
    dimnames(pred) <- list(object$lambda.pen, object$lambda.curv,
                           object$rownames)
  }else{ ## if we have *new* data object (!!!add some test-functions later!!!)
    x <- newdata
    n <- NROW(x); p <- NCOL(x)

    pred <- array(0, dim = c(nr.lambda.pen, nr.lambda.curv, n))
    dimnames(pred) <- list(object$lambda.pen, object$lambda.curv, rownames(x))
    lowdiff <- highdiff <- outrange <-
      matrix(FALSE, nrow = NROW(newdata), ncol = p)

    ## Perform *linear* extrapolation

    ## Create the 'cut' design matrix: Sets values of x which exceed the
    ## training range to their corresponding boundaries (given by the range
    ## of the training data)
    x.cut  <- matrix(0, nrow = n, ncol = p)
    m.list <- m.deriv.list <- list()
    length(m.list) <- length(m.deriv.list) <- p

    lower.slopes <- upper.slopes <- list()

    for(j in 1:p){
      ## Extract current predictor
      x.current     <- x[,j]
      x.cut[,j]     <- x.current

      ## Current basis matrix (B-spline basis)
      bas           <- object$splineobj$basis[[j]]

      ## Get the smallest and largest knots (=boundaries)
      lower.end     <- bas$rangeval[1]
      upper.end     <- bas$rangeval[2]

      ## Which observations are exceeding the boundary-values?
      ind.lower     <- x.current < lower.end
      ind.upper     <- x.current > upper.end

      lowdiff[ind.lower,j]  <- (x.current - lower.end)[ind.lower]
      highdiff[ind.upper,j] <- (x.current - upper.end)[ind.upper]

      x.cut[ind.lower,j] <- lower.end
      x.cut[ind.upper,j] <- upper.end

      ## Get the slopes at the boundaries
      m.list[[j]] <- eval.basis(x.cut[,j], bas)
      deriv.info  <- eval.basis(c(lower.end, upper.end), bas, Lfdobj = 1)

      lower.slopes[[j]] <- deriv.info[1,]
      upper.slopes[[j]] <- deriv.info[2,]
    }
    m     <- matrix(unlist(m.list), nrow = n) ## Gets filled *columnwise*
    index <- object$splineobj$index ## without intercept

    for(l in 1:nr.lambda.pen){
      ##cat(l, "\n")
      for(u in 1:nr.lambda.curv){
        beta <- coef(object)[l,u,]

        pred.pre <- cbind(1, m) %*% coef(object)[l,u,]

        ## Put the design matrix of the first derivates (lower.slopes,
        ## upper.slopes) into one long vector each (same length as index) and
        ## multiply with beta vector and take the sum. I.e. perform the matrix
        ## operation in a bit a special form.
        ## The result are the derivatives at the left- and the right-hand side
        ## boundaries of the training range (of the fitted object with the
        ## current coefficients)
        slopes.left  <- rowsum(unlist(lower.slopes) * beta[-1],
                               group = index)
        slopes.right <- rowsum(unlist(upper.slopes) * beta[-1],
                               group = index)

        ## Now we have to multiply the derivatives with the difference
        ## in the x-values (contained in lowdiff and highdiff)
        ## lowdiff and highdiff are matrices with dimensions n x p, i.e. the
        ## dimension of the newdata object.
        ## Each column of lowdiff and highdiff is multiplied with the slope
        ## value. The result will be what we have to add beyond the boundaries.
        ## add.left and add.right will also have dimension n x p.

        ## 'as.array' is here to force a warning message if recycling would
        ## take place (see help file of sweep)
        add.left  <- sweep(lowdiff, MARGIN = 2, STATS = as.array(slopes.left),
                           FUN = "*")
        add.right <- sweep(highdiff, MARGIN = 2, STATS = as.array(slopes.right),
                           FUN = "*")

        ## Calculate the final prediction:
        ## Take the prediction of the 'cut-down' matrix and add the linear
        ## extrapolation part (add.left + add.right). We have to take the sum
        ## in each row of the linear extrapolation part (add.left + add.right)
        pred[l,u,] <- pred.pre + rowSums(add.left + add.right)
      }
    }
  }

  ## Transform to original scale if necessary
  pred <- switch(type,
                 link     = pred,
                 response = object$model@invlink(pred))

  dims <- dim(pred)

  ## If we have only one combination of lpen/lcurv we allow drop = TRUE.
  ## In all other cases
  if(dims[1] == 1 & dims[2] == 1)
    pred <- pred[1,1,,drop = TRUE]

  attr(pred, "lambda.pen")  <- object$lambda.pen
  attr(pred, "lambda.curv") <- object$lambda.curv

  pred
  ##out <- list(pred = pred,
  ##            lambda.pen  = object$lambda.pen,
  ##            lambda.curv = object$lambda.curv)
  ##out
}

plot.penGAM <- function(x, which = NULL, ask = TRUE && dev.interactive(),
                        nrgrid = 100, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  6 May 2008, 12:08

  n      <- NROW(x$splineobj$m[[1]])
  index  <- x$splineobj$index

  p      <- index[length(index)]

  dims   <- dim(coef(x))
  if(dims[1] > 1 | dims[2] > 1)
    stop("Use plot(fit[i,j]) to select the desired penalty parameters")

  coefs <- coef(x)[1, 1,][-1] ## without intercept!

  if(is.null(which))
    which <- 1:p
  if(ask){
    op <- par(ask = TRUE)
    on.exit(par(op)) ## use original settings when we leave the function
  }

  for(j in 1:length(which)){
    bas       <- x$splineobj$basis[[j]]
    rng       <- bas$range
    x.current <- seq(rng[1], rng[2], length = nrgrid)
    m.current <- eval.basis(x.current, basisobj = bas)
    ind       <- index == which[j]
    fn        <- scale(m.current %*% coefs[ind], scale = FALSE)
    plot(sort(x.current), fn[order(x.current)], type = "l",
         xlab = x$colnames[which[j]], ylab = bquote(f[.(which[j])]), ...)
  }
}

"[.penGAM" <- function(x, i, j){
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  6 May 2008, 12:08

  ## First get dimensions of the original object x
  dim <- dim(coef(x))
  nrlambda.pen  <- dim[1]
  nrlambda.curv <- dim[2]

  if(missing(i))
    i <- 1:nrlambda.pen
  if(missing(j))
    j <- 1:nrlambda.curv

  ## Subset the object
  fit.red <- x

  fit.red$coefficients <- coef(x)[i,j,,drop = FALSE]

  ## We do not allow the situation where everything is removed
  if(length(fit.red$coefficients) == 0)
    stop("Not allowed to remove everything!")

  fit.red$lambda.max   <- x$lambda.max[j]
  fit.red$lambda.pen   <- x$lambda.pen[i]
  fit.red$lambda.curv  <- x$lambda.curv[j]

  fit.red
}

getDesign <- function(x, knots = 10, norder = 4)
{
  ## Purpose: Get design matrix and penalty matrix in the basis of
  ##          B-splines
  ## ----------------------------------------------------------------------
  ## Arguments: x:       Design matrix, *without* intercept column
  ##            nrknots: Number of quantiles at which to put knots
  ##                     (these will be evenly spaced on the quantile scale
  ##                      between 0 and 1)
  ##            norder:  Set norder = 4 for piecewise cubic polynomials
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  4 Mar 2008, 14:49

  p <- NCOL(x)

  if(is.null(colnames(x)))
    colnames(x) <- paste("x", 1:p, sep = "")

  ## If knots is only a number, we have to create the knots first for each
  ## predictor
  if(is.numeric(knots)){
    knots.use <- list()
    for(j in 1:p){
      x.current <- x[,j]
      ## Make sure that we do not place multiple knots at the same
      ## x value (if we have replicates etc.)
      ux.current <- unique(x.current)

      if(knots > length(ux.current)){
        nrknots.use <- length(ux.current)
        cat("Too many knots in predictor", j, "...\n")
        cat("Using as many knots as there are unique predictor values.\n")
      }else{
        nrknots.use <- knots
      }
      knots.use[[j]] <- quantile(ux.current, seq(0, 1, length = nrknots.use))
    }
  }else if(is.list(knots)){
    if(length(list) != p)
      stop("knots has wrong length")

    knots.use <- knots
  }

  ## Create the design and the penalty matrix
  m     <- list()
  pen   <- list()
  basis <- list()
  length(m) <- length(pen) <- length(basis) <- p

  df    <- numeric(p)

  for(j in 1:p){
    x.current <- x[,j]

    ## Create knot sequence
    breaks  <- knots.use[[j]]
    stopifnot(all.equal(range(breaks), range(x.current)))

    ## Get basis-object for the j-th spline
    basis[[j]] <- create.bspline.basis(rangeval = range(breaks),
                                       breaks = breaks, norder = norder)
    ## removed range(x.current)

    ## Get the design matrix for the j-th spline (with the basis of above)
    m[[j]] <- eval.basis(x.current, basis[[j]])

    ## Get the degrees of freedom for the current predictor
    df[j] <- NCOL(m[[j]])

    colnames(m[[j]]) <- paste(colnames(x)[j], ".", 1:df[j], sep = "")

    ## Get penalty matrix ## !!! DOUBLE-CHECK THIS !!!
    pen[[j]] <- getbasispenalty(basis[[j]], Lfdobj = 2)
  }
  list(m = m, pen = pen, basis = basis, index = rep(1:p, times = df))
}

.onAttach <- function(libname, pkgname){
  cat("----------------------------------------------------------------------",
      "Please note that this is an early test release of package 'penGAM'.",
      "It should only be used for experimental reasons. Use at your own risk!",
      "----------------------------------------------------------------------",
      sep = "\n")
}
