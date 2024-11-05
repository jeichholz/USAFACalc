#' Plot a function of one or two variables
#' @description
#' USAFACalc::plotFun is just a few bugfixes on mosaic::plotFun.  It also defaults npts to
#' 200 for contour plots.  The bug fixes allow you to set line color, width, and type on contour plots.
#'
#' @inheritParams mosaic::plotFun
#' @export
plotFun<-function (object, ..., plot = lattice::trellis.last.object(), add = NULL,
          under = FALSE, xlim = NULL, ylim = NULL, npts = NULL, ylab = NULL,
          xlab = NULL, zlab = NULL, filled = TRUE, levels = NULL,
          nlevels = 10, labels = TRUE, surface = FALSE, groups = NULL,
          col = lattice::trellis.par.get("superpose.line")$col, col.regions = topo.colors,
          type = "l", lwd = lattice::trellis.par.get("superpose.line")$lwd,
          lty = lattice::trellis.par.get("superpose.line")$lty, alpha = NULL,
          discontinuities = NULL, discontinuity = 1, interactive = mosaic::rstudio_is_available())
{
  if (is.function(object)) {
    formula <- f(x) ~ x
    formula[[2]] <- as.call(list(substitute(object), quote(x)))
    object <- formula
  }
  if (is.null(add))
    add <- under
  if (discontinuity <= 0 | length(discontinuity) != 1)
    stop("discontinuity should be a positive number or Inf.")
  dots <- list(...)
  ..f.. <- do.call(mosaicCore::makeFun, c(object, dots, strict.declaration = FALSE),
                   envir = parent.frame())
  if (is.vector(col.regions))
    col.regions <- makeColorscheme(col.regions)
  if (length(unique(levels)) == 1)
    levels = c(levels, Inf)
  vars <- formals(..f..)
  rhsVars <- all.vars(mosaicCore::rhs(object))
  otherVars <- setdiff(names(vars), rhsVars)
  ndims <- length(rhsVars)
  cleanDots <- dots
  for (v in otherVars) {
    cleanDots[[v]] <- NULL
  }
  pgArgs <- list()
  otherGroups <- if (length(otherVars) > 0)
    paste(otherVars, "", sep = "")
  else c()
  if (length(setdiff(otherGroups, names(dots))) > 0) {
    stop(paste("Cannot plot with free parameters; try setting",
               paste(setdiff(otherGroups, names(dots)), collapse = ", ")))
  }
  for (param in otherVars) pgArgs[[param]] <- dots[[paste(param,
                                                          "", sep = "")]]
  paramGrid <- do.call("expand.grid", pgArgs)
  fList <- if (nrow(paramGrid) < 1)
    list(..f..)
  else list()
  for (r in mosaicCore::rows(paramGrid)) {
    rowAsList <- as.list(paramGrid[r, ])
    names(rowAsList) <- otherVars
    fList <- c(fList, do.call(mosaicCore::makeFun, c(object, c(dots,
                                                     rowAsList, strict.declaration = FALSE, envir = parent.frame()))))
  }
  if (add) {
    if (ndims == 1) {
      rlang::check_installed("latticeExtra")
      return(plot + latticeExtra::layer(do.call(panel.usafacalc.plotFun1,
                                                c(list(..f.. = fList, npts = npts, discontinuity = discontinuity,
                                                       discontinuities = discontinuities, filled = filled,
                                                       levels = levels, nlevels = nlevels, surface = surface,
                                                       col.regions = col.regions, type = type, alpha = alpha,
                                                       col = col, lty = lty, lwd = lwd), dots)),
                                        data = as.list(environment()), ..., under = under))
    }
    rlang::check_installed("latticeExtra")
    return(plot + latticeExtra::layer(do.call(panel.usafacalc.plotFun,
                                              c(list(object = object, npts = npts, discontinuity = discontinuity,
                                                     discontinuities = discontinuities, filled = filled,
                                                     levels = levels, nlevels = nlevels, suface = surface,
                                                     col.regions = col.regions, type = type, alpha = alpha,
                                                     lwd = lwd, lty = lty,col=col), dots)), data = as.list(environment()),
                                      under = under))
  }
  limits <- inferArgs(dots = dots, vars = rhsVars, defaults = list(xlim = xlim,
                                                                   ylim = ylim))
  if (ndims == 1) {
    npts <- ifelse(is.null(npts), 5000, npts)
    if (is.null(ylab))
      ylab <- deparse(mosaicCore::lhs(object))
    if (is.null(xlab))
      xlab <- rhsVars
    if (is.null(limits$xlim) || length(limits$xlim) < 2) {
      zeros <- c()
      tryCatch(zeros <- mosaic::findZeros(object, nearest = 6,
                                  ...)[[1]], error = function(e) {
                                    e
                                  }, warning = function(e) {
                                  })
      limits$xlim <- switch(as.character(length(zeros)),
                            `0` = c(0, 1), `1` = c(-1.5, 1.5) * (zeros +
                                                                   ifelse(zeros == 0, 1, 0)), c(-0.1, 0.1) *
                              diff(range(zeros)) + range(zeros))
    }
    limits$xlim <- range(limits$xlim)
    if (length(limits$xlim) != 2)
      stop("Invalid limits.")
    .f. <- fList[[1]]
    .xvals <- if ("h" %in% type)
      seq(base::min(limits$xlim), base::max(limits$xlim),
          length.out = npts)
    else mosaic::adapt_seq(base::min(limits$xlim), base::max(limits$xlim),
                   f = function(xxqq) {
                     .f.(xxqq)
                   }, length.out = npts, quiet = TRUE)
    .yvals <- c()
    for (.f. in fList) .yvals <- c(.yvals, sapply(.xvals,
                                                  .f.))
    localData <- data.frame(x = .xvals, y = .yvals, component = contintuous_components(.xvals,
                                                                                       .yvals, adjust = discontinuity))
    localData$alone <- mosaicCore::ediff(localData$component) * mosaicCore::ediff(localData$component,
                                                          pad = "tail") != 0
    localData <- localData[!localData$alone, ]
    localData <- data.frame(x = range(localData$x, na.rm = TRUE,
                                      fintie = TRUE), y = range(localData$y, na.rm = TRUE,
                                                                finite = TRUE))
    if (length(limits$ylim) != 2) {
      thePlot <- do.call(lattice::xyplot, c(list(y ~ x,
                                                 ..f.. = fList, data = localData, groups = substitute(groups),
                                                 xlim = limits$xlim, xlab = xlab, ylab = ylab,
                                                 panel = "panel.usafacalc.plotFun1", npts = npts, discontinuity = discontinuity,
                                                 discontinuities = discontinuities, col = col,
                                                 lwd = lwd, lty = lty), cleanDots))
    }
    else {
      thePlot <- do.call(lattice::xyplot, c(list(y ~ x,
                                                 ..f.. = fList, data = localData, groups = substitute(groups),
                                                 xlim = limits$xlim, ylim = limits$ylim, xlab = xlab,
                                                 ylab = ylab, panel = "panel.usafacalc.plotFun1", npts = npts,
                                                 discontinuity = discontinuity, discontinuities = discontinuities,
                                                 col = col, lty = lty, lwd = lwd, alpha = alpha),
                                            cleanDots))
    }
    return((thePlot))
  }
  if (ndims == 2) {

    if (surface){
      npts <- ifelse(is.null(npts), 40, npts)
      alpha=1
    }
    else{
      npts <- ifelse(is.null(npts), 200, npts)
    }
    if (length(ylab) == 0)
      ylab <- rhsVars[2]
    if (length(xlab) == 0)
      xlab <- rhsVars[1]
    if (length(zlab) == 0)
      zlab <- deparse(mosaicCore::lhs(object))
    if (is.null(limits$xlim) || length(limits$xlim) < 2)
      limits$xlim <- c(0, 1)
    else limits$xlim <- range(limits$xlim)
    if (is.null(limits$ylim) || length(limits$ylim) < 2)
      limits$ylim <- c(0, 1)
    else limits$ylim <- range(limits$ylim)
    .xvals <- seq(base::min(limits$xlim), base::max(limits$xlim),
                  length = npts)
    .yvals <- seq(base::min(limits$ylim), base::max(limits$ylim),
                  length = npts)
    zvals <- tryCatch(outer(.xvals, .yvals, function(x,
                                                     y) {
      ..f..(x, y)
    }), warning = function(w) {
    })
    grid <- expand.grid(.xvals, .yvals)
    grid$height <- c(zvals)
    localData <- grid
    names(localData) <- c(rhsVars, ".height")
    if (surface) {
      if (add) {
        stop("Should not get here, but no add option for surface plots yet anyway.")
        return(invisible(NULL))
      }
      zcuts = pretty(grid$height, 50)

      zcolors = col.regions(length(zcuts), alpha =1)
      if (interactive && mosaic::rstudio_is_available()) {

        #Keep the original behavior, just in case.
        fig0=manipulate::manipulate(lattice::wireframe(height ~
                                                         Var1 + Var2, xlab = xlab, ylab = ylab, zlab = list(zlab,
                                                                                                            rot = 90), data = grid, groups = eval(substitute(groups),
                                                                                                                                                  localData), drape = filled, shade = FALSE,
                                                       colorkey = FALSE, scales = list(arrows = FALSE),
                                                       screen = c(z = rot, x = elev - 90), distance = dist,
                                                       at = zcuts, col.regions = zcolors), rot = manipulate::slider(min = -180,
                                                                                                                    max = 180, step = 5, initial = 35, label = "Rotation"),
                                    elev = manipulate::slider(min = -90, max = 90,
                                                              step = 5, initial = 30, label = "Elevation"),
                                    dist = manipulate::slider(min = 0, max = 1,
                                                              step = 0.01, initial = 0.2, label = "Distance"))
        #Create a plotly surface plot and print it, it looks better. By doing it last, if this is working, then the Viewer tab should be
        #the one that the student sees by default.
        fig <- plotly::plot_ly(x=.xvals,y=.yvals,z=zvals,
                               contours = list(
                                 z = list(
                                   show=TRUE,
                                   usecolormap=TRUE,
                                   highlightcolor="black",
                                   project=list(z=TRUE)
                                 )
                               )) %>% plotly::add_surface()
        print(fig)


        #return(fig)
        #return the original, just in case that is important in some way.
        return(fig0);
      }
      else {
        return((lattice::wireframe(height ~ Var1 + Var2, xlab = xlab,
                          ylab = ylab, zlab = list(zlab, rot = 90),
                          data = grid, groups = eval(substitute(groups),
                                                     localData), drape = filled, shade = FALSE,
                          colorkey = FALSE, col.regions = zcolors, at = zcuts,
                          ...)))
      }
    }
    else {
      funPlot.draw.contour <- function(x, y, z, ncontours = 6,
                                       at = NULL, filled = TRUE, col.regions = topo.colors,
                                       labels = TRUE, showlabels = TRUE, contours = TRUE,
                                       groups = NULL, label = TRUE, xlab = "", ylab = "",
                                       ...) {
        if (is.null(at))
          at = pretty(z, ncontours)
        argsToPass <- list(z ~ x * y, at = at, xlab = xlab,
                           ylab = ylab, panel = panel.usafacalc.levelcontourplot,
                           groups = substitute(groups), col.regions = col.regions(60),
                           contour = contours, labels = labels, colorkey = FALSE,
                           retion = TRUE, filled = filled, ...)
        argsToPass[["k"]] <- NULL
        return((do.call(lattice::levelplot, argsToPass)))
      }
      if (is.null(alpha))
        alpha <- 1
      fillcolors <- col.regions(length(levels) + 2, alpha = alpha)
      if (all(is.logical(zvals))) {
        fillcolors <- col.regions(4, alpha = alpha)
        nlevels <- 2
      }
      if (add) {
        return(mosaic::ladd(panel.usafacalc.levelcontourplot(grid$Var1,
                                           grid$Var2, grid$height, subscripts = 1, at = levels,
                                           labels = labels, filled = filled, groups = eval(substitute(groups),
                                                                                           localData), col.regions = col.regions, contour = TRUE,
                                           region = FALSE, lwd=lwd, col=col,lty=lty, ...)))
      }
      return((funPlot.draw.contour(grid$Var1, grid$Var2,
                                   grid$height, xlab = xlab, ylab = ylab, at = levels,
                                   filled = filled, labels = labels, groups = eval(substitute(groups),
                                                                                   localData), col.regions = col.regions, lwd=lwd,col=col, lty=lty, ...)))
    }
  }
  stop("Bug alert: You should not get here.  Please report.")
}

#' @export
guess_discontinuities <- function( f, xlim, ..., resolution = 10000, adjust = 1) {
  x <- seq(xlim[1], xlim[2], length.out = resolution)
  F <- Vectorize(f)
  y <- F(x, ...)
  discontinuity_at(x, y, adjust = adjust)
}

#' @export
discontinuity_at <- function(x, y, adjust = 1) {
  y[!is.finite(y)] <- NA
  y_scaled <- (y - base::min(y, na.rm = TRUE)) / diff(range(y[is.finite(y)]))
  x_scaled <- (x - base::min(x, na.rm = TRUE)) / diff(range(x[is.finite(x)]))
  slope  <- diff(y_scaled) / diff(x_scaled)
  jump <- diff(y_scaled)
  left_jump <- c(0, jump)
  right_jump <- c(jump, 0)
  left_slope  <- c(slope[1], slope)
  right_slope  <- c(slope, slope[length(slope)])
  x[which( is.na(y) |
             !is.finite(y) |
             (
               (abs(atan(left_slope) - atan(right_slope)) > adjust / 2) &
                 (abs(atan(left_slope)) > adjust/2 | abs(atan(right_slope)) > adjust/2)
             )
  )]
}

#' @export
contintuous_components <- function(x, y, adjust = 1) {
  y[!is.finite(y)] <- NA
  y_scaled <- (y - base::min(y, na.rm = TRUE)) / diff(range(y[is.finite(y)]))
  x_scaled <- (x - base::min(x, na.rm = TRUE)) / diff(range(x[is.finite(x)]))
  slope  <- diff(y_scaled) / diff(x_scaled)
  jump <- diff(y_scaled)
  left_jump <- c(0, jump)
  right_jump <- c(jump, 0)
  left_slope  <- c(slope[1], slope)
  right_slope  <- c(slope, slope[length(slope)])
  discontinuities <- x[which( is.na(y) |
                                !is.finite(y) |
                                (
                                  (abs(atan(left_slope) - atan(right_slope)) > adjust /2) &
                                    (abs(atan(left_slope)) > adjust /2 | abs(atan(right_slope)) > adjust /2)
                                )
  )]
  j <- sapply(x, function(x) sum(x > discontinuities))
  k <- sapply(x, function(x) sum(x >= discontinuities))
  j + k
}

# lengths of points on one continuous branch of a function.
#' @export
branch_lengths <- function(x, y, discontinuities = NULL, discontinuity = 1) {
  if (is.null(discontinuities)) discontinuities <- discontinuity_at(x, y, adjust = discontinuity)
  if (length(discontinuities) < 1L) return( length(x) )
  # check <  and <= in case a grid point is exactly a discontinuity.
  j <- sapply(x, function(x) sum(x > discontinuities))
  k <- sapply(x, function(x) sum(x >= discontinuities))
  diff( c( 0, which(diff(j+k) > 0), length(x) ) )
}


#' @export
panel.usafacalc.plotFun1 <- function( ..f.., ...,
                            x, y,
                            type="l",
                            lwd = lattice::trellis.par.get("superpose.line")$lwd,
                            lty = lattice::trellis.par.get("superpose.line")$lty,
                            col = lattice::trellis.par.get('superpose.line')$col,
                            npts=NULL,
                            zlab=NULL,
                            filled=TRUE,
                            levels=NULL,
                            nlevels=10,
                            surface=FALSE,
                            alpha=NULL,
                            discontinuity = NULL,
                            discontinuities = NULL) {
  dots <- list(...)

  if (is.function(..f..) ) ..f.. <- list(..f..)
  if (! is.list(..f..) || length(..f..) < 1) {
    print(str(..f..))
    stop("Empty or malformed list of functions.")
  }

  if (is.numeric(col)) {
    message('converting numerical color value into a color using lattice settings')
    col <- lattice::trellis.par.get('superpose.line')$col[col]
  }


  # funny names (like ..f..) are to avoid names that might be used by the user
  # not sure whether this precaution is necessary in current implementation

  # perhaps use environment(object)?
  # ..f.. <- do.call( "makeFun", c(object, dots, strict.declaration=FALSE), envir= parent.frame())
  idx <- 0

  for ( .f. in ..f..) {
    idx <- idx + 1
    # print(.f.)
    vars <- formals(.f.)
    #  rhsVars <- all.vars(rhs(object))
    #  ndims <- length(rhsVars)

    parent.xlim <- lattice::current.panel.limits()$xlim
    parent.ylim <- lattice::current.panel.limits()$ylim

    npts <- ifelse( is.null(npts), 200, npts)

    if (! missing(x) && !missing(y) && length(x) >= npts && length(y) == length(x)) {
      .xvals <- x
      .yvals <- y
    } else {
      # Evaluate the function on appropriate inputs to help figure out y limits.
      .xvals <-  if ('h' %in% type)
        seq(base::min(parent.xlim), base::max(parent.xlim), length.out=npts)
      else
        mosaic::adapt_seq(base::min(parent.xlim), base::max(parent.xlim),
                  f=function(xxqq){ .f.(xxqq) },
                  length.out=npts,
                  quiet=TRUE)

      .yvals <- suppressWarnings( sapply( .xvals, .f. ) )
    }

    # need to strip out any components of ... that are in the object so they
    # don't get passed to the panel function.
    cleandots = list(...)
    # cleandots[ intersect(names(cleandots), all.vars(object)) ] <- NULL
    cleandots[c('x','y','type','alpha','col')] <- NULL
    # use do.call to call the panel function so that the cleandots can be put back in

    if (type == "l") {
      if (is.null(discontinuities))
        discontinuities <- discontinuity_at(.xvals, .yvals, adjust = discontinuity)
      do.call(
        grid::grid.polyline,
        c(list(x=.xvals, y=.yvals, default.units = "native",
               gp=do.call(
                 grid::gpar,
                 c(list(alpha=alpha, col=.getColor(idx, col),
                        lty = lty, lwd = lwd),
                   cleandots)
               ),
               id.lengths = branch_lengths(.xvals, .yvals, discontinuities)
        ))
      )
    } else {
      do.call(
        lattice::panel.xyplot,
        c(list(x=.xvals, y=.yvals, type=type,  alpha=alpha,
               lwd = lwd, lty = lty, col=.getColor(idx,col)),  cleandots)
      )
    }
  }
}



#' @export
panel.usafacalc.plotFun1a <- function( ..f.., ...,
                             x, y,
                             type="l",
                             col = lattice::trellis.par.get('superpose.line')$col,
                             lwd = lattice::trellis.par.get("superpose.line")$lwd,
                             lty = lattice::trellis.par.get("superpose.line")$lty,
                             npts=NULL,
                             zlab=NULL,
                             filled=TRUE,
                             levels=NULL,
                             nlevels=10,
                             surface=FALSE,
                             alpha=NULL,
                             discontinuity = NULL,
                             discontinuities = NULL) {
  dots <- list(...)

  if (is.function(..f..) ) ..f.. <- list(..f..)
  if (! is.list(..f..) || length(..f..) < 1) {
    print(str(..f..))
    stop("Empty or malformed list of functions.")
  }

  if (is.numeric(col)) {
    message('converting numerical color value into a color using lattice settings')
    col <- lattice::trellis.par.get('superpose.line')$col[col]
  }


  # funny names (like ..f..) are to avoid names that might be used by the user
  # not sure whether this precaution is necessary in current implementation

  # perhaps use environment(object)?
  # ..f.. <- do.call( "makeFun", c(object, dots, strict.declaration=FALSE), envir= parent.frame())

  parent.xlim <- lattice::current.panel.limits()$xlim
  parent.ylim <- lattice::current.panel.limits()$ylim
  points_data <- data.frame()
  groups_data <- data.frame()
  id_lengths <- integer(0)
  idx <- 0L

  for (.f. in ..f..) {
    idx <- idx + 1
    # print(.f.)
    vars <- formals(.f.)
    #  rhsVars <- all.vars(rhs(object))
    #  ndims <- length(rhsVars)

    npts <- ifelse( is.null(npts), 200, npts)

    if (! missing(x) && !missing(y) && length(x) >= npts && length(y) == length(x)) {
      .xvals <- x
      .yvals <- y
    } else {
      # Evaluate the function on appropriate inputs to get y values.
      .xvals <-  if ('h' %in% type)
        seq(base::min(parent.xlim), base::max(parent.xlim), length.out=npts)
      else
        mosaic::adapt_seq(base::min(parent.xlim), base::max(parent.xlim),
                  f=function(xxqq){ .f.(xxqq) },
                  length.out=npts,
                  quiet=TRUE)

      tryCatch(.yvals <- sapply( .xvals, .f. ), warning=function(w) {} )
    }

    if (is.null(discontinuities))
      discontinuities <- discontinuity_at(.xvals, .yvals, adjust = discontinuity)
    groups_data <-
      dplyr::bind_rows(
        groups_data,
        data.frame(
          length = branch_lengths(.xvals, .yvals, discontinuities),
          id = idx
        )
      )

    points_data <-
      dplyr::bind_rows(points_data,
                       data.frame(x = .xvals, y = .yvals, id=idx))

  }

  # points_data is now a data frame with x, y, and an id showing which function
  # is being described.

  # need to strip out any components of ... that are in the object so they
  # don't get passed to the panel function.
  cleandots = list(...)
  # cleandots[ intersect(names(cleandots), all.vars(object)) ] <- NULL
  cleandots[c('x','y','type','alpha','col')] <- NULL
  # use do.call to call the panel function so that the cleandots can be put back in

  if (type == "l") {
    do.call(
      grid::grid.polyline,
      c(list(x=points_data$x, y=points_data$y, default.units = "native",
             gp=do.call(
               grid::gpar,
               c(list(alpha=alpha, col=.getColor(groups_data$id, col),
                      lwd = lwd, lty = lty), cleandots)
             ),
             id.lengths = groups_data$length
      ))
    )
  } else {
    do.call(
      lattice::panel.xyplot,
      c(list(x=points_data$x, y=points_data$y, type=type,
             lwd = lwd, lty = lty,
             alpha=alpha, col=.getColor(points_data$id, col)),  cleandots)
    )
  }
}




#' @export
panel.usafacalc.plotFun <- function( object, ...,
                           type="l",
                           npts=NULL,
                           zlab=NULL,
                           filled=TRUE,
                           levels=NULL,
                           nlevels=10,
                           surface=FALSE,
                           col.regions =topo.colors,
                           lwd = lattice::trellis.par.get("superpose.line")$lwd,
                           lty = lattice::trellis.par.get("superpose.line")$lty,
                           col = lattice::trellis.par.get("superpose.line")$col,
                           alpha=NULL,
                           discontinuity = NULL,
                           discontinuities = NULL ) {
  dots <- list(...)
  if ( is.function(object) ) {
    formula <- f(x) ~ x
    formula[[2]] <- as.call( list(substitute(object), quote(x)))
    object <- formula
  }

  if ( is.vector(col.regions ) ) col.regions  <- makeColorscheme(col.regions )

  plot.line <- lattice::trellis.par.get('plot.line')
  superpose.line <- lattice::trellis.par.get('superpose.line')

  # funny names (like ..f..) are to avoid names that might be used by the user
  # not sure whether this precaution is necessary in current implementation

  # perhaps use environment(object)?
  ..f.. <- do.call( mosaicCore::makeFun, c(object, dots, strict.declaration=FALSE), envir= parent.frame())

  vars <- formals(..f..)
  rhsVars <- all.vars(mosaicCore::rhs(object))
  ndims <- length(rhsVars)

  parent.xlim <- lattice::current.panel.limits()$xlim
  parent.ylim <- lattice::current.panel.limits()$ylim

  if( ndims > 2 || ndims < 1 )
    stop("Formula must provide 1 or 2 independent variables (right hand side).")

  if( ndims == 1 ){
    npts <- ifelse( is.null(npts), 200, npts)


    # Evaluate the function on appropriate inputs.
    .xvals <-
      if ('h' %in% type)
        seq(base::min(parent.xlim), base::max(parent.xlim), length.out=npts)
    else
      mosaic::adapt_seq(base::min(parent.xlim), base::max(parent.xlim),
                f=function(xxqq){ ..f..(xxqq) },
                length.out=npts,
                quiet=TRUE)
    .yvals <- tryCatch( sapply( .xvals, ..f.. ),
                        warning=function(w) {} )


    # need to strip out any components of ... that are in the object so they
    # don't get passed to the panel function.
    cleandots = list(...)
    cleandots[ names(cleandots) %in% all.vars(object) ] <- NULL
    # use do.call to call the panel function so that the cleandots can be put back in
    if (type == "l") {
      idx <- 1
      return(
        do.call(
          grid::grid.polyline,
          c(list(x=.xvals, y=.yvals, default.units = "native",
                 gp=do.call(
                   grid::gpar,
                   c(list(alpha=alpha, col=.getColor(idx, col),
                          lty = lty, lwd = lwd), cleandots)
                 ),
                 id.lengths = branch_lengths(.xvals, .yvals, discontinuities)
          ))
        )
      )
    } else {
      return(do.call(
        lattice::panel.xyplot,
        c(list(x=.xvals, y=.yvals, type=type, alpha=alpha, lty = lty, lwd = lwd),
          cleandots)))
    }
  }

  if (ndims == 2 ) {
    if( surface ) {
      stop('no add option for surface plots yet.')
      return(NULL)
    }
    # if we get here, surface == FALSE & ndims=2
    npts <- ifelse( is.null(npts), 200, npts)
    # create a function of those two variables

    if( length(zlab) == 0 ) zlab <- deparse(mosaicCore::lhs(object) )

    .xvals <- seq(base::min(parent.xlim),base::max(parent.xlim),length=npts)
    .yvals <- seq(base::min(parent.ylim),base::max(parent.ylim),length=npts)
    zvals <- tryCatch(
      outer(.xvals, .yvals, function(x,y){..f..(x,y)} ),
      warning=function(w) {} )
    grid <- expand.grid( .xvals, .yvals )
    grid$height <- c(zvals)

    zcuts <- pretty(grid$height,50)
    zcolors <- col.regions (length(zcuts),alpha=.5)
    # print(zcolors)
    if( is.null(alpha) ) alpha<-.4

    if( all(is.logical(zvals)) ) {  # it's a constraint function
      #nlevels <- 2
      levels <- c(0.0,Inf)
    }
    fillcolors <- col.regions (length(levels) + 2, alpha=alpha)
    if(is.null(levels)) levels=pretty(grid$height, nlevels)

    return( panel.usafacalc.levelcontourplot(x = grid$Var1, y = grid$Var2, z = grid$height,
                                   subscripts = 1:nrow(grid),
                                   at = levels,
                                   col.regions = fillcolors,
                                   filled=filled, lwd=lwd,lty=lty,col=col,
                                   ...

    )
    )
  }
  stop("Bug alert: You should not get here.  Please report.")
}


#' @export
.getColor <- function( n=1, col=lattice::trellis.par.get('superpose.line')$col) {
  if (length(col) < n) col <- rep(col, length.out=n)
  col[n]
}
# Infer arguments
#'
# The primary purpose is for inferring argument settings from names derived from variables
# occurring in a formula.  For example, the default use is to infer limits for variables
# without having to call them `xlim` and `ylim` when the variables in the formula
# have other names.  Other uses could easily be devised by specifying different `variants`.
#
# @param vars a vector of variable names to look for
# @param dots a named list of argument values
# @param defaults named list or alist of default values for limits
# @param variants a vector of optional postfixes for limit-specifying variable names
# @return a named list or alist of limits.  The names are determined by the names in `defaults`.
#
# If multiple `variants` are matched, the first is used.
# @examples
# inferArgs(c('x','u','t'), list(t=c(1,3), x.lim=c(1,10), u=c(1,3), u.lim=c(2,4)))
# inferArgs(c('x','u'), list(u=c(1,3)), defaults=list(xlim=c(0,1), ylim=NULL))

#' @export
inferArgs <- function( vars, dots, defaults=alist(xlim=, ylim=, zlim=), variants=c('.lim','lim') ) {
  limNames <- names(defaults)
  if (length(vars) > length(limNames))
    stop( paste("Can only infer",
                paste(limNames, sep=""),
                "; but you provided",
                length(vars)),
          " variables to search for.",
          sep="")
  result <-defaults

  for (slot in 1:length(vars)) {
    for (variant in rev(variants)){
      var.lim <- paste(vars[slot],variant,sep="")
      if (! is.null(dots[[var.lim]]) ) {
        result[[limNames[slot]]] <- dots[[var.lim]]
      }
    }
  }
  return(result)
}

#' =============================
#' Create an environment for storing axis limits, etc.
#' .plotFun.envir = new.env(parent=baseenv())
#' =====
#'
#' Lattice plot that draws a filled contour plot
#'
#' Used within plotFun
#'
#' @param x x on a grid
#' @param y y on a grid
#' @param z zvalues for the x and y
#' @param subscripts which points to plot
#' @param at cuts for the contours
#' @param shrink what does this do?
#' @param labels draw the contour labels
#' @param label.style where to put the labels
#' @param contour logical draw the contours
#' @param region logical color the regions
#' @param col color for contours
#' @param lty type for contours
#' @param lwd width for contour
#' @param border type of border
#' @param ... dots additional arguments
#' @param col.regions  a vector of colors or a function (`topo.colors` by default) for generating such
#' @param filled whether to fill the contours with color
#' @param alpha.regions transparency of regions
#' @importFrom lattice trellis.par.get
#' @importFrom lattice panel.levelplot

#' @export
panel.usafacalc.levelcontourplot <- function(x, y, z, subscripts=1,
                                   at, shrink, labels = TRUE,
                                   label.style = c("mixed","flat","align"),
                                   contour = FALSE,
                                   region = TRUE,
                                   col = add.line$col,
                                   lty = add.line$lty,
                                   lwd = add.line$lwd,
                                   border = "transparent", ...,
                                   col.regions = regions$col,
                                   filled=TRUE,
                                   alpha.regions = regions$alpha
)
{
  add.line <- lattice::trellis.par.get('add.line')
  regions <- lattice::trellis.par.get('regions')

  if(filled) panel.levelplot(x, y, z, subscripts,
                             at = pretty(z,5*length(at)), shrink,
                             labels = FALSE,
                             label.style = label.style,
                             contour = FALSE,
                             region = TRUE,
                             border = border, ...,
                             col.regions = col.regions #,
                             #                  alpha.regions = regions$alpha
  )
  if( all(is.logical(z)) ) {  # it's a constraint function
    at <- c(0,1)
    labels <- FALSE
  }
  panel.levelplot(x, y, z, subscripts,
                  at = at, shrink, labels = labels,
                  label.style = label.style,
                  contour = TRUE,
                  region = FALSE, lty=lty,
                  col = col, lwd=lwd,
                  border = border, ...)
}

#' Create a color generating function from a vector of colors
#'
#' @param col a vector of colors
#' @return a function that generates a vector of colors interpolated among the colors in `col`
#'
#' @examples
#' cs <- makeColorscheme( c('red','white','blue') )
#' cs(10)
#' cs(10, alpha=.5)
#'


#' @export
makeColorscheme <- function(col) {
  result <- function(n, alpha=1, ...)  {
    idx <-  0.5 + (0:n)/(n+0.001) * (length(colorList))
    return( apply(col2rgb(col[ round(idx) ]), 2,
                  function(x) { rgb(x[1],x[2],x[3], alpha=round(alpha*255), maxColorValue=255) } ) )
  }
  e <- new.env()
  e[['colorList']] <- col
  environment(result) <- e
  return(result)
}
