#' Run Euler's method on an ODE using data.table for improved performance.
#'
#' @param dydt The definition of the ODE. This may be just the right-hand side
#' of the ODE, or you may specify an equation using _t to denote derivative. See examples.
#' @param tlim a list containing the start and end of the timeframe to integrate over
#' @param ic the initial condition. If this is a scalar ODE then you can do
#' ic=4. If this is a system then you must provide a list of values in the correct order as determined by the
#' listing of variables in dydt. You may also provide a list with names, in that
#' case order doesn't matter.
#' @param stepSize length of timestep to take.
#' @examples
#' # See original function documentation for examples
#' @returns A data.table containing the approximations of the solutions and the
#' corresponding times.
#' @export
Euler_dt = function(dydt, tlim, ic, stepSize = (tlim[[2]] - tlim[[1]]) / 10, ...) {
  # Check if data.table package is available
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package is required. Please install it with install.packages('data.table')")
  }
  
  t0 = tlim[[1]]
  t_final = tlim[[2]]
  
  # Handle dots parameters
  dots = list(...)
  if (length(dots) > 0) {
    # Check to see if this is specifying a derivative
    derivstr = names(dots)[1]
    if (substr(derivstr, nchar(derivstr) - 1, nchar(derivstr)) == "_t") {
      # Here is the formula for the derivative
      dform = deparse(mosaic::lhs(dots[[derivstr]]))
      # Here are the variables
      dformvars = deparse(mosaic::rhs(dots[[derivstr]]))
      # Put this in the form c(y_t=t*y)~y
      dydt = stats::as.formula(paste("c(", derivstr, "=", dform, ")~", dformvars))
    }
  }
  
  # Parse the formula to get state variables
  dydtstr = as.character(dydt)
  
  # Get the variable declarations
  RHS = dydtstr[[3]]
  
  # Split at the & and delete whitespace
  allInputVars = strsplit(gsub("[[:space:]]", "", RHS), split = "&")[[1]]
  
  # Make a list of the state variables
  stateVarsInputOrder = c()
  isAutonomous = TRUE
  for (i in 1:length(allInputVars)) {
    if (allInputVars[[i]] != "t") {
      stateVarsInputOrder = c(stateVarsInputOrder, allInputVars[[i]])
    } else {
      isAutonomous = FALSE
    }
  }
  
  # If the function was declared without t as a variable, add t to the list
  if (isAutonomous) {
    dydt = stats::as.formula(paste(dydtstr[[2]], "~t&", dydtstr[[3]], collapse = " "))
  }
  
  # Process initial conditions
  if (length(ic) != length(stateVarsInputOrder)) {
    cat(paste("Error. You need to supply ", length(stateVarsInputOrder), 
              " initial conditions. You supplied ", length(ic), "."))
    return()
  }
  
  # If no names given in IC, assume IC is in the "correct" order
  if (is.null(names(ic))) {
    names(ic) = stateVarsInputOrder
  }
  
  # Check that the names given in the initial condition match the state variables
  if (length(setdiff(names(ic), stateVarsInputOrder)) > 0 ||
      length(setdiff(stateVarsInputOrder, names(ic))) > 0) {
    cat(paste("Error. You must give initial conditions with names ", 
              paste(stateVarsInputOrder, collapse = ","), " you gave initial conditions with names ",
              paste(names(ic), collapse = ",")))
    return()
  }
  
  # Make sure the initial condition is in the correct order
  y0 = ic[stateVarsInputOrder]
  
  # Make a function out of dydt
  dydtfunc = mosaicCore::makeFun(dydt)
  
  # Make a list of what all the derivatives should be
  derivnames = paste0(stateVarsInputOrder, "_t")
  
  # Test the function at the initial condition
  dydtvalue = do.call(dydtfunc, args = as.list(c(t = t0, y0)))
  
  # Is it returning the correct number of outputs?
  if (length(dydtvalue) != length(stateVarsInputOrder)) {
    cat(paste("Error. You supplied ", length(stateVarsInputOrder), 
              " state variables, but gave ", length(dydtvalue), "differential equations."))
    return()
  }
  
  # Ensure proper naming of outputs
  basedydtfunc = dydtfunc
  if (length(names(dydtvalue)) == 0) {
    dydtfunc = function(...) {
      output = basedydtfunc(...)
      names(output) = derivnames
      return(output)
    }
  }
  
  dydtvalue = do.call(dydtfunc, args = as.list(c(t = t0, y0)))
  
  # Check if derivative names are correct
  if (length(setdiff(names(dydtvalue), derivnames)) != 0 ||
      length(setdiff(derivnames, names(dydtvalue))) != 0) {
    cat(paste("Error, your ODEs should have ", paste(derivnames, collapse = ","), 
              " on the left-hand side. Instead they have ", paste(names(dydtvalue), collapse = ",")))
    return()
  }
  
  # Ensure outputs are in the same order as inputs
  basedydtfunc2 = dydtfunc
  if (!all(names(dydtvalue) == derivnames)) {
    dydtfunc = function(...) {
      output = basedydtfunc2(...)
      output = output[derivnames]
      return(output)
    }
  }
  
  # Use data.table for solution storage
  numSteps = ceiling((t_final - t0) / stepSize)
  
  # Pre-allocate the data.table
  solution = data.table::data.table(t = numeric(numSteps + 1))
  for (var in stateVarsInputOrder) {
    solution[, (var) := numeric(numSteps + 1)]
  }
  
  # Set the first row to initial conditions
  solution[1, "t" := t0]
  for (var in stateVarsInputOrder) {
    solution[1, (var) := y0[var]]
  }
  
  # Create a list to hold values during each iteration (for performance)
  yi = y0
  ti = t0
  
  # Use data.table's set() for fast row updates during iteration
  for (i in 2:(numSteps + 1)) {
    yi = yi + stepSize * do.call(dydtfunc, args = as.list(c(t = ti, yi)))
    ti = ti + stepSize
    
    # Update the row using data.table's set() - much faster than [i,]
    data.table::set(solution, i, "t", ti)
    for (j in seq_along(stateVarsInputOrder)) {
      data.table::set(solution, i, stateVarsInputOrder[j], yi[j])
    }
  }
  
  # Return the data.table
  return(solution)
}