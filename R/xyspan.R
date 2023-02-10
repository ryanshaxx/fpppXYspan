#' @title xyspan
#'
#' @description Provides an overview table for the time and scope conditions of
#'     a data set
#'
#' @param dat A data set object
#' @param id Scope (e.g., country codes or individual IDs)
#' @param time Time (e.g., time periods are given by years, months, ...)
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(CanadianWeather)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
#' @importFrom dplyr "%>%"

goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=161, norder=8)
LfdobjNonDur <- int2Lfd(4)
argvals = seq(1919,2000,len=length(nondurables))
logNondurSm <- smooth.basisPar(argvals,
                               y=log10(nondurables), fdobj=goodsbasis,
                               Lfdobj=LfdobjNonDur, lambda=1e-11)
phaseplanePlot(1964, logNondurSm$fd)


#find x and y span of phaseplane plots, takes the arguement of a smooth fd object
xyspan <- function (smooth_fd_object) {

  #range for evalarg
  evalarg_rng <- c(smooth_fd_object$basis$rangeval[1]:smooth_fd_object$basis$rangeval[2])
  #number of years
  n_col <- dim(smooth_fd_object$coefs)[2]

  #empty vectors for x and y span
  y_span_vec <- c()
  x_span_vec <- c()

  for (j in 1:n_col){
    temp <- phaseplanePlot(evalarg = evalarg_rng, fdobj =  smooth_fd_object[j], returnMatrix = TRUE)

    #get y span (acceleration)
    convert_diff <-diff(sign(temp[,1]))
    x_intercept_indices <- which(convert_diff != 0)
    #use these indices to get the values of y (acceleration) at these indices
    y_vals <- temp[,2][x_intercept_indices]
    #get min and max of this list to get y span
    y_vals_min <- min(y_vals)
    y_vals_max <- max(y_vals)
    y_span <- y_vals_max - y_vals_min

    #append these values to a vector
    y_span_vec <- c(y_span_vec, y_span)

    #save the points taken for y span
    x_point_min <- temp[,1][which(temp[,2] %in% y_vals_min)]
    x_point_max <- temp[,1][which(temp[,2] %in% y_vals_max)]
    # print(x_point_max)
    # print(y_vals_max)

    #get x span (velocity)
    convert_diff <-diff(sign(temp[,2]))
    y_intercept_indices <- which(convert_diff != 0)
    #use these indices to get the values of x (velocity) at these indices
    x_vals <- temp[,1][y_intercept_indices]
    #get min and max of this list to get y span
    x_vals_min <- min(x_vals)
    x_vals_max <- max(x_vals)
    x_span <- x_vals_max - x_vals_min

    #append these values to a vector
    x_span_vec <- c(x_span_vec, x_span)

    #save the points taken for x span
    y_point_min <- temp[,2][which(temp[,1] %in% x_vals_min)]
    y_point_max <- temp[,2][which(temp[,1] %in% x_vals_max)]
    # print(x_vals_min)
    # print(y_point_min)


  }

  df_span<-data.frame(x_span_vec, y_span_vec)

  #add the years column
  year_col <- smooth_fd_object$fdnames[[2]]

  df_span <- data.frame(year_col, df_span)
  df_span$year_col <- as.integer(df_span$year_col)

  return(df_span)

}
