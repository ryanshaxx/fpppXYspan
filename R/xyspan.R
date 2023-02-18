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
yoyo  =phaseplanePlot(c(1964,1965), logNondurSm$fd, returnMatrix = T)

log_fd = logNondurSm$fd


#find x and y span of phaseplane plots, takes the arguement of a smooth fd object
xyspan <- function (eval_period, smooth_fd_object) {

  if (is.fd(smooth_fd_object)){

    if (is.numeric(eval_period)) {

      if (length(eval_period) == 1) {
        if (eval_period[1] == round(eval_period[1])) {
          arg_vals = seq(eval_period[1], eval_period[1]+1, length=181)
          period = eval_period[1]
        } else {
          print("Error: Please provide a valid wholenumber numeric value")
          return(NULL)
        }
      } else if (length(eval_period) == 2) {
        if ((eval_period[1] == round(eval_period[1])) & (eval_period[2] == round(eval_period[2]))) {
          arg_vals = seq(eval_period[1], eval_period[2], length=181)
          period = paste(eval_period[1], eval_period[2], sep = "-")
        } else {
          print("Error: Please provide a valid wholenumber numeric value(s)")
          return(NULL)
        }

      } else {
        print("Error: Please provide a valid numeric range")
        return(NULL)
      }

      #number of sets of data
      n_col <- dim(smooth_fd_object$coefs)[2]

      #empty vectors for x and y span
      y_span_vec <- c()
      x_span_vec <- c()
      coef_set <- c()

      for (j in 1:n_col){

        temp <- phaseplanePlot(evalarg = arg_vals, fdobj =  smooth_fd_object[j], returnMatrix = TRUE)

        #get y span (acceleration)
        convert_diff <-diff(sign(temp[,1]))
        x_intercept_indices <- which(convert_diff != 0)

        if (length(x_intercept_indices) < 2) {
          y_span = Inf
        } else {
          #use these indices to get the values of y (acceleration) at these indices
          y_vals <- temp[,2][x_intercept_indices]
          #get min and max of this list to get y span
          y_vals_min <- min(y_vals)
          y_vals_max <- max(y_vals)
          y_span <- y_vals_max - y_vals_min
        }


        #append these values to a vector
        y_span_vec <- c(y_span_vec, y_span)

        #get x span (velocity)
        convert_diff <-diff(sign(temp[,2]))
        y_intercept_indices <- which(convert_diff != 0)
        if (length(y_intercept_indices) < 2) {
          x_span = Inf
        } else {
          #use these indices to get the values of x (velocity) at these indices
          x_vals <- temp[,1][y_intercept_indices]
          #get min and max of this list to get y span
          x_vals_min <- min(x_vals)
          x_vals_max <- max(x_vals)
          x_span <- x_vals_max - x_vals_min
        }

        #append these values to a vector
        x_span_vec <- c(x_span_vec, x_span)

        #append which set of coefs this is for
        coef_set <- c(coef_set, j)

        }

      df_span<-data.frame(x_span_vec, y_span_vec)

      #add the range of arguments used in the evaluation column
      evalarg <- toString(period)

      df_span <- data.frame(coef_set, evalarg, df_span)
      #df_span$year_col <- as.integer(df_span$year_col)

      return(df_span)

    } else {
      print("Error: The first argument must be of class 'numeric'")
      return(NULL)
    }

  } else {
    print("Error: Second argument must be a functional data object of class 'fd'")
    return(NULL)
  }

}

check = xyspan(c(1920,2000),  log_fd)

