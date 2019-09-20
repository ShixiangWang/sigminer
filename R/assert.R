assert_class <- function(object, class) {
  if (inherits(object, what = class)) {
    TRUE
  } else {
    stop("Input should be a ", class, " object!")
  }
}
