show_rank_survey = function(object, left_y, right_y) {
  stopifnot(class(object) == "Survey" | is.data.frame(object))
}
