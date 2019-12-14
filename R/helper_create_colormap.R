helper_create_colormap <- function(x, n) {
  if (length(x) > 7) {
    stop("Not support so many colors now!")
  }
  # brewer = c("GnBu", "Reds", "Oranges", "Purples", "Greens", "Blues", "Greys")
  brewer_light <- c(
    "lightcyan", "mistyrose", "lightyellow", "plum",
    "lightgreen", "lightblue", "lightgray"
  )
  brewer <- c(
    "cyan", "red", "yellow", "purple",
    "green", "blue", "black"
  )
  maps <- c()
  for (i in seq_along(x)) {
    # map_append = colorRampPalette(RColorBrewer::brewer.pal(3, name = brewer[i]))(n[i])
    map_append <- colorRampPalette(c(brewer_light[i], brewer[i]))(n[i])
    maps <- c(maps, map_append)
  }
  return(maps)
}
