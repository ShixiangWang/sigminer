library(hexSticker)

sticker(
  package = "sigminer",                     # package name to display on sticker
  p_size = 24,                          # size of package name
  p_y = 1.5,                            # y of package name
  #p_color = "#C9B128",                  # color of package name
  subplot = "inst/figures/barchart.png",          # sticker feature
  s_x = 1,                          # x of feature
  s_y = .8,                           # y of feature
  s_width = .6,                        # width of feature - maintains aspect ratio
  h_size = 2,                           # border
  h_color = "#C9B128",                  # color of border
  h_fill = "#f9690e",                     # color of background
  url = "github.com/ShixiangWang/sigminer",   # url at the bottom
  u_color = "white",                    # color of url at the bottom
  u_size = 3.5,                         # size of url at the bottom
  filename = "inst/figures/sigminer.png"                 # location to save the image
)
