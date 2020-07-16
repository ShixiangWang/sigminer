logo_data2 <- function(seqs, idor = NULL, method = "bits", stack_width = 0.95, rev_stack_order = F,
                      font, seq_group = 1, seq_type = "auto", namespace = NULL) {

  font_df <- ggseqlogo:::get_font(font)
  if (method == "bits") {
    hh <- ggseqlogo:::bits_method(seqs,
      decreasing = rev_stack_order,
      seq_type = seq_type, namespace = namespace
    )
  }
  else if (method == "probability") {
    hh <- ggseqlogo:::probability_method(seqs,
      decreasing = rev_stack_order,
      seq_type = seq_type, namespace = namespace
    )
  }
  else if (method == "custom") {
    if (seq_type == "auto") {
      seq_type <- ggseqlogo:::guessSeqType(rownames(seqs))
    }
    hh <- ggseqlogo:::matrix_to_heights(seqs, seq_type, decreasing = rev_stack_order)
  }
  else {
    stop("Invalid method!")
  }

  ## Modified by Shixiang
  ## TO support custom indicator
  if (!is.null(idor)) {
    hh$letter2 <- hh$letter
    hh$letter = as.character(idor[hh$letter])
  }
  ## END

  ff <- merge(font_df, hh, by = "letter")
  x_pad <- stack_width / 2
  ff$x <- ggseqlogo:::newRange(ff$x, ff$position - x_pad, ff$position +
    x_pad)
  ff$y <- ggseqlogo:::newRange(ff$y, ff$y0, ff$y1)

  if (!is.null(idor)) {
    cls <- c("x", "y", "letter", "position", "order", "letter2")
  } else {
    cls <- c("x", "y", "letter", "position", "order")
  }

  ff <- as.data.frame(ff)[, cls]
  ff$seq_group <- seq_group
  attr(ff, "seq_type") <- attr(hh, "seq_type")
  ff
}


ggseqlogo2 <- function (data, facet = "wrap", scales = "free_x", ncol = NULL,
          nrow = NULL, idor = NULL, ...)
{
  p = ggplot() + geom_logo2(data = data, idor = idor, ...) + ggseqlogo::theme_logo()
  if (!"list" %in% class(data))
    return(p)
  facet_opts = c("grid", "wrap")
  pind = pmatch(facet, facet_opts)
  facet = facet_opts[pind]
  if (is.na(facet))
    stop("facet option must be set to 'wrap' or 'grid'")
  if (facet == "grid") {
    p = p + facet_grid(~seq_group, scales = scales)
  }
  else if (facet == "wrap") {
    p = p + facet_wrap(~seq_group, scales = scales, nrow = nrow,
                       ncol = ncol)
  }
  return(p)
}

geom_logo2 <- function (data = NULL, method = "bits", seq_type = "auto", namespace = NULL,
          font = "roboto_medium", stack_width = 0.95, rev_stack_order = F,
          col_scheme = "auto", low_col = "black", high_col = "yellow",
          na_col = "grey20", plot = T, idor = NULL, ...)
{
  if (!"ggseqlogo" %in% .packages()) {
    attachNamespace("ggseqlogo")
  }
  if (stack_width > 1 | stack_width <= 0)
    stop("\"stack_width\" must be between 0 and 1")
  if (is.null(data))
    stop("Missing \"data\" parameter!")
  if (!is.null(namespace))
    seq_type = "other"
  all_methods = c("bits", "probability", "custom")
  pind = pmatch(method, all_methods)
  method = all_methods[pind]
  if (is.na(method))
    stop("method must be one of 'bits' or 'probability', or 'custom'")
  if (is.character(data) | is.matrix(data))
    data = list(`1` = data)
  if (is.list(data)) {
    if (is.null(names(data)))
      names(data) = seq_along(data)
    lvls = names(data)
    data_sp = lapply(names(data), function(n) {
      curr_seqs = data[[n]]
      logo_data2(seqs = curr_seqs, idor = idor, method = method, stack_width = stack_width,
                rev_stack_order = rev_stack_order, seq_group = n,
                seq_type = seq_type, font = font, namespace = namespace)
    })
    data = do.call(rbind, data_sp)
    data$seq_group = factor(data$seq_group, levels = lvls)
  }
  if (!plot)
    return(data)
  seq_type = attr(data, "seq_type")
  cs = ggseqlogo:::get_col_scheme(col_scheme, seq_type)
  legend_title = attr(cs, "cs_label")

  ## Modified by Shixiang
  if (!is.null(idor)) {
    data = merge(data, cs, by.x = "letter2", by.y = "letter", all.x = T)
  } else {
    data = merge(data, cs, by = "letter", all.x = T)
  }
  ## END
  data = data[order(data$order), ]
  colscale_gradient = is.numeric(cs$group)
  colscale_opts = NULL
  if (colscale_gradient) {
    colscale_opts = scale_fill_gradient(name = legend_title,
                                        low = low_col, high = high_col, na.value = na_col)
  }
  else {
    tmp = cs[!duplicated(cs$group) & !is.na(cs$group), ]
    col_map = unlist(split(tmp$col, tmp$group))
    colscale_opts = scale_fill_manual(values = col_map, name = legend_title,
                                      na.value = na_col)
  }
  guides_opts = NULL
  if (identical(cs$letter, cs$group))
    guides_opts = guides(fill = F)
  y_lim = NULL
  extra_opts = NULL
  if (method == "tsl") {
    y_lab = "Depleted    Enriched"
    tmp = max(abs(data$y))
    row_a = row_b = data[1, ]
    row_a$y = -tmp
    row_b$y = tmp
    data = rbind(data, row_a, row_b)
    data$facet = factor(data$y > 0, c(T, F), c("Enriched",
                                               "Depleted"))
    extra_opts = NULL
  }
  else if (method == "custom") {
    y_lab = ""
  }
  else {
    y_lab = method
    substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
  }
  data$group_by = with(data, interaction(seq_group, letter,
                                         position))
  data$x = data$x
  logo_layer = layer(stat = "identity", data = data, mapping = aes_string(x = "x",
                                                                          y = "y", fill = "group", group = "group_by"), geom = "polygon",
                     position = "identity", show.legend = NA, inherit.aes = F,
                     params = list(na.rm = T, ...))
  breaks_fun = function(lim) {
    1:floor(lim[2]/1.05)
  }
  list(logo_layer, scale_x_continuous(breaks = breaks_fun,
                                      labels = identity), ylab(y_lab), xlab(""), colscale_opts,
       guides_opts, coord_cartesian(ylim = y_lim), extra_opts)
}
