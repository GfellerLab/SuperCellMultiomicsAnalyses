
# ------------------------------------------------------------
# Rank-biserial correlation from two vectors
# ------------------------------------------------------------

rbc_from_vectors <- function(x, y) {
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  n1 <- length(x)
  n2 <- length(y)
  
  if (n1 == 0 || n2 == 0) {
    return(NA_real_)
  }
  
  # Pairwise differences:
  # positive -> x > y
  # negative -> x < y
  # zero     -> x = y
  differences <- outer(x, y, "-")
  
  n_greater <- sum(differences > 0)
  n_smaller <- sum(differences < 0)
  
  # Rank-biserial correlation
  RBC <- (n_greater - n_smaller) / (n1 * n2)
  
  RBC
}


HeatmapWilcox <- function(results_plot){
  ggplot(
    results_plot,
    aes(
      x = group2,
      y = group1,
      fill = RBC,
      size = logP
    )
  ) +
    
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.3
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Rank-biserial\ncorrelation"
    ) +
    
    scale_size_continuous(
      range = c(1.5, 6),
      name = expression(-log[10]("p-value"))
    ) +
    
    facet_wrap(
      ~ metric
    ) +
    
    coord_fixed() +
    
    theme_bw() +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      
      axis.title = element_blank(),
      
      strip.background = element_rect(
        fill = "grey95"
      ),
      
      legend.position = "right"
    )
}



ComputePairwiseWilcox <- function(
    bench.results,
    cols = c(
      "method",
      "purity",
      "multimodal_compactness",
      "multimodal_separation"
    )
) {
  
  stopifnot("method" %in% cols)
  
  ## Check columns
  missing.cols <- setdiff(cols, colnames(bench.results))
  if (length(missing.cols) > 0) {
    stop(
      "Missing columns: ",
      paste(missing.cols, collapse = ", ")
    )
  }
  
  ## Keep requested columns
  dat <- bench.results[, cols, drop = FALSE]
  dat$metacell <- rownames(dat)
  
  ## Long format
  dat.long <- dat %>%
    tidyr::pivot_longer(
      cols = -c(metacell, method),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      method = as.character(method),
      metric = as.character(metric)
    )
  
  ## Pairwise Wilcoxon tests
  results <- dat.long %>%
    dplyr::group_by(metric) %>%
    dplyr::group_modify(~ {
      
      df <- .x
      
      methods <- unique(df$method)
      
      pairs <- combn(
        methods,
        2,
        simplify = FALSE
      )
      
      pair.results <- purrr::map_dfr(pairs, function(pair) {
        
        group1 <- pair[1]
        group2 <- pair[2]
        
        x <- df$value[df$method == group1]
        y <- df$value[df$method == group2]
        
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
        
        wt <- wilcox.test(
          x,
          y,
          exact = FALSE
        )
        
        tibble::tibble(
          group1 = group1,
          group2 = group2,
          n1 = length(x),
          n2 = length(y),
          W = unname(wt$statistic),
          RBC = rbc_from_vectors(x, y),
          p = wt$p.value
        )
        
      })
      
      pair.results %>%
        dplyr::mutate(
          padj = p.adjust(p, method = "BH"),
          logP = -log10(pmax(p, .Machine$double.xmin)),
          logFDR = -log10(
            pmax(padj, .Machine$double.xmin)
          )
        )
      
    }) %>%
    dplyr::ungroup()
  
  return(results)
}
