plot_pca <- function(types = NULL, components = c(1, 2)) {
  df <- self$pca(types = types)
  g <- ggplot(df,
         aes_string(
          x = paste0("PC", components[1]),
          y = paste0("PC", components[2])
        )
    ) +
    geom_point()
  g <- g + self$my_theme
  return(g)
}
