library(ggplot2)

write.img <- function(plt, path, width=6.5, height=4,
                     fn=png, ext=".png") {
  # Width x Height: 6.5x4" at 300 dpi
  fn(paste0(path, ext), width=width*300, height=height*300)
  tryCatch(
    print(plt),
    finally = dev.off()
  )
}
theme_std <- function(size=32,...) {
  theme_bw(size) +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(
        fill="white",
        color="white"
      ),
      ...
    )
}
theme_vert <- function(size=32,...) {
  theme_std(
    size=size,
    axis.text.x = element_text(
      angle=90, vjust=0.5, hjust=1
    ),
    ...
  )
}
twrap <- scale_x_discrete(labels = function(x) str_wrap(str_to_title(x), width=15))
