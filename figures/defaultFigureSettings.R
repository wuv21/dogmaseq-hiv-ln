library(ggplot2)
library(patchwork)
library(png)
library(ggrastr)
library(glue)
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))

BASEPTAXISFONTSIZE <- 8
BASEPTFONTSIZE <- 8
BASEFONTSIZE <- BASEPTFONTSIZE / ggplot2:::.pt
BASELINEWIDTHSIZE <- 12/14


plotTagTheme <- theme(
  plot.tag = element_text(size = 10, face = "plain", family = "sans"))

subplotTheme <- theme(
  plot.margin = unit(c(0,0,0,0), "pt"),
)

pValSymnum <- function(x, showNs = TRUE) {
  tmp <- sapply(x, function(y) {
    if (is.na(y)) {
      return(NA)
    }
    
    if (y < 0.001) {
      return("***")
    } else if (y >= 0.001 & y < 0.01) {
      return("**")
    } else if (y >= 0.01 & y < 0.05) {
      return("*")
    } else {
      return(ifelse(showNs, "ns", ""))
    }
  })
  
  return(tmp)
}

umapPlotThemeNoLeg <- list(
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))
)

unclipStripText <- function(g) {
  g2 <- patchworkGrob(g)
  for (i in which(grepl("strip-r", g2$layout$name))){
    
    g2$grobs[[i]]$grobs[[1]]$layout$clip <- "off"
  }
  
  return(g2)
}

# draw positive error bar only
mean_se_onlyPos <- function(x, mult = 1) {
  x <- stats::na.omit(x)
  se <- mult * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  return(mean + se)
  # data_frame(ymin = mean, ymax = mean)
}

mean_se_onlyPosVertical <- function(x, mult = 1) {
  x <- stats::na.omit(x)
  se <- mult * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  data_frame(ymin = mean, ymax = mean + se)
}


saveFinalFigure <- function(
  plot,
  prefixDir = "../outs",
  fn,
  devices = c("png", "pdf", "pdf_base", "Cairo"),
  addTimestamp = FALSE,
  gheight,
  gwidth
) {
  
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  
  for (d in devices) {
    if (!dir.exists(glue("{prefixDir}/{d}"))) {
      dir.create(glue("{prefixDir}/{d}"))
    }
    if (addTimestamp) {
      currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      if (d == "pdf_base") { 
        gfn <- glue("{prefixDir}/{d}/{fn}_{currentTime}.pdf")
      } else {
        gfn <- glue("{prefixDir}/{d}/{fn}_{currentTime}.{d}")
      } 
    } else {
      if (d == "pdf_base") { 
        gfn <- glue("{prefixDir}/{d}/{fn}.pdf")
      } else if (d == "Cairo") {
        gfn <- glue("{prefixDir}/{d}/{fn}.pdf")
      } else {
        gfn <- glue("{prefixDir}/{d}/{fn}.{d}")
      }
    }
    if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (d == "pdf") {
      ggsave(gfn, plot = plot, dpi = "retina", device = cairo_pdf, width = gwidth, height = gheight, units = "in")
    } else if (d == "pdf_base") {
      ggsave(gfn, plot = plot, dpi = "retina", device = "pdf", width = gwidth, height = gheight, units = "in")
    } else if (d == "Cairo") {
      CairoPDF(width=gwidth, height=gheight, file=gfn, family="sans")
      print(plot)
      dev.off()
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, units = "in")  
    }
  }
}