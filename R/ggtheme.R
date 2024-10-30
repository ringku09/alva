library(grid)
library(gridExtra)
library(ggthemes)
#library(extrafont)
library(latex2exp)
library(scales)
# font_import()
# loadfonts(device = "win", quiet = TRUE)
# fonts()

common_legend<-function(gplot){
  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


theme_Publication <- function(base_size=14,
                              base_family="Helvetica",
                              title_size = rel(1),
                              axix_size = rel(1),
                              lab_size = rel(1),
                              leg_size = rel(0.76),
                              y_lab = element_text(angle=90,vjust =2),
                              x_lab = element_text(vjust = -0.2),
                              tag_size = rel(0.7),
                              top_mar = 5,
                              bot_mar = 5,
                              ...) {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(title_size), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            plot.tag = element_text(face = "bold", size = tag_size),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = axix_size),
            axis.title.y = y_lab,
            axis.title.x = x_lab,
            axis.text = element_text(size = lab_size),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size= leg_size),
            legend.key.size= unit(0.4, "cm"),
            legend.margin = margin(0,unit = "cm"),
            legend.title = element_text(face="italic",size = leg_size),
            plot.margin=unit(c(top_mar,3,bot_mar,3),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold", size= rel(0.6))
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill",
                 "Publication",
                 manual_pal(values = c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
                                       "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")),
                 ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour",
                 "Publication",
                 manual_pal(values = c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
                                                "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")),
                 ...)

}



