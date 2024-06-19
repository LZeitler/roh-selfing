fread <- data.table::fread
require(tidyr)
require(dplyr)
require(ggplot2)
require(cowplot)
require(grid)
# require(gridExtra)
require(viridis)
require(colorspace)
# require(reshape2)
require(ggrepel)
require(RColorBrewer)
options(max.print=999)

theme_set(theme_cowplot())

## a function to save plots quickly in pdf and png
saveplot <- function(x, name, width=10, height=6){
    ggsave(paste0(name, ".pdf"),
           x, width = width, height = height, limitsize = F)
    ggsave(paste0(name, ".png"),
           x, width = width, height = height, limitsize = F,
	   type="cairo",
           bg = "white")
}

mytimestamp <- function() as.character(system('date +%Y%m%d%H%M%S',T))

