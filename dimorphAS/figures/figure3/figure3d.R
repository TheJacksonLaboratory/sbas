## Create Figure 3d -- the 3d data (NOTE- the Python script is working better at this time and we are not using this R script!)
# Removal all variables from workspace
rm(list=ls())

data <- read.table("3dtable.txt", header = TRUE, sep = "\t")

library(plotly)
library(magrittr)

mypal = pal_npg("npg", alpha = 0.7)(9)


# there is one outlier skipped that is far away from everybody 
plot_ly(data, x = ~Inclusion, y = ~Skip, z = ~Expression, alpha = 0.8, 
        marker = list(color = ~Sex, colorscale = c("blue", "red"), showscale=TRUE,symbol = "circle", sizemode = "diameter", size = 5)) %>%
  layout(scene = list(xaxis = list(title = 'Inclusion', range = c(20, 50)),
                      yaxis = list(title = 'Skip', range = c(0, 120)),
                      zaxis = list(title = 'Expression')))


