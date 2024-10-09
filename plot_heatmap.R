library(ggplot2)
library(scales) # for muted function

## dataMQ <- read.csv("HeatMap_data.csv")
files <- list.files(".",".csv")
dataMQ <- read.csv(files[1])

########################

MQ.resistance.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                             y = qq,
                                                             fill = cincres)) +
  geom_tile() +
  ylab(label = "API") +
  # xlab(label = "Number of days to access to health care after symptoms of malaria appear") +
  #   facet_grid(~ AWT)
  facet_grid(~ AWT2, switch = "y", scales = "free_y", space = "free_y") +
  scale_fill_gradient(name = "Cumulative resistance cases",
                      low = "#000030",
                      high = "#66CCFF") +
  # low = "#66CCFF",
  #    high = "#000030") +
  theme(strip.placement = "outside",
        strip.text = element_text(size=8,lineheight=5.0),
        axis.text.x = element_text(angle=90, hjust = 0.5,vjust=1,size = 7,face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 8)) +
  labs(title = "Impact of Medicine Quality",
       subtitle = "Access to health care after symptoms of malaria appear (days)",
       x = "Basic reproduction number (R0)")
#  ggtitle(label = "Impact of Medicine Quality")
MQ.resistance.heatmap   

MQ.incidence.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                            y = qq,
                                                            fill = cinc)) +
  geom_tile() +
  ylab(label = "API") +
  # xlab(label = "Number of days to access to health care after symptoms of malaria appear") +
  #   facet_grid(~ AWT)
  facet_grid(~ AWT2, switch = "y", scales = "free_y", space = "free_y") +
  scale_fill_gradient(name = "Incidence cases",
                      low = "#000030",
                      high = "#66CCFF") +
  # low = "#66CCFF",
  #    high = "#000030") +
  theme(strip.placement = "outside",
        strip.text = element_text(size=8,lineheight=5.0),
        axis.text.x = element_text(angle=90, hjust = 0.5,vjust=1,size = 7,face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 8)) +
  labs(title = "Impact of Medicine Quality",
       subtitle = "Access to health care after symptoms of malaria appear (days)",
       x = "Basic reproduction number (R0)")
#  ggtitle(label = "Impact of Medicine Quality")

MQ.incidence.heatmap   

########################





MQ.clincidence.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                              y = qq,
                                                              fill = clinc)) +
  geom_tile() +
  ylab(label = "API") +
  # xlab(label = "Number of days to access to health care after symptoms of malaria appear") +
  #   facet_grid(~ AWT)
  facet_grid(~ AWT2, switch = "y", scales = "free_y", space = "free_y") +
  scale_fill_gradient(name = "Inifcted Clinical (per 1000)",
                      low = "#000030",
                      high = "#66CCFF") +
  # low = "#66CCFF",
  #    high = "#000030") +
  theme(strip.placement = "outside",
        strip.text = element_text(size=8,lineheight=5.0),
        axis.text.x = element_text(angle=90, hjust = 0.5,vjust=1,size = 7,face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 8)) +
  labs(title = "Impact of Medicine Quality",
       subtitle = "Access to health care after symptoms of malaria appear (days)",
       x = "Basic reproduction number (R0)")
#  ggtitle(label = "Impact of Medicine Quality")

MQ.clincidence.heatmap   

########################








MQ.incidence.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                            y = qq,
                                                            fill = cinc)) +
  geom_tile() +
  ylab(label = "API") +
  xlab(label = "Number of days to access to health care") +
  facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(colours = terrain.colors(10))
  
  scale_fill_gradient(name = "Reported incidence",
                      low = "palegreen",
                      high = "red") +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Impact of Medicine Quality")

MQ.incidence.heatmap   


MQ.resistance.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                             y = qq,
                                                             fill = cincres)) +
  geom_tile() +
  ylab(label = "API") +
  xlab(label = "Number of days to access to health care") +
  facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(colours = terrain.colors(10))
  
  scale_fill_gradient(name = "Resistance cases",
                      low = "palegreen",
                      high = "red") +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Impact of Medicine Quality")

MQ.resistance.heatmap   



MQ.resistance.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                             y = qq,
                                                             fill = cincres)) +
  geom_tile() +
  ylab(label = "API") +
  xlab(label = "Number of days to access to health care") +
  facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(colours = terrain.colors(10))
  
  scale_fill_gradient(name = "Resistance cases",
                      low = "palegreen",
                      high = "red") +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Impact of Medicine Quality")

MQ.resistance.heatmap   






g# dataMQ$log.cincres <- log(dataMQ$cincres)

MQ.resistance.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                             y = qq,
                                                             fill = cincres)) +
  geom_tile() +
  xlab(label = "R0") +
  facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(colours = terrain.colors(10))
  
  scale_fill_gradient(name = "Cumulative Incidence",
                      low = "#CCFF77",
                      high = "#FF3300")

MQ.resistance.heatmap   




MQ.incidence.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                            y = qq,
                                                            fill = cinc)) +
  geom_tile() +
  ylab(label = "API") +
  xlab(label = "Number of days to access to health care") +
  facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
  # scale_fill_gradientn(colours = terrain.colors(10))
  
  scale_fill_gradient(name = "Reported incidence",
                      low = "#66CCFF",
                      high = "#000030") +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Impact of Medicine Quality")

MQ.incidence.heatmap   


MQ.resistance.heatmap <- ggplot(data = dataMQ, mapping = aes(x = R0,
                                                             y = qq,
                                                             fill = cincres)) +
  geom_tile() +
  ylab(label = "API") +
  xlab(label = "Number of days to access to health care after infection") +
  facet_grid(~ AWT)
# facet_grid(~ AWT, switch = "x", scales = "free_x", space = "free_x") +
# scale_fill_gradientn(colours = terrain.colors(10))

scale_fill_gradient(name = "Resistance cases",
                    low = "#66CCFF",
                    high = "#000030") +
  theme(strip.placement = "outside",
        axis.text.x = element_text(angle=90, hjust = 0.5,vjust=1,size = 7,face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle(label = "Impact of Medicine Quality")

MQ.resistance.heatmap   
