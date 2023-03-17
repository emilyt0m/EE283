library(tidyverse)
install.packages("nycflights13")
library(nycflights13)
flights

library(ggplot2)

P1 <- ggplot(flights, aes(x=distance,y=arr_delay)) +
  geom_point(size=1) + 
  geom_smooth() 

temp_flights <- flights %>%
  group_by(carrier) %>%
  summarize(m_arr_delay = mean(arr_delay,na.rm=TRUE))

P2 <- ggplot(temp_flights, aes(x=carrier,y= m_arr_delay)) +
  geom_bar(stat="identity") +
  theme(axis.text.y = element_text(size = 4)) +
  theme(axis.text.x = element_text(size = 4)) +
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

P3 <- ggplot(flights, aes(x=carrier, y=arr_delay)) + 
  geom_boxplot(outlier.size=.1) +
  theme(axis.text.y = element_text(size = 4)) +
  theme(axis.text.x = element_text(size = 3)) +
  theme(axis.title = element_text(size = 8))

P4 <- ggplot(flights, aes(x=arr_delay)) +
  geom_histogram() +
  theme(axis.text.y = element_text(size = 4)) +
  theme(axis.title = element_text(size = 8))

library(gridExtra)
library(grid)

grid.arrange(P1,P2,P3,P4,ncol=2)

lay <- rbind(c(1,1,1,2),
             c(1,1,1,3),
             c(1,1,1,4))

tiff("figure1.tiff", width = 7, height = 6, units = "in", res=600)
grid.arrange(P1,P2,P3,P4, layout_matrix = lay)
graphics.off()

