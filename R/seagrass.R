library(here)
library(tidyverse)
library(patchwork)

# Seagrass losses in upper Harbor (1988-2022)
dat <- read.csv(here("data-raw/CH.seagrass.csv"))
years <- na.omit(dat)$Year
p1 <- ggplot( dat, aes(y=CH.Seagrass.ha / 1000, x=Year) ) + 
       geom_bar(stat="identity", fill = '#035172', alpha = 0.7) +
       labs(
              x = NULL,
              y = "Coverage (1000 x hectares)", 
              subtitle = "Charlotte Harbor Seagrass Coverage"
       ) +
       theme_minimal(base_size = 16) +
       scale_x_continuous(breaks = years, labels = years ) +
       scale_y_continuous(expand = c(0, 0)) + 
       theme( 
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.x = element_blank()
       )

# Recent seagrass losses in upper Harbor
dat2 <- read.csv(here("data-raw/CH.seagrass_2022.csv"))
dat2$Year <- factor( dat2$Year, levels = c("2018","2020","2022") )
dat2$Stratum <- factor( dat2$Stratum, levels = c("TPR","TMR","CHW","CHE","CHL","CHZ","LB") )
p2 <- ggplot(dat2, aes(fill=Year, y=Hectares / 1000, x=Stratum)) + 
       geom_bar(position="dodge", stat="identity", alpha = 0.7) +
       geom_text( aes( y=Hectares / 1000, x=Stratum,
                       label=sub("NA%","",paste0(Loss_perc_2018,"%")) ),
                   position = position_dodge(width=1),
                   size  = 3, vjust = -0.2) +
       scale_fill_manual(values = c(rgb(0.1,0.5,0.2,1),
                                     rgb(0.2,0.6,0.1,1),rgb(0.4,0.75,0.1,1))) +
       
       labs(
              x = NULL, 
              y = 'Coverage (1000 x hectares)', 
              subtitle = 'Recent change in upper harbor'
       ) +
       theme_minimal(base_size = 16) +
       scale_y_continuous(expand = c(0, 0)) + 
       theme(
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.x = element_blank(), 
              legend.position = 'bottom'
       )

p <- p1 + p2 + plot_layout(ncol = 1, axis_titles = 'collect_y')

png(here('figs/seagrass.png'), width = 8, height = 6, units = 'in', res = 300)
print(p)
dev.off()

