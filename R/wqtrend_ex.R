library(wqtrends)
library(tidyverse)
library(patchwork)
library(here)

data(rawdat)

tomod <- rawdat %>% 
  filter(station %in% 30) %>% 
  filter(param %in% 'chl')

# model 
mod <- anlz_gam(tomod, trans = 'log10')

# get metric predictions
doystr <- 213
doyend <- 304
mets <- anlz_metseason(mod, doystr = doystr, doyend = doyend)

# add thrs info
thrsv <- 0.75
toplo1 <- mets |> 
  mutate(
    thrs = case_when(
      met - se > thrsv ~ 'Above',
      met + se < thrsv ~ 'Below',
      TRUE ~ 'Near'
    ), 
    thrs = factor(thrs, levels = c('Below', 'Near', 'Above'))
  )

maxslo <- 0.1
minslo <- -1 * maxslo
ymin <- 0 # 0.25 * min(toplo1$met - toplo1$se)
ylab <- bquote(log[10] ~ 'Chl-a (μg/L)')
yrstr1 <- 1990
yrend1 <- 1999
idx <- 0

slos <- NULL

while(yrend1 <= max(tomod$yr)){

  cat(yrend1, '\t')

  trns <- anlz_mixmeta(mets, yrstr = yrstr1, yrend = yrend1)

  addv <- tibble(
      yrend = yrend1, 
      slope = coefficients(trns)[[2]], 
      pval = coefficients(summary(trns))[[8]] 
    ) |> 
    mutate(
      dirv = case_when(
        pval < 0.05 & slope > 0 ~ 'Inc',
        pval < 0.05 & slope < 0 ~ 'Dec',
        TRUE ~ 'None'
      ),
      dirv = factor(dirv, levels = c('Inc', 'None', 'Dec'))
    )

  slos <- rbind(slos, addv)

  fits <- data.frame(
    yr = seq(yrstr1, yrend1, length = 50)
    ) %>% 
    dplyr::mutate( 
      met = predict(trns, newdata = data.frame(yr = yr)), 
      se = predict(trns, newdata = data.frame(yr = yr), se = T)[, 2]
    )

  toplo1sub <- toplo1 |> 
    dplyr::filter(yr <= yrend1)

  p1 <- ggplot(toplo1, aes(x = yr, y = met)) + 
    geom_hline(yintercept = thrsv, linetype = 'dashed', color = 'gray', linewidth = 1) +
    geom_point(col = "deepskyblue3") +
    geom_errorbar(aes(ymin = met - se, ymax = met + se), col = "deepskyblue3") +
    geom_ribbon(data = fits, ggplot2::aes(ymin = met - se, ymax = met + se), fill = 'pink', alpha = 0.4) +
    geom_line(data = fits, aes(y = met), color = 'pink') +
    geom_point(data = toplo1sub, y = ymin, aes(fill = thrs), 
      color = 'black', pch = 22, size = 3, show.legend = T) +
    scale_fill_manual(
      values = c('Below' = 'white', 'Near' = 'grey', 'Above' = 'black'), 
      drop = F
    ) +
    coord_cartesian(
      ylim = c(ymin, NA)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.025, 0.025))) +
    theme_minimal() + 
    theme(
      axis.text.x = element_blank(), 
      legend.position = 'bottom'
    ) +
    labs(
      x = NULL, 
      y = ylab,
      fill = 'Threshold'
    )

  p2 <- ggplot(slos, aes(x = yrend, y = slope)) + 
    geom_hline(yintercept = 0, linetype = 'solid', color = 'gray') +
    geom_point(size = 3, aes(fill = dirv), color = 'black', pch = 21, show.legend = T) +
    theme_minimal() + 
    theme(legend.position = 'bottom') +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
    scale_fill_manual(
      values = c('Inc' = 'tomato1', 'None' = 'white', 'Dec' = 'deepskyblue3'), 
      drop = F
    ) +
    labs(
      x = NULL, 
      y = bquote(log[10] ~ 'slope'),  
      fill = 'Trend'
    ) +
    coord_cartesian(
      xlim = range(mets$yr),
      ylim = c(minslo, maxslo)
    )

  p <- p1 + p2 + plot_layout(ncol = 1, heights = c(1, 1), guides = 'collect') & 
    theme(
      text = element_text(size = 13), 
      panel.grid.minor.y = element_blank(), 
      legend.margin = margin(0, 0, 0, 0), 
      legend.position = 'bottom'
    )

  idx <- idx + 1
  outfl <- here(paste0('figs/wqtrend_ex/p', sprintf('%02d', idx), '.png'))
  png(outfl, width = 9, height = 4.25, units = 'in', res = 300)
  print(p)
  dev.off()

  yrstr1 <- yrstr1 + 1
  yrend1 <- yrend1 + 1

}

system('magick convert -delay 50 -loop 0 figs/wqtrend_ex/*.png figs/wqtrend_ex/wqtrend.gif')
