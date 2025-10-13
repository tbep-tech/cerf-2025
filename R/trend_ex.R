library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(mgcv)      # for GAM
library(segmented) # for segmented regression
library(trend)     # for Seasonal Mann-Kendall test

# Set seed for reproducibility
set.seed(123)

# Generate synthetic water quality data (e.g., chlorophyll-a)
# 5 years of monthly observations with seasonal pattern and subtle trend
dates <- seq(as.Date("2019-01-01"), as.Date("2023-12-01"), by = "month")
n <- length(dates)
time_numeric <- as.numeric(dates - min(dates)) / 365.25

# Create data with seasonal cycle, slight non-linear trend, and noise
seasonal <- 3 * sin(2 * pi * time_numeric)
trend <- 0.3 * time_numeric + 0.15 * time_numeric^2
noise <- rnorm(n, 0, 1.5)
chlorophyll <- 8 + seasonal + trend + noise
chlorophyll <- pmax(chlorophyll, 0.1)  # ensure positive values

df <- data.frame(
  date = dates,
  time_numeric = time_numeric,
  chlorophyll = chlorophyll,
  month = as.numeric(format(dates, "%m"))
)

# 1. Raw Data Plot
p1 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  theme_minimal(base_size = 14) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 2. Linear Regression (simple trend line)
lm_model <- lm(chlorophyll ~ time_numeric, data = df)
df$lm_pred <- predict(lm_model)
lm_pval <- summary(lm_model)$coefficients[2, 4]
lm_slope <- coef(lm_model)[2]

p2 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  geom_line(aes(y = lm_pred), color = "red", linewidth = 1.2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 3. GAM (Generalized Additive Model) - captures non-linear trends
poly_model <- lm(chlorophyll ~ poly(time_numeric, 2, raw = TRUE), data = df)
df$poly_pred <- predict(poly_model)

p3 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  geom_line(aes(y = poly_pred), color = "purple", linewidth = 1.2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 4. Polynomial Regression (overfitting example)
# Fit a high-degree polynomial (degree 8) - too flexible!
poly_model <- lm(chlorophyll ~ poly(time_numeric, 8, raw = TRUE), data = df)
df$poly_pred <- predict(poly_model)

p4 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  geom_line(aes(y = poly_pred), color = "darkgreen", linewidth = 1.2) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 5. Seasonal Mann-Kendall Test (non-parametric, accounts for seasonality)
mk_results <- smk.test(ts(df$chlorophyll, frequency = 12))

p5 <- ggplot(df, aes(x = date, y = chlorophyll, color = factor(month))) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "darkgreen", linewidth = 1.2, linetype = "dashed") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)",
       color = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none") +
  coord_cartesian(ylim = c(2, 19))

# 6. Deseasonalized Trend (removes seasonal component first)
# Fit seasonal component
seasonal_model <- lm(chlorophyll ~ sin(2*pi*time_numeric) + cos(2*pi*time_numeric), 
                     data = df)
df$deseasonalized <- df$chlorophyll - 
  (predict(seasonal_model) - mean(df$chlorophyll))

# Trend on deseasonalized data
deseas_trend <- lm(deseasonalized ~ time_numeric, data = df)
df$deseas_pred <- predict(deseas_trend)

p6 <- ggplot(df, aes(x = date)) +
  geom_point(aes(y = deseasonalized), size = 3, alpha = 0.4, color = "steelblue") +
  geom_line(aes(y = deseas_pred), color = "darkred", linewidth = 1.2) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 7. Segmented Regression (detects change points)
seg_model <- tryCatch({
  lm_base <- lm(chlorophyll ~ time_numeric, data = df)
  segmented(lm_base, seg.Z = ~time_numeric, psi = 2.5)
}, error = function(e) NULL)

df$seg_pred <- predict(seg_model)
breakpoint <- seg_model$psi[2]
breakpoint_date <- min(df$date) + breakpoint * 365.25
  
p7 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  geom_line(aes(y = seg_pred), color = "orange", linewidth = 1.2) +
  geom_vline(xintercept = as.numeric(breakpoint_date), 
              linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))

# 8. Uncertainty visualization - what if measurements have error?
# Add realistic measurement uncertainty (±10-20% typical for chlorophyll)
df$uncertainty <- df$chlorophyll * runif(n, 0.10, 0.20)
df$lower <- df$chlorophyll - df$uncertainty
df$upper <- df$chlorophyll + df$uncertainty

p8 <- ggplot(df, aes(x = date, y = chlorophyll)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0, color = "steelblue", alpha = 0.5) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", expand = c(0.05, 0.05)) +
  geom_point(size = 3, alpha = 0.4, color = "steelblue") +
  labs(x = NULL, y = "Chlorophyll-a (μg/L)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(2, 19))


ht <- 3.5
wd <- 8
dp <- 150
ggsave(here("figs/trend_ex/p1.png"), p1, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p2.png"), p2, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p3.png"), p3, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p4.png"), p4, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p5.png"), p5, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p6.png"), p6, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p7.png"), p7, width = wd, height = ht, dpi = dp)
ggsave(here("figs/trend_ex/p8.png"), p8, width = wd, height = ht, dpi = dp)
