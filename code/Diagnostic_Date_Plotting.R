# Required packages
require(tidyverse)
require(knitr)
require(janitor)
require(mistr)
require(ggplot2)

dat.diag <- dat %>% filter(!(CODE %in% code_exclude),
                           CODE %in% code_diagnostic)

diag.art <- dat.diag %>%
  group_by(SITE_ID, BEGDATE, ENDDATE, CODE, IDENTIFY) %>%
  select(N, START_DEPTH, END_DEPTH) %>%
  summarise(
    fq = length(N),
    nisp = sum(N),
    start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 2),
    end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Straight-Florentino Farmstead Site ----
# SITE_ID == 1

mysite <- 7

diagnostics <- diag.art %>% filter(SITE_ID == mysite)

diagdates <-
  mapply(c, diagnostics$BEGDATE, diagnostics$ENDDATE, SIMPLIFY = F)

diag_mix_eq <-
  mixdist(rep("unif", length(diagdates)),
          diagdates,
          weights = rep(1 / length(diagdates), length(diagdates)))

diag_mix <-
  mixdist(rep("unif", length(diagdates)),
          diagdates,
          weights = diagnostics$fq / sum(diagnostics$fq))

mindate <- min(diagnostics$BEGDATE) - 10
maxdate <- max(diagnostics$ENDDATE) + 10
# mindate <- 1800
# maxdate <- 2000

date_seq <- seq(mindate, maxdate, length.out = maxdate - mindate)

plot(
  date_seq,
  # d(diag_mix, date_seq),
  d(diag_mix_eq, date_seq),
  type = 'l',
  lwd = 2,
  col = 'gray',
  xlab = "Year",
  ylab = "Density",
  main = paste(
    "Distribution of Diagnostic Artifact Dates for ",
    filter(X581450_Sites, Site_ID == mysite) %>% pull(Site),
    ".",
    sep = ""
  ),
  xlim = c(mindate, maxdate),
  # ylim = c(0, 0.012)
)

# Random sampling of PDF
lines(
  density(replicate(1000, r(
      diag_mix_eq, sum(diagnostics$fq)
      # diag_mix, sum(diagnostics$fq)
    )),
    kernel = "biweight",bw="nrd"),
  lwd = 2,
  lty = 2,
  col = 'red'
)

# Compare against mean diagnostic dates
lines(
  density((diagnostics$BEGDATE + diagnostics$ENDDATE) / 2,
          kernel = "biweight",
          bw = "nrd"
  ),
  lwd = 2,
  lty = 3,
  col = "blue"
)

legend(
  "topright",
  legend = c("Empirical Density", "Estimated Density", "Mean Date"),
  col = c("gray", "red", "blue"),
  lwd = 2,
  lty = 1:3,
  cex = 0.8
)

autoplot(diag_mix_eq)
