# Mixed model estimation of diagnostic dates
require(plyr)
require(mistr)

diagnostics <-
  ddply(
    droplevels(subset(dat, !is.na(BEGDATE + ENDDATE))),
    ~ CODE + IDENTIFY + BEGDATE + ENDDATE,
    summarise,
    fq = sum(N) ,
    n.obs = length(N)
  )

diagdates <-
  mapply(c, diagnostics$BEGDATE, diagnostics$ENDDATE, SIMPLIFY = F)

diag_mix <-
  mixdist(rep("unif", length(diagdates)),
          diagdates,
          weights = diagnostics$n.obs / sum(diagnostics$n.obs)
          )

mindate <- min(diagnostics$BEGDATE)
maxdate <- max(diagnostics$ENDDATE)

date_seq <- seq(mindate, maxdate, length.out = maxdate - mindate)

plot(
  density(r(diag_mix, sum(diagnostics$n.obs))),
  xlim = c(mindate, maxdate),
  type = 'l',
  lwd = 1,
  col = 'red',
  xlab = "Year",
  ylab = "Probability Density",
  main = "Diagnostics Mixed Model PDF"
)

# Compare against mean diagnostic dates

# Un-weighted density
lines(
  density((
    diagnostics$BEGDATE + diagnostics$ENDDATE
  ) / 2),
  type = 'l',
  lty = 2,
  lwd = 1,
  col = "green"
)

# Frequency weighted density
lines(
  density((diagnostics$BEGDATE + diagnostics$ENDDATE) / 2,
          weights = diagnostics$n.obs / sum(diagnostics$n.obs)
  ),
  type = 'l',
  lty = 2,
  lwd = 1,
  col = "blue"
)
