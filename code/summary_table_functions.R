# Functions to summarize artifacts 

require(tidyverse)
require(janitor)

summary_count_depth <- function(x, grp_sort = FALSE, total = FALSE) {
  count_depth_sum <- x %>%
    group_by(.add = TRUE) %>%
    select(UNIT_ID,
           PROVENIENCE_ID,
           N,
           START_DEPTH,
           END_DEPTH) %>%
    summarise(
      n.units = length(unique(UNIT_ID)),
      n.provs = length(unique(PROVENIENCE_ID)),
      n.obj = sum(N),
      start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
      end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1)
    ) %>%
    arrange(desc(n.units),
            desc(n.provs),
            desc(n.obj),
            .by_group = grp_sort)
  if (total) {
    count_depth_sum <- count_depth_sum %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        n.obj
      )
  }
  return(count_depth_sum)
}

summary_date_depth <- function(x, grp_sort = FALSE, total = FALSE) {
  date_depth_sum <- x %>%
    group_by(.add = TRUE) %>%
    select(UNIT_ID,
           PROVENIENCE_ID,
           N,
           START_DEPTH,
           END_DEPTH,
           BEGDATE,
           ENDDATE) %>%
    summarise(
      n.units = length(unique(UNIT_ID)),
      n.provs = length(unique(PROVENIENCE_ID)),
      n.obj = sum(N),
      start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
      end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1),
      beg.date = round(weighted.mean(BEGDATE, N, na.rm = TRUE), 0),
      end.date = round(weighted.mean(ENDDATE, N, na.rm = TRUE), 0)
    ) %>%
    arrange(beg.date,
            end.date,
            desc(n.units),
            desc(n.provs),
            desc(n.obj),
            .by_group = grp_sort)
  if (total) {
    date_depth_sum <- date_depth_sum %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        n.obj
      )
  }
  return(date_depth_sum)
}

# Cross-tabulations -----
