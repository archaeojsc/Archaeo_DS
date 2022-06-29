# Geocomputation with R
# https://geocompr.robinlovelace.net/index.html

# Spatial Data Science
# with applications in R
# https://keen-swartz-3146c4.netlify.app/

# Required packages ----
require(tidyverse)
require(sf)
require(terra)
require(knitr)
require(kableExtra)
require(janitor)
require(mistr)

# Data setup ----

## Load spatial data ----

# Shapefile locations

APE.filename <- "maps/APE.shp"
SITE.filename <- "maps/Site boundry.shp"
Testing_Loc.filename <- "maps/604828_testing_pt.shp"
Units.filename <- "maps/Phase II Units.shp"

# Load spatial data as simple features

spatial_APE <- st_read(dsn = APE.filename)
spatial_Site <- st_read(dsn =SITE.filename)
spatial_Test_loc <- st_read(dsn = Testing_Loc.filename)
spatial_Units <- st_read(dsn = Units.filename)


# spatial_APE_s4 <- terra::vect(APE.filename)

## Old version using rgdal

# statial_APE <-
#   readOGR(dsn = path.expand("maps"),
#           layer = "APE",
#           stringsAsFactors = TRUE)

# Site boundary



# spatial_Site <-
#   readOGR(dsn = path.expand("maps"),
#           layer = "Site boundry",
#           stringsAsFactors = TRUE)

# Testing

# spatial_test_loc <-
#   readOGR(dsn = path.expand("maps"),
#           layer = "604828_testing_pt",
#           stringsAsFactors = TRUE)
# 
# spatial_units <- 
#   readOGR(dsn = path.expand("maps"),
#           layer = "Phase II Units",
#           stringsAsFactors = TRUE)



## Combine data frames from each project phase

Provenience <-
  bind_rows(Provenience_P1, Provenience_P2, .id = "PHASE")

Catalog <- bind_rows(Catalog_P1, Catalog_P2, .id = "PHASE")

## Create data frame of Units

## Extract unique units from provenience data, clean labeling of unit types,
## assign unique unit id's, and join with spatial data

# Units <-
#   Provenience %>%
#   select(SITE_ID, PHASE, PROJECT_PHASE, STP, UNIT) %>%
#   distinct() %>%
#   mutate(UNIT_TYPE = case_when(is.na(UNIT) ~ "STP",
#                                is.na(STP) ~ "EU"),
#          .before = STP) %>%
#   mutate(UNIT_LABEL = case_when(is.na(UNIT) ~ STP, is.na(STP) ~ UNIT),
#          +.before = STP) %>%
#   select(-STP, -UNIT) %>%
#   inner_join(
#     spatial_test_loc@data,
#     by = c(
#       'PROJECT_PHASE' = 'PHASE',
#       'UNIT_TYPE' = 'U_TYPE',
#       'UNIT_LABEL' = 'U_LABEL'
#     )
#   ) %>%
#   select(-fid_1,-layer,-path) %>%
#   relocate(PROJECT, .before = PHASE)
# 

### Testing joins with sf spatial data.frames

Units <- Provenience %>%
  select(SITE_ID, PHASE, PROJECT_PHASE, STP, UNIT) %>%
  distinct() %>%
  mutate(UNIT_TYPE = case_when(is.na(UNIT) ~ "STP",
                               is.na(STP) ~ "EU"),
         .before = STP) %>%
  mutate(UNIT_LABEL = case_when(is.na(UNIT) ~ STP,
                                is.na(STP) ~ UNIT),
         .before = STP) %>%
  select(-STP,-UNIT) %>%
  group_by(PHASE, UNIT_TYPE, UNIT_LABEL) %>%
  mutate(UNIT_ID = cur_group_id(),
         .before =UNIT_TYPE)

test <- spatial_Test_loc %>%
  inner_join(Units,
             by = c(
               'PHASE' = 'PROJECT_PHASE' ,
               'U_TYPE'= 'UNIT_TYPE',
               'U_LABEL' = 'UNIT_LABEL'
             )) %>%
  select(-fid_1,-layer,-path, -POINT_X, -POINT_Y) 

# %>%
#   inner_join(
#     spatial_Test_loc,
#     by = c(
#       'PROJECT_PHASE' = 'PHASE',
#       'UNIT_TYPE' = 'U_TYPE',
#       'UNIT_LABEL' = 'U_LABEL'
#     )
#   ) %>%
#   select(-fid_1, -layer, -path) %>% 
#   relocate(PROJECT, .before = PHASE)
# 
  
## Merge provenience, catalog, assemblage, and code list fields to new data
## frame

dat <-
  right_join(
    select(Provenience, all_of(
      c(
        'PHASE',
        'PROJECT_PHASE',
        'SITE_ID',
        'PROVENIENCE_ID',
        'STP',
        'UNIT',
        'LEVEL',
        'STRATUM',
        'START_DEPTH',
        'END_DEPTH'
      )
    )),
    left_join(
      select(Catalog, all_of(
        c('PHASE',
          'PROVENIENCE_ID',
          'CATALOG_ID',
          'N',
          'CODE')
      )),
      right_join(
        code_assemblage,
        select(code_list, all_of(
          c(
            'CODE',
            'MATNAME',
            'TYPENAME',
            'CATNAME',
            'IDENTIFY',
            'BEGDATE',
            'ENDDATE'
          )
        )),
        by = 'CODE',
        keep = FALSE
      ),
      by = 'CODE',
      keep = FALSE
    ),
    by = c('PHASE', 'PROVENIENCE_ID'),
    keep = FALSE
  )

## Data cleaning ----

# Reorganize columns, combine unit labels into single column

dat <-
  dat %>%
  relocate(SITE_ID, PROVENIENCE_ID, CATALOG_ID, .after = PHASE) %>%
  mutate(UNIT_TYPE = case_when(is.na(UNIT) ~ "STP", is.na(STP) ~ "EU"),
         .before = STP) %>%
  mutate(UNIT_LABEL = case_when(is.na(UNIT) ~ STP, is.na(STP) ~ UNIT),
         .before = STP) %>%
  select(-STP, -UNIT) %>%
  
  # Assign unique Unit IDs
  group_by(PHASE, UNIT_TYPE, UNIT_LABEL) %>%
  mutate(UNIT_ID = cur_group_id(),
         .before = PROVENIENCE_ID)

# Diagnostic codes ----

code_diagnostic <- code_list %>%
  filter(!is.na(BEGDATE) | !is.na(ENDDATE)) %>%
  pull(CODE)

# Excludes ----

file.name_excl <- "604828_site_excl.md"

if (file.exists(file.name_excl)) {
  file.remove(file.name_excl)
}

## Excludes for all project phases ----

cat(
  paste0("\n",
         "# Artifacts excluded from analysis",
         "\n"),
  sep = "\n",
  file = file.name_excl,
  append = TRUE
)

cat(
  dat %>%
    filter(CODE %in% code_exclude) %>%
    group_by(IDENTIFY) %>%
    summary_count_depth(total = TRUE) %>%
    kable(
      format = 'markdown',
      row.names = FALSE,
      digits = 1,
      caption ="Artifacts excluded from analysis, all project phases.",
      align = "lrrrrr"
    ), 
  sep = "\n",
  file = file.name_excl,
  append = TRUE
)


## Excludes by Project Phase ----

cat(
  paste0("\n",
         "## Artifacts excluded from analysis by project phase"
         ),
  sep = "\n",
  file = file.name_excl,
  append = TRUE
)

for (p in unique(dat$PHASE)) {
  if (nrow(filter(dat, PHASE == p, CODE %in% code_exclude)) > 0) {
    cat("\n", file = file.name_excl, append = TRUE)
    cat(
      dat %>%
        filter(PHASE == p, CODE %in% code_exclude) %>%
        group_by(IDENTIFY) %>%
        summary_count_depth(total = TRUE) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          align = "lrrrrr",
          caption = paste(
            "Artifacts excluded from analysis for ",
            unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
            ".",
            sep = ""
          )
        ),
      sep = "\n",
      file = file.name_excl,
      append = TRUE
    )
  }
}


# Assemblages ----

file.name_assemblages <- "604828_site_assemblages.md"

if (file.exists(file.name_assemblages)) {
  file.remove(file.name_assemblages)
}

## Assemblages for all project phases ----

cat(
  paste0("\n",
         "# Artifact assemblages",
         "\n"),
  sep = "\n",
  file = file.name_assemblages,
  append = TRUE
)

cat(
  dat %>%
    filter(!(CODE %in% code_exclude)) %>%
    group_by(ASSEMBLAGE) %>%
    summary_count_depth(total = TRUE) %>%
    kable(
      format = 'markdown',
      row.names = FALSE,
      digits = 1,
      caption ="Artifact assemblages, all project phases.",
      align = "lrrrrr"
    ), 
  sep = "\n",
  file = file.name_assemblages,
  append = TRUE
)

## Assemblages by project phase ----



# Artifact summaries by assemblage for each site

for (s in Site$SITE_ID) {
  cat(
    paste(
      "\n",
      "# Finds by Assemblages for ",
      filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
      "\n",
      sep = ""
    ),
    sep = "\n",
    file = file.name_assemblages,
    append = TRUE
  )
  
  cat(
    dat %>%
      filter(SITE_ID == s, !(CODE %in% code_exclude)) %>% 
      group_by(ASSEMBLAGE) %>%
      select(PROVENIENCE_ID, N, START_DEPTH, END_DEPTH) %>%
      summarise(
        n.provs = length(unique(PROVENIENCE_ID)),
        n.obj = sum(N),
        start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
        end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1),
        .groups = 'drop'
      ) %>%
      arrange(desc(frequency), desc(NISP), .by_group = FALSE) %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        NISP
      ) %>%
      kable(
        format = 'markdown',
        row.names = FALSE,
        digits = 1,
        caption = paste(
          "Summary of assemblages for ",
          filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
          ".",
          sep = ""
        ), 
        align = "lrrrr"
      ),
    sep = "\n",
    file = file.name_assemblages,
    append = TRUE
  )
  
  for (a in distinct(filter(dat, SITE_ID == s), ASSEMBLAGE) %>% pull()) {
    cat("\n",
        file = file.name_assemblages,
        append = TRUE)
    cat(
      dat %>%
        filter(SITE_ID == s, ASSEMBLAGE == a , !(CODE %in% code_exclude)) %>%
        group_by(TYPENAME, IDENTIFY) %>%
        select(PROVENIENCE_ID, N, START_DEPTH, END_DEPTH) %>%
        summarise(
          n.provs = length(unique(PROVENIENCE_ID)),
          n.obj = sum(N),
          start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
          end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1),
          .groups = 'drop_last'
        ) %>%
        # filter(frequency == 1) %>%
        mutate(art_type_sum = sum(frequency)) %>%
        arrange(
          desc(art_type_sum),
          desc(frequency),
          desc(NISP),
          .by_group = FALSE
        ) %>%
        select(-art_type_sum) %>%
        adorn_totals(
          .,
          where = "row",
          fill = "--",
          na.rm = TRUE,
          name = "Total",
          NISP
        ) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          caption = paste(
            "Summary of ",
            a,
            " assemblage finds for ",
            filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
            ".",
            sep = ""
          ), 
          align = "llrrrr"
        ),
      file = file.name_assemblages,
      append = TRUE,
      sep = "\n"
    )
  }
}

# Artifact prevalence -----------------------------------------------------

file.name_prevalence <- "604828_site_prevalence.md"

if (file.exists(file.name_prevalence)) {
  file.remove(file.name_prevalence)
}

for (s in Site$SITE_ID) {
  cat(
    paste(
      "\n",
      "# Artifact Prevalence for ",
      filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
      "\n",
      sep = ""
    ),
    sep = "\n",
    file = file.name_prevalence,
    append = TRUE
  )
  cat(
    dat %>%   filter(SITE_ID == s, !(CODE %in% code_exclude)) %>%
      group_by(TYPENAME, IDENTIFY) %>%
      select(PROVENIENCE_ID, N, START_DEPTH, END_DEPTH) %>%
      summarise(
        n.provs = length(unique(PROVENIENCE_ID)),
        n.obj = sum(N),
        start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
        end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1)
      ) %>%
      filter(frequency > 1) %>%
      mutate(art_assem_sum = sum(frequency)) %>%
      arrange(desc(frequency), desc(NISP), .by_group = FALSE) %>%
      arrange(desc(art_assem_sum), .by_group = FALSE) %>%
      select(-art_assem_sum) %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        NISP
      ) %>%
      kable(
        format = 'markdown',
        row.names = FALSE,
        digits = 1,
        caption = paste(
          "Most prevalent artifacts for ",
          filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
          ".",
          sep = ""
        )
      ),
    sep = "\n",
    file = file.name_prevalence,
    append = TRUE
  )
}

# Artifact types ----

file.name_type <- "604828_site_type.md"

if (file.exists(file.name_type)) {
  file.remove(file.name_type)
}

for (s in Site$SITE_ID) {
  # cat(
  #   paste(
  #     "\n",
  #     "# Artifact types by Prevalence for ",
  #     filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
  #     "\n",
  #     sep = ""
  #   ),
  #   sep = "\n",
  #   file = file.name_type,
  #   append = TRUE
  # )
  # cat(
  print(
    dat %>%
      # filter(SITE_ID == s,!(CODE %in% code_exclude)) %>%
      filter(
        SITE_ID == s,
        !(CODE %in% code_exclude),
        CODE %in% code_precontact
      ) %>% #Pre-Contact Only
      # group_by(ASSEMBLAGE, TYPENAME, CATNAME) %>%
      group_by(TYPENAME, CATNAME, IDENTIFY) %>% #Pre-Contact only
      select(PROVENIENCE_ID, N, START_DEPTH, END_DEPTH) %>%
      summarise(
        n.provs = length(unique(PROVENIENCE_ID)),
        n.obj = sum(N),
        start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 2),
        end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 2),
        .groups = "drop_last"
      ) %>%
      # filter(frequency > 1) %>%
      mutate(art_type_sum = sum(n.provs)) %>%
      arrange(desc(art_type_sum), desc(n.provs), desc(n.obj), .by_group = FALSE) %>%
      # arrange(.by_group = FALSE) %>%
      select(-art_type_sum) %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        frequency,
        nisp
      ) %>%
      kable(
        format = 'markdown',
        row.names = FALSE,
        digits = 1,
        caption = paste(
          "Summary of prevalent artifact types for ",
          filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
          ".",
          sep = ""
        )
        #   ),
        # sep = "\n",
        # file = file.name_type,
        # append = TRUE
      )
  )
}

# Historical diagnostics ----

file.name_diag <- "604828_site_hist_diag.md"

if (file.exists(file.name_diag)) {
  file.remove(file.name_diag)
}

for (s in Site$SITE_ID) {
  # Check if site has diagnostics:
  if (nrow(filter(
    dat,
    SITE_ID == s,!(CODE %in% code_exclude),
    CODE %in% code_diagnostic
  )) > 0) {
    cat(
      paste(
        "\n",
        "# Diagnostic historical Artifacts for ",
        filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
        "\n",
        sep = ""
      ),
      sep = "\n",
      file = file.name_diag,
      append = TRUE
    )
    cat(
      dat %>%
        filter(
          SITE_ID == s,!(CODE %in% code_exclude),
          CODE %in% code_diagnostic
        ) %>%
        group_by(ASSEMBLAGE, IDENTIFY) %>%
        select(PROVENIENCE_ID, BEGDATE, ENDDATE, N, START_DEPTH, END_DEPTH) %>%
        summarise(
          BEGDATE = median(BEGDATE),
          ENDDATE = median(ENDDATE),
          n.provs = length(unique(PROVENIENCE_ID)),
          n.obj = sum(N),
          start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
          end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1)
        ) %>%
        arrange(desc(frequency), desc(NISP),
                .by_group = FALSE) %>%
        adorn_totals(
          .,
          where = "row",
          fill = "--",
          na.rm = TRUE,
          name = "Total",
          NISP
        ) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          caption = paste(
            "Historical diagnostic artifacts for ",
            filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
            ".",
            sep = ""
          )
        ),
      sep = "\n",
      file = file.name_diag,
      append = TRUE
    )
  }
}

# Pre-contact Artifacts ----

file.name_precontact <- "604828_site_precontact.md"

if (file.exists(file.name_precontact)) {
  file.remove(file.name_precontact)
}

for (s in Site$SITE_ID) {
  # Check if site has pre-contact:
  if (nrow(filter(dat,
                  SITE_ID == s, CODE %in% code_precontact)) > 0) {
    cat(
      paste(
        "\n",
        "# Pre-contact artifacts for ",
        filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
        "\n",
        sep = ""
      ),
      sep = "\n",
      file = file.name_precontact,
      append = TRUE
    )
    cat(
      dat %>%   filter(SITE_ID == s, CODE %in% code_precontact) %>%
        group_by(TYPENAME, CATNAME, IDENTIFY) %>%
        select(PROVENIENCE_ID, N, START_DEPTH, END_DEPTH) %>%
        summarise(
          n.provs = length(unique(PROVENIENCE_ID)),
          n.obj = sum(N),
          start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
          end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1)
        ) %>%
        mutate(art_type_sum = sum(frequency)) %>%
        arrange(desc(frequency), desc(NISP), .by_group = FALSE) %>%
        arrange(desc(art_type_sum), .by_group = FALSE) %>%
        select(-art_type_sum) %>%
        adorn_totals(
          .,
          where = "row",
          fill = "--",
          na.rm = TRUE,
          name = "Total",
          NISP
        ) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          caption = paste(
            "Pre-contact artifacts for ",
            filter(Site, SITE_ID == s) %>% pull(SITE_NAME),
            ".",
            sep = ""
          )
        ),
      sep = "\n",
      file = file.name_precontact,
      append = TRUE
    )
  }
}

# Artifacts by Stratum ----

## Assemblage by Stratum ----


for (s in Site$SITE_ID) {
  print(
    dat %>%
      filter(CODE %in% code_precontact) %>%
      group_by(STRATUM, TYPENAME) %>%
      summarise(
        n.provs = length(unique(PROVENIENCE_ID)),
        n.obj = sum(N),
        start.depth = round(weighted.mean(START_DEPTH, N, na.rm = TRUE), 1),
        end.depth = round(weighted.mean(END_DEPTH, N, na.rm = TRUE), 1),
        mean.depth = (start.depth + end.depth) / 2,
      ) %>%
      mutate(stratum.mean = mean(mean.depth)) %>%
      arrange(stratum.mean, mean.depth) %>%
      select(-mean.depth,-stratum.mean) %>%
      adorn_totals(
        .,
        where = "row",
        fill = "--",
        na.rm = TRUE,
        name = "Total",
        NISP
      ) %>%
      kable(
        format = "markdown",
        caption = "Artifacts by stratum.",
        row.names = FALSE,
        digits = 1
      )
  )
}


## Type by Stratum ----

for (s in Site$SITE_ID) {
  for (typ in distinct(filter(dat, SITE_ID == s,!(CODE %in% code_exclude)),
                       TYPENAME) %>% pull()) {
    print(
      dat %>%
        filter(TYPENAME == typ,!(CODE %in% code_exclude)) %>%
        group_by(STRATUM) %>%
        summarise(
          start.depth = weighted.mean(START_DEPTH, N, na.rm = TRUE),
          end.depth = weighted.mean(END_DEPTH, N, na.rm = TRUE),
          mean.depth = mean(start.depth, end.depth),
          n.provs = length(unique(PROVENIENCE_ID)),
          n.obj = sum(N)
        ) %>%
        arrange(mean.depth, desc(frequency), .by_group = FALSE) %>%
        select(-mean.depth) %>%
        kable(
          format = "markdown",
          caption = paste(typ, "by stratum."),
          row.names = FALSE,
          digits = 1
        )
    )
  }
}
