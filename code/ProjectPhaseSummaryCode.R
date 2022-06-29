require(tmap)

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

cat(
  paste(
    "\n",
    "## Finds by Assemblages",
    "\n",
    sep = ""
  ),
  sep = "\n",
  file = file.name_assemblages,
  append = TRUE
)

for (a in unique(filter(dat,!(CODE %in% code_exclude)) %>%
                 pull(ASSEMBLAGE))) {
  if (nrow(filter(dat,
                  ASSEMBLAGE == a, !(CODE %in% code_exclude))) > 0) {
    cat("\n",
        file = file.name_assemblages,
        append = TRUE)
    
    cat(
      dat %>%
        filter(ASSEMBLAGE == a ,!(CODE %in% code_exclude)) %>%
        group_by(TYPENAME, CATNAME) %>%
        summary_count_depth(total = TRUE) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          caption = paste(
            "Summary of ",
            a,
            " assemblage artifact types and categories.",
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


## Assemblages by project phase ----

cat(
  paste0("\n",
         "# Artifact assemblages by project phase"
  ),
  sep = "\n",
  file = file.name_assemblages,
  append = TRUE
)


for (p in unique(dat$PHASE)) {
  if (nrow(filter(dat, PHASE == p, !(CODE %in% code_exclude))) > 0) {
    cat("\n", file = file.name_assemblages, append = TRUE)
    cat(
      dat %>%
        filter(PHASE == p, !(CODE %in% code_exclude)) %>%
        group_by(ASSEMBLAGE) %>%
        summary_count_depth(total = TRUE) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          align = "lrrrrr",
          caption = paste(
            "Artifact assemblages for ",
            unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
            ".",
            sep = ""
          )
        ),
      sep = "\n",
      file = file.name_assemblages,
      append = TRUE
    )
  }
}

# Artifact prevalence -----------------------------------------------------

file.name_prevalence <- "604828_site_prevalence.md"

if (file.exists(file.name_prevalence)) {
  file.remove(file.name_prevalence)
}

## Artifact prevalence for all project phases -----------------------------

cat(
  paste0("\n",
         "# Artifact prevalence",
         "\n"),
  sep = "\n",
  file = file.name_prevalence,
  append = TRUE
)

cat(
  dat %>%
    filter(!(CODE %in% code_exclude)) %>%
    group_by(ASSEMBLAGE, TYPENAME, CATNAME) %>%
    summary_count_depth(total = TRUE) %>%
    filter(n.units > 1) %>%
    kable(
      format = 'markdown',
      row.names = FALSE,
      digits = 1,
      caption ="Artifact prevalence, all project phases.",
      align = "lllrrrrr"
    ), 
  sep = "\n",
  file = file.name_prevalence,
  append = TRUE
)

## Artifact prevalence by project phase ----

cat(
  paste0("\n",
         "## Artifact prevalence by project phase"
  ),
  sep = "\n",
  file = file.name_prevalence,
  append = TRUE
)


for (p in unique(dat$PHASE)) {
  if (nrow(filter(dat, PHASE == p, !(CODE %in% code_exclude))) > 0) {
    cat("\n", file = file.name_prevalence, append = TRUE)
    cat(
      dat %>%
        filter(PHASE == p, !(CODE %in% code_exclude)) %>%
        group_by(ASSEMBLAGE, TYPENAME, CATNAME) %>%
        summary_count_depth(total = TRUE) %>%
        filter(n.units > 1) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          align = "lllrrrrr",
          caption = paste(
            "Artifact prevalence for ",
            unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
            ".",
            sep = ""
          )
        ),
      sep = "\n",
      file = file.name_prevalence,
      append = TRUE
    )
  }
}

# Precontact artifacts -------------------------------------

file.name_precontact <- "604828_site_precontact.md"

if (file.exists(file.name_precontact)) {
  file.remove(file.name_precontact)
}

## Precontact artifacts for all project phases -----------------------------

cat(
  paste0("\n",
         "# Pre-Contact artifacts",
         "\n"),
  sep = "\n",
  file = file.name_precontact,
  append = TRUE
)

cat(
  dat %>%
    filter(CODE %in% code_precontact, !(CODE %in% code_exclude)) %>%
    group_by(CATNAME, IDENTIFY) %>%
    summary_count_depth(total = TRUE) %>%
    # filter(n.units > 1) %>%
    kable(
      format = 'markdown',
      row.names = FALSE,
      digits = 1,
      caption ="Pre-Contact artifacts, all project phases.",
      align = "llrrrrr"
    ), 
  sep = "\n",
  file = file.name_precontact,
  append = TRUE
)

## Precontact artifacts by project phase ----

cat(
  paste0("\n",
         "## Pre-Contact artifacts by project phase"
  ),
  sep = "\n",
  file = file.name_precontact,
  append = TRUE
)


for (p in unique(dat$PHASE)) {
  if (nrow(filter(
    dat,
    PHASE == p,
    CODE %in% code_precontact,
    !(CODE %in% code_exclude)
  )) > 0) {
    cat("\n", file = file.name_precontact, append = TRUE)
    cat(
      dat %>%
        filter(
          PHASE == p,
          CODE %in% code_precontact,
          !(CODE %in% code_exclude)
        ) %>%
        group_by(CATNAME, IDENTIFY) %>%
        summary_count_depth(total = TRUE) %>%
        # filter(n.units > 1) %>%
        kable(
          format = 'markdown',
          row.names = FALSE,
          digits = 1,
          align = "llrrrrr",
          caption = paste(
            "Pre-Contact artifacts for ",
            unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
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

# Assemblages by project phase ----


# Artifact summaries by assemblage for each site

for (p in unique(dat$PHASE)) {
  cat(
    paste(
      "\n",
      "# Finds by Assemblages for ",
      unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
      "\n",
      sep = ""
    ),
    sep = "\n",
    file = file.name_assemblages,
    append = TRUE
  )
  
  cat(
    dat %>%
      filter(PHASE == p,!(CODE %in% code_exclude)) %>%
      group_by(ASSEMBLAGE) %>%
      summary_count_depth(total = TRUE) %>%
      kable(
        format = 'markdown',
        row.names = FALSE,
        digits = 1,
        caption = paste(
          "Summary of assemblages for ",
          unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
          ".",
          sep = ""
        ),
        align = "lrrrr"
      ),
    sep = "\n",
    file = file.name_assemblages,
    append = TRUE
  )
  
  for (a in unique(filter(dat, PHASE == 1, 
                          !(CODE %in% code_exclude)) %>% 
                   pull(ASSEMBLAGE))) {
    if (nrow(filter(dat,
                    PHASE == p,
                    ASSEMBLAGE == a,!(CODE %in% code_exclude))) > 0) {
      cat("\n",
          file = file.name_assemblages,
          append = TRUE)
      
      cat(
        dat %>%
          filter(PHASE == p, ASSEMBLAGE == a , !(CODE %in% code_exclude)) %>%
          group_by(TYPENAME, CATNAME) %>%
          summary_count_depth(total = TRUE) %>%
          kable(
            format = 'markdown',
            row.names = FALSE,
            digits = 1,
            caption = paste(
              "Summary of ",
              a,
              " assemblage artifact types and categories for ",
              unique(filter(dat, PHASE == p) %>% pull(PROJECT_PHASE)),
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
}

# Pre-Contact Pivot tables ----

unit.precontact <-  dat %>%
  filter(ASSEMBLAGE == "Pre-Contact", !(CODE %in% code_exclude)) %>%
  group_by(PHASE,
           PROJECT_PHASE,
           UNIT_TYPE,
           UNIT_LABEL,
           TYPENAME) %>%
  summarise(nisp = sum(N)) %>%
  pivot_wider(names_from = TYPENAME,
              values_from = nisp,
              values_fill = 0) %>%
  adorn_totals("col",
               name = "Precontact Total") %>%
  clean_names("all_caps")

unit.debitage <- dat %>%
  filter(TYPENAME == "Debitage", !(CODE %in% code_exclude)) %>%
  group_by(PHASE, PROJECT_PHASE, UNIT_TYPE, UNIT_LABEL, CATNAME) %>%
  summarise(nisp = sum(N)) %>%
  pivot_wider(names_from = CATNAME,
              values_from = nisp,
              values_fill = 0) %>%
  adorn_totals(where = "col", name = "Total Debitage") %>%
  clean_names("all_caps")

unit.debitage.spatial <-
  spatial_Test_loc %>% inner_join(
    unit.debitage,
    by = c(
      'P_TYPE' = 'PROJECT_PHASE' ,
      'U_TYPE' = 'UNIT_TYPE',
      'U_LABEL' = 'UNIT_LABEL'
    )
  ) %>% select(-fid_1, -layer, -path,-POINT_X,-POINT_Y)


unit.debitage.flake <- dat %>%
  filter(CATNAME == "Flake",!(CODE %in% code_exclude)) %>%
  group_by(PHASE, PROJECT_PHASE, UNIT_TYPE, UNIT_LABEL, IDENTIFY) %>%
  summarise(nisp = sum(N)) %>%
  pivot_wider(names_from = IDENTIFY,
              values_from = nisp,
              values_fill = 0) %>%
  adorn_totals("col", name = "Total Flake") %>%
  clean_names("all_caps")


# Assemblage pivot tables ----

unit.assemblages <- dat %>%
  filter(!(CODE %in% code_exclude)) %>%
  group_by(PHASE, PROJECT_PHASE, UNIT_TYPE, UNIT_LABEL, ASSEMBLAGE) %>%
  summarise(nisp = sum(N)) %>%
  pivot_wider(names_from = ASSEMBLAGE,
              values_from = nisp,
              values_fill = 0) %>%
  adorn_totals("col", name = "Total CM") %>%
  clean_names("all_caps")

unit.assemblages.spatial <-
  spatial_Test_loc %>% inner_join(
    unit.assemblages,
    by = c(
      'PHASE' = 'PROJECT_PHASE' ,
      'U_TYPE' = 'UNIT_TYPE',
      'U_LABEL' = 'UNIT_LABEL'
    )
  ) %>% select(-fid_1, -layer, -path,-POINT_X,-POINT_Y)



# Distribution bubble maps ----

## Increase map's bounding box default size ---- 
# from https://www.jla-data.net/eng/adjusting-bounding-box-of-a-tmap-map/

ape_bbox <- st_bbox(spatial_APE)
bb_xrange <- ape_bbox$xmax - ape_bbox$xmin
bb_yrange <- ape_bbox$ymax - ape_bbox$ymin

ape_bbox[1] <- ape_bbox[1] - (1.25 * bb_xrange) # xmin - left
ape_bbox[3] <- ape_bbox[3] + (0.25 * bb_xrange) # xmax - right
# ape_bbox[2] <- ape_bbox[2] - (0.5 * bb_yrange) # ymin - bottom
# ape_bbox[4] <- ape_bbox[4] + (0.5 * bb_yrange) # ymax - top

ape_bbox <- ape_bbox %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon

## Background orthoimage ---

ortho <- terra::rast("./maps/604828_ortho.tif")

# set projection
crs(ortho) <- "epsg:32117"

ortho_crop <- crop(ortho, ape_bbox)


## Pre-Contact assemblage ----

map_distPC <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.assemblages.spatial) +
  tm_bubbles(
    size = "PRE_CONTACT",
    scale = sqrt(2),
    size.lim = c(0, 50),
    alpha = 0.5,
    col = "red",
    title.size = "Pre-Contact"
  )

## Domestic assemblage ----
map_distDOM <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.assemblages.spatial) +
  tm_bubbles(
    size = "DOMESTIC",
    scale = sqrt(2),
    size.lim = c(0, 50),
    alpha = 0.3,
    col = "blue",
    title.size = "Domestic"
  )

## Architectural assemblage ----

map_distARC <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.assemblages.spatial) +
  tm_bubbles(
    size = "ARCHITECTURAL",
    scale = sqrt(2),
    size.lim = c(0, 50),
    alpha = 0.3,
    col = "darkolivegreen",
    title.size = "Architectural"
  )

## Other assemblage ----

map_distOTH <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.assemblages.spatial) +
  tm_bubbles(
    size = "OTHER",
    scale = sqrt(2),
    # size.lim = c(0, 50),
    alpha = 0.3,
    col = "yellow",
    title.size = "Other"
  )

map_dist <-
  tmap_arrange(map_distPC, map_distDOM, map_distARC, map_distOTH, nrow = 2)

tmap_save(tm = map_dist, "./report/604828_assemblage_map.png", dpi = 300)

## Flake distribution ----

map_distFLAKE <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.debitage.spatial) +
  tm_bubbles(
    size = "FLAKE",
    scale = sqrt(2),
    size.lim = c(0, 50),
    alpha = 0.5,
    col = "red",
    title.size = "Flakes"
  )

map_distDEBITAGE <-
  tm_shape(ortho, bbox = ape_bbox) +
  tm_rgb(alpha = 0.4) +
  tm_shape(shp = spatial_APE,
           name = "APE",
           unit = "m", ) +
  tm_lines(col = "red", lwd = 1.5, lty = 2) +
  tm_shape(unit.debitage.spatial) +
  tm_bubbles(
    size = "TOTAL_DEBITAGE",
    scale = sqrt(2),
    size.lim = c(0, 50),
    alpha = 0.5,
    col = "red",
    title.size = "Lithic Debitage"
  )
