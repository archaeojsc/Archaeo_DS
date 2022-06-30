# Archaeological Data Science

A "day in the life" repository of using data science as an archaeologist...
samples from some of my work projects and experiments (mainly in R). Much of
this work revolves around simple ETL and EDA, summarization and descriptive
statistics, manipulating and summarizing spatial data, and a few more "advanced"
experiments with probability distribution modeling.

The code in this repository is very rough for now. I'll clean them up into
formal scripts as I have time.

Currently available:

* Summarization utility functions for artifact counts, find depths, and
  temporally diagnostic dates using `tidyverse`. ([R
  Script](code/summary_table_functions.R))
* Experiments with creating mixture distributions (using `mixdist`) for
  temporally diagnostic artifact dates comprised of the independent uniform
  distributions for each artifact type. Working on bootstrap estimates. ([R
  Script](code/Diagnostic_Date_Plotting.R))
* Example automation code for generating standardized summaries for different
  artifact categories and spatial units for use in reports, piped to markdown
  tables using `tidyverse`. Incorporate spatial data pulled form ESRI shapefiles
  with `terra` and `sf`. ([R Script](code/ProjectPhaseSummaryCode.R))
