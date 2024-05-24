#### Extracts a stable snapshot from the database to make the project independent ####

library(magrittr)

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

leaching_get_data_database(
  default.file = "~/my.cnf",
  file = "derived_data/dpeatdecomposition_snapshot.rds"
)
