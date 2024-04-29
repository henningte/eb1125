#### Extracts a stable snapshot from the database to make the project independent ####

library(magrittr)

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)

leaching_get_data_database(
  username = "root", # ---todo: adjust or get from config file
  password = "coucou",
  host = "mariadb",
  file = "derived_data/dpeatdecomposition_snapshot.rds"
)
