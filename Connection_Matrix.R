#' Calculate the Connection Matrix Between Locations Based on Common genotypes
#'
#' This function calculates a connection matrix between environments (locations)
#' based on the number of genotype common to each pair of environments. It accepts 
#' a dataset that contains information on locations, genotype, and additional factors
#' related to the experimental design.
#'
#' @param data A data frame containing the experimental data, with at least one column 
#' representing locations and one column representing genotypes.
#' @param location_col A character string specifying the name of the column that contains the locations.
#' @param gen_col A character string specifying the name of the column that contains the genotypes.
#' @return A matrix with the number of common genotypes between each pair of locations.
#' @importFrom reshape2 melt
#' @examples
#' set.seed(123)
#' data_example <- expand.grid(
#'   Local = paste0("L", 1:5), 
#'   Clone = paste0("C", 1:10), 
#'   Repetition = 1:5, 
#'   Plot = 1:3, 
#'   Plant = 1:5, stringsAsFactors = FALSE
#' )
#' # Randomly removing some clones from certain locations to create imbalance
#' data_example_desbalanceado <- data_example |> 
#'   group_by(Local, Clone) |> 
#'   filter(runif(1) > 0.2) |>  # 20% chance of removing a clone from a location
#'   ungroup()
#' # Calculate the connection matrix
#' connection_matrix <- calculate_connection_matrix(data_example_desbalanceado, "Local", "Clone")
#' print(connection_matrix)
calculate_connection_matrix <- function(data, location_col, gen_col) {
  
  # Check if the 'reshape2' package is installed, if not, install it
  if (!require("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
  }
  
  # Extract the unique locations from the specified location column
  locations <- unique(data[[location_col]])
  
  # Internal function to calculate the number of genotypes common between two locations
  common_clones_between_locations <- function(local1, local2) {
    gen_local1 <- unique(data[data[[location_col]] == local1, gen_col, drop = TRUE])  # Vector of genotype in local1
    gen_local2 <- unique(data[data[[location_col]] == local2, gen_col, drop = TRUE])  # Vector of genotype in local2
    length(intersect(gen_local1, gen_local2))  # Number of common genotype
  }
  
  # Create an empty connection matrix
  connection_matrix <- matrix(0, nrow = length(locations), ncol = length(locations),
                              dimnames = list(locations, locations))
  
  # Fill the matrix with the number of common clones between each pair of locations
  for (i in 1:length(locations)) {
    for (j in i:length(locations)) {
      connection_matrix[i, j] <- common_clones_between_locations(locations[i], locations[j])
      connection_matrix[j, i] <- connection_matrix[i, j]  # The matrix is symmetric
    }
  }
  
  return(connection_matrix)
}
