#--------------------------- Function for get breeding values ------------------#

BV_BreedR <- function(data = data,
                      Model = NULL,
                      NamesID = "ID2",
                      Pedigree = Pedigree,
                      dominanceName = NULL,
                      ColumnFamilies = "Prog",
                      Rank = T,
                      NamesReturnDataset = c("ID2", "Mae", "Bloco", "Arvore", "DAP", "ALT"),
                      model_type = "std",
                      random_effect = NULL,
                      fixed_effect_av = "Intercept") {
  if (class(Model)[1] != "breedR") {
    stop("This is not a breedR adjusted model")
  }
  
  # Map original individuals' ID
  MapPed <- attr(Pedigree, "map")
  
  # Get the original individuals' ID
  labelPed <- match(get_pedigree(Model)@label, MapPed)
  nParents <-
    as.matrix(summary(as.data.frame(get_pedigree(Model))$dam))["Max.", 1]
  nId <-
    as.matrix(summary(as.data.frame(get_pedigree(Model))$self))["Max.", 1] - nParents
  # Average
  mu <- mean(fixef(Model)[[fixed_effect_av]])
  
  # Standard and AR1 model
  if (model_type == "std" | model_type == "ar1") {
    if (model_type == "std" & !"genetic" %in% names(ranef(Model))) {
      stop("model_type is not 'std' model, please, set properly!")
    }
    if (model_type == "ar1" & !"genetic" %in% names(ranef(Model))) {
      stop("model_type is not 'ar1' model, please, set properly!")
    }
    # BV Parents
    BV_Parents <-
      Model[["ranef"]][["genetic"]][[1]][1:nParents, ] |>
      `rownames<-`(labelPed[1:nParents]) |>
      mutate(`u+a` = mu + value) |>
      rownames_to_column(var = "Family") |>
      `colnames<-`(c("Family", "a", "s.e", "u+a"))
    # BV Individuals
    BV_Ind <-
      Model[["ranef"]][["genetic"]][[1]][nParents + 1:nId, ] |>
      `rownames<-`(labelPed[nParents + 1:nId]) |>
      mutate(`u+a` = mu + value) |>
      rownames_to_column(var = NamesID) |>
      `colnames<-`(c(NamesID, "a", "s.e", "u+a"))
    # Genetic variance
    Va <- summary(Model)$var["genetic", "Estimated variances"]
    
    if (is.null(random_effect)) {
      BV_random_effect = NULL
    } else {
      # BV random effect
      BV_random_effect <- Model[["ranef"]][[random_effect]][[1]]
    }
  } else {
    # Comp and AR1-Comp model
    if (model_type == "Comp" | model_type == "ar1Comp") {
      if (model_type == "Comp" &
          !"genetic_direct" %in% names(ranef(Model))) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if (model_type == "ar1" &
          !"genetic_direct" %in% names(ranef(Model))) {
        stop("model_type is not 'ar1Comp' model, please, set properly!")
      }
      # Direct effect
      # BV Parents
      BV_Parents <-
        Model[["ranef"]][["genetic_direct"]][[1]][1:nParents, ] |>
        `rownames<-`(labelPed[1:nParents]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = "Family") |>
        `colnames<-`(c("Family", "a", "s.e", "u+a"))
      # BV Individuals
      BV_Ind <-
        Model[["ranef"]][["genetic_direct"]][[1]][nParents + 1:nId, ] |>
        `rownames<-`(labelPed[nParents + 1:nId]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = NamesID) |>
        `colnames<-`(c(NamesID, "a", "s.e", "u+a"))
      
      # Competition effect
      BV_Competition_Parents <-
        Model[["ranef"]][["genetic_competition"]][[1]][1:nParents, ] |>
        `rownames<-`(labelPed[1:nParents]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = "Family") |>
        `colnames<-`(c(
          "Family",
          "competition",
          "s.e_competition",
          "u+competition"
        ))
      # BV Individuals
      BV_Competition_Ind <-
        Model[["ranef"]][["genetic_competition"]][[1]][nParents + 1:nId, ] |>
        `rownames<-`(labelPed[nParents + 1:nId]) |>
        mutate(`u+a` = mu + value) |>
        rownames_to_column(var = NamesID) |>
        `colnames<-`(c(
          NamesID,
          "competition",
          "s.e_competition",
          "u+competition"
        ))
      
      # Genetic variance
      Va <-
        summary(Model)$var$genetic["genetic_direct", "genetic_direct"]
      if (is.null(random_effect)) {
        BV_random_effect = NULL
      } else {
        # BV random effect
        BV_random_effect <- Model[["ranef"]][[random_effect]][[1]]
      }
    }
  }
  # Accuracy from family model
  Accuracy.Family = mean(sqrt(1 - BV_Parents$s.e ^ 2 / Va))
  # Accuracy from ID model
  Accuracy.Ind <- mean(sqrt(1 - BV_Ind$s.e ^ 2 / Va), na.rm = T)
  # Change BV_Ind class
  class(BV_Ind[[NamesID]]) <- class(data[, NamesID])
  Data_Total <- left_join(data[, NamesReturnDataset],
                          BV_Ind, by = NamesID)
  
  # attr(BV_random_effect, "class") <- c("data.frame", "BV_BreedR")
  # attr(BV_Parents, "class") <- c("data.frame", "BV_BreedR")
  # attr(BV_Ind, "class") <- c("data.frame", "BV_BreedR")
  # attr(Data_Total, "class") <- c("data.frame", "BV_BreedR")
  #
  Results <- list(
    BV_random_effect = BV_random_effect,
    BV_Parents = BV_Parents,
    BV_Ind = BV_Ind,
    Data_Total = Data_Total,
    Accuracy.Family = Accuracy.Family,
    Accuracy.Ind = Accuracy.Ind
  )
  return(Results)
}