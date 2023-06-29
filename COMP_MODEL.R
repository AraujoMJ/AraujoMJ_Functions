#---------------------------- Function for model comparison --------------------------------------#
## 5.1. AIC and -2LogLik values
COMP_MODEL <- function(Traits = c("DAP", "ALT", "VOL"),
                       #Site_Column = "Local",
                       #Age_Column = "Idade",
                       Age = NULL,
                       # Models
                       standard_model = NULL,
                       spatial_model = NULL,
                       competition_model = NULL,
                       competition_spatial_model = NULL) {
  if (is.null(Age)) {
    standard_model2 <- list()
    spatial_model2 <- list()
    competition_model2 <- list()
    competition_spatial_model2 <- list()
    for (i in names(standard_model)) {
      age = "1"
      # Creating lists
      
      # Filling lists
      standard_model2[[i]][[age]] <- standard_model[[i]]
      spatial_model2[[i]][[age]] <- spatial_model[[i]]
      competition_model2[[i]][[age]] <- competition_model[[i]]
      competition_spatial_model2[[i]][[age]] <-
        competition_spatial_model[[i]]
    }
    # Replace lists
    standard_model <- standard_model2
    spatial_model <- spatial_model2
    competition_model <- competition_model2
    competition_spatial_model <- competition_spatial_model2
    # Remove lists
    rm(
      standard_model2,
      spatial_model2,
      competition_model2,
      competition_spatial_model2
    )
    Age = age
    
  }
  # Names models
  models <-
    c(
      "standard_model",
      "spatial_model",
      "competition_model",
      "competition_spatial_model"
    )
  # Create list with models
  MODELS <- as.list(mget(models))
  
  # if (class(Model)[1] != "breedR") {
  #       stop("This is not a breedR adjusted model")
  # }
  
  # Remove non-existent models
  TrueModel <- lapply(MODELS, function(x)
    !is.null(unlist(x)))
  TrueModel <- names(TrueModel[which(TrueModel == T)])
  
  # Extract
  AIC_LIST <- list()
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        AIC_LIST[[i]][[j]][[k]] <- list(
          AIC = MODELS[[i]][[j]][[k]][["fit"]][["AIC"]],
          Loglik = MODELS[[i]][[j]][[k]][["fit"]][["-2logL"]],
          Df = summary(MODELS[[i]][[j]][[k]])$var |>
            unlist() |>
            unique() |>
            length()
        )
      }
    }
  }
  
  DIF_LogLik <- list()
  Df_MODELS <- list()
  for (i in Traits) {
    for (j in Age) {
      # Extract Loglik
      DIF_LogLik[[i]][[j]] <-
        lapply(AIC_LIST, function(x)
          x[[i]][[j]][["Loglik"]]) |>
        unlist(use.names = T) %>%
        outer(., ., FUN = "-") |>
        abs()
      # Extract Degree of freedom
      Df_MODELS <-
        lapply(AIC_LIST, function(x)
          x[[i]][[j]][["Df"]]) |>
        unlist()
    }
  }
  
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        # LRT value
        AIC_LIST[[i]][[j]][[k]][["LRT_value"]] <-
          DIF_LogLik[[j]][[k]]["standard_model", ][i]
        # LRT p.value
        AIC_LIST[[i]][[j]][[k]][["LRT_p.value"]] <-
          pchisq(
            DIF_LogLik[[j]][[k]]["standard_model", ][i],
            df = abs(Df_MODELS["standard_model"] - Df_MODELS[i]),
            lower.tail = F
          )
      }
    }
  }
  
  MODEL_COMPARISON <- list()
  for (i in TrueModel) {
    for (j in Traits) {
      for (k in Age) {
        MODEL_COMPARISON[[i]][[j]][[k]] <-
          data.frame(AIC_LIST[[i]][[j]][[k]]) |>
          rownames_to_column(var = "Model") |>
          mutate(Traits = j,
                 Age = k)
      }
    }
  }
  
  #rm(M_LAPPLY)
  M_LAPPLY <- lapply(MODEL_COMPARISON, function(x)
    lapply(x, function(y)
      bind_rows(y))) |>
    lapply(function(x)
      bind_rows(x)) |>
    bind_rows()
  
  return(M_LAPPLY)
}
