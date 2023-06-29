#------------ Function for perform deviance analysis from remlf90 function -----#

Deviance_BreedR <-
  function(Trait = Trait,
           Model = breedR_model,
           Data = Progenies,
           Pedigree = NULL,
           Method = "EM",
           model_type = "std_animal") {
    if (sum(stringr::str_count(names(Model[["mf"]]), Trait)) < 1 &
        sum(stringr::str_count(names(Model[["mf"]]), "get")) < 1) {
      stop(paste(
        Trait,
        "trait not found. Please check the evaluated trait in the model"
      ))
    }
    if (!model_type %in% c("random", "std_animal")) {
      stop("model_type is not appropriate. Please, choose between random or std_animal model")
    }
    
    # Detect fixed effects
    fixed_effect <- names(Model[["fixed"]])
    # Detect random effects
    random_effect <- names(Model[["ranef"]])
    # random effect without genetic effect
    random_effect2 <-
      random_effect[str_detect(random_effect, "genetic", negate = T)]
    
    # Function to transform variables in factors
    to.factors <- function(df, variables) {
      for (variable in variables) {
        df[[variable]] <- as.factor(df[[variable]])
      }
      return(df)
    }
    
    # Convert effects in the model in factors
    Data <- to.factors(Data, c(fixed_effect, random_effect2))
    
    if (is.null(Pedigree)) {
      ped <<- get(Model[["call"]][["genetic"]][["pedigree"]])
    } else {
      ped <<- Pedigree
    }
    
    
    DEV <- list()
    MOD <- list()
    #FIT.2logL <- list()
    
    for (i in random_effect2) {
      # loop for model_type
      # if model_type == 'random':
      if (model_type == "random") {
        # Update models
        Model2 <- remlf90(
          as.formula(Model[["call"]][["fixed"]]),
          random = update.formula(Model[["call"]][["random"]],
                                  formula(paste("~ . ", "-", i))),
          dat = Data,
          method = Method
        )
        
      } else {
        # if model_type == 'animal':
        if (sum(stringr::str_count(names(Model[["mf"]]), Trait)) < 1 &
            sum(stringr::str_count(names(Model[["mf"]]), "get")) > 0) {
          model_fix <-  paste(Trait, paste("~", sub("get(+$)", "", Model[["call"]][["fixed"]])[3]), collapse = "")
        } else {
          model_fix <- Model[["call"]][["fixed"]]
        }
        
        
        Model2 <- remlf90(
          as.formula(model_fix),
          random = update.formula(Model[["call"]][["random"]],
                                  formula(paste("~ . ", "-", i))),
          genetic = list(
            model = "add_animal",
            pedigree = ped,
            id = Model[["call"]][["genetic"]][["id"]]
          ),
          dat = Data,
          method = Method
        )
        # Model without genetic effect
        
        Model3 <- remlf90(
          as.formula(model_fix),
          random = as.formula(Model[["call"]][["random"]]),
          dat = Data,
          method = Method
        )
      }
      # Extract -2logL
      TestStat <-
        Model2[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]]
      # Extract parameters differences between models
      DifGL <-
        length(names(Model[["ranef"]])) - length(names(Model2[["ranef"]]))
      # Obtain the p-value
      p.value <- pchisq(TestStat, df = DifGL, lower.tail = FALSE)
      # Building the deviance table
      DataDev <-
        data.frame(Efeito = i,
                   LRT = TestStat,
                   p.valor = p.value)
      # FIT.2logL[[i]] <- Model2[["fit"]][["-2logL"]]
      DEV[[i]] <- DataDev
      
    }
    Deviance <- bind_rows(DEV)
    
    DevGenetic <-
      c(
        "genetic",
        Model3[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]],
        pchisq(
          Model3[["fit"]][["-2logL"]] - Model[["fit"]][["-2logL"]],
          df = length(names(Model[["ranef"]])) - length(names(Model3[["ranef"]])),
          lower.tail = FALSE
        )
      )
    #DF <- length(names(Model[["ranef"]])) - length(names(Model3[["ranef"]]))
    Deviance[3, ] <- DevGenetic
    return(Deviance)
  }