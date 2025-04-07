#' Fit Linear Mixed-Effects Models for Multiple Trials and Extract Genetic Parameters
#'
#' This function fits linear mixed-effects models for multiple trials and extracts genetic parameters.
#'
#' @param data A data frame containing the data to be analyzed. Default is `DataDiag`.
#' @param model A formula specifying the model to be fitted. Default is `"VOL36 ~ REP + (1 | PROG) + (1 | PARC)"`.
#' @param factor_for_mean A character string specifying the factor for calculating the overall mean. Default is `"TRIAL"`.
#' @param genetic_factor A character string specifying the genetic factor. Default is `"PROG"`.
#' @param Plot_column A character string specifying the plot column. Default is `"PARC"`.
#' @param Rep_column A character string specifying the replicate column. Default is `"REP"`.
#' @param Trial_column A character string specifying the trial column. Default is `"TRIAL"`.
#' @param random A character string specifying additional random effects. Default is `NULL`.
#' @param survival_calc A logical value indicating whether to calculate survival. Default is `TRUE`.
#' @param survival_column A character string specifying the survival column. Default is `"SR36"`.
#' @param survival_values A numeric vector specifying the survival values. Default is `c(0, 1)`.
#' @return A list of results for each trial, including genetic parameters and BLUPs.
#' @examples
#' \dontrun{
#' results <- lmer_func(data = my_data, model = "VOL36 ~ REP + (1 | PROG) + (1 | PARC)")
#' }
#' @export
# Function to fit the model for multiple trials and extract genetic parameters
lmer_func <- function(data = DataDiag,
                      model = "VOL36 ~ REP + (1 | PROG) + (1 | PARC)",
                      factor_for_mean = "TRIAL",
                      genetic_factor = "PROG",
                      Plot_column = "PARC",
                      Rep_column = "REP",
                      Trial_column = "TRIAL",
                      random = NULL,
                      survival_calc = TRUE,
                      survival_column = "SR36",
                      survival_values = c(0, 1),
                      is_progeny_trial = FALSE) {
  # Prepare data
  ## Clean the string
  cleaned_str <- gsub("\\)", "", gsub("\\(1 \\| ", "", model))
  ## get the terms
  terms <- unlist(strsplit(cleaned_str, "[~+]"))
  ## Trim whitespace and remove response variable (first element)
  terms  <- trimws(terms[-1])
  
  # Fit the model
  ## Formula
  model_formula = as.formula(model)
  
  # Get a loop for multiple trials
  ## Get trial names
  trial_names <- unique(data[[Trial_column]])
  Results <- list()
  # Make all variables as character
  data <- data |> 
    mutate_at(
      all_of(terms), as.character
    )
  for (i in trial_names) {
    #i = "EP03EEAA"
    # Fit model
    mod_lmer <- lmer(model_formula, 
                     data = data |> 
                       # Filter by trial
                       filter(
                         get(Trial_column) == i
                       ) |> 
                       mutate(
                         # Convert terms to factors
                         across(
                           all_of(terms), as.factor
                         )
                       ))
    # Extract trait
    trait <- trimws(gsub("~.*", "", model))
    # Obtain the overall mean
    overall_mean <- mean(data[[trait]], na.rm = T)
    # Extract genetic parameter
    varcomp_lmer <- summary(mod_lmer)$varcor |>
      data.frame() |>
      column_to_rownames(var = "grp")
    
    n_obs <- data |>
      filter(get(Trial_column) == i) |> 
      group_by(get(genetic_factor)) |>
      summarise(n = sum(!is.na(get(trait))))
    # Load psych package to calculate the harmonic mean
    require(psych)
    n_obs <-  psych::harmonic.mean(n_obs$n, zero = FALSE)
    
    #weighted.mean(n_obs$n)
    # Extract variance components
    VarComp <- as.data.frame(VarCorr(mod_lmer)) %>%
      dplyr::select(!c("var1", "var2")) %>%
      `colnames<-`(c("FV", "Variance", "Standard Deviation")) %>%
      `rownames<-`(., .$FV)
    
    Plot <- ifelse(is.null(Plot_column), NA, Plot_column)
    nRep <- length(unique(subset(data, get(Trial_column) == i)[[Rep_column]]))
    
    
    if (is_progeny_trial) {
      GenPar <- tibble(
        VFam = VarComp[genetic_factor, "Variance"],
        Vplot = ifelse(!is.null(Plot_column), VarComp[Plot, "Variance"], Plot),
        Vrandom = ifelse(is.null(random), NA, VarComp[random, "Variance"]),
        Vres = VarComp["Residual", "Variance"],
        Vpheno = sum(VFam, Vplot, Vrandom, Vres, na.rm = T),
        h2Fam_mean = VFam / sum(VFam, (Vplot / nRep), (Vres / n_obs), na.rm = T),
        h2a = (4 * VFam) / ( 4 * VFam + sum(Vplot, Vrandom, Vres, na.rm = T)),
        c2p = ifelse(is.na(Vplot), NA, Vplot / Vpheno)
      ) %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`("Estimates") %>%
        rownames_to_column(var = "Components")
    } else {
      GenPar <- tibble(
        Vg = VarComp[genetic_factor, "Variance"],
        Vplot = ifelse(!is.null(Plot_column), VarComp[Plot, "Variance"], Plot),
        Vrandom = ifelse(is.null(random), NA, VarComp[random, "Variance"]),
        Vres = VarComp["Residual", "Variance"],
        Vpheno = sum(Vg, Vplot, Vrandom, Vres, na.rm = T),
        `h2g` = Vg / Vpheno,
        `h2m` = Vg / sum(Vg, (Vplot / nRep), (Vres / n_obs), na.rm = T),
        `c2p` = ifelse(is.na(Vplot), NA, Vplot / Vpheno)
      ) %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`("Estimates") %>%
        rownames_to_column(var = "Components")
    }
    
    ## Extract standard errors for BLUPs
    se_blups <- arm::se.ranef(mod_lmer)[[genetic_factor]] |>
      as.data.frame() |>
      rownames_to_column(genetic_factor) |>
      rename(SE = `(Intercept)`)
    
    ## Extract BLUPs
    blups <- ranef(mod_lmer)[[genetic_factor]] |>
      `colnames<-`("solution") |>
      data.frame() |>
      rownames_to_column(genetic_factor)
    ## Combine BLUPs and SEs
    blups_with_se <- left_join(blups, se_blups, by = genetic_factor)
    ## Extract blups
    if (is_progeny_trial) {
      blups <- blups_with_se |>
        mutate(
          # Calculate the new mean BV + adjusted mean
          `BLUP_plus_mean (a+u)` = solution + overall_mean,
          # calculate PEV
          PEV = SE^2,
          # Calculate accuracy
          accuracy = sqrt(1 - PEV / GenPar[which(GenPar[["Components"]] == "VFam"), "Estimates"]),
          `overall_mean (u)` = overall_mean
        ) |>
        rename(`BLUP (a)` = solution) |>
        mutate_at(all_of(genetic_factor), as.character())
    } else {
      blups <- blups_with_se |>
        mutate(
          # Calculate the new mean BV + adjusted mean
          `BLUP_plus_mean (g+u)` = solution + overall_mean,
          # calculate PEV
          PEV = SE^2,
          # Calculate accuracy
          accuracy = sqrt(1 - PEV / GenPar[which(GenPar[["Components"]] == "Vg"), "Estimates"]),
          `overall_mean (u)` = overall_mean
        ) |>
        rename(`BLUP (g)` = solution) |>
        mutate_at(all_of(genetic_factor), as.character())
    }
    
    if (survival_calc) {
      # get the survival of clones
      survival <- data |>
        filter(get(Trial_column) == i) |> 
        group_by(get(genetic_factor)) |>
        summarise(
          Dead = sum(get(survival_column) == survival_values[1]),
          Alive = sum(get(survival_column) == survival_values[2]),
          .groups = "drop"
        ) |>
        mutate(Survival = round(Alive / (Dead + Alive), 2)) |>
        dplyr::select(-c(Dead, Alive)) |>
        `colnames<-`(c(genetic_factor, "Survival"))
      # Merge with blups
      if (is_progeny_trial) {
        blups <- blups |>
          left_join(survival, by = genetic_factor) |> 
          arrange(desc(`BLUP (a)`))
      } else {
        blups <- blups |>
          left_join(survival, by = genetic_factor) |> 
          arrange(desc(`BLUP (g)`)) 
      }
    }
    Results[[i]] <- list(
      GenPar = GenPar,
      blups = blups
    )
    
  }
  
  return(
    Results
  )
}
