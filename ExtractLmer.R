#--------------------------- Function for extract estimates from lmer model ---------------#
# Function to fit the model and extract genetic parameters
  Extract_lmer <- function(data = DataDiag |>
                             filter(TRIAL == "EP03EBWK"),
                           factor_for_adj_mean = NULL,
                           genetic_factor = "PROG",
                           Plot_column = NULL,
                           Rep_column = "REP",
                           random = NULL,
                           trait = "VOL36",
                           mod_lmer = mod_lmer_EP03EBWK,
                           survival_calc = TRUE,
                           survival_column = "SR36",
                           survival_values = c(0, 1),
                           GEI = TRUE,
                           Trial_column = "TRIAL",
                           GEI_term = "TRIAL:PROG") {
    if (is.null(factor_for_adj_mean)) {
      adjusted_means <- NULL
    } else {
      require(emmeans)
      # Obtain the adjusted means for the 'mean' factor
      adjusted_means <- emmeans(mod_lmer, as.formula(paste("~", factor_for_adj_mean))) |>
        as.data.frame()
    }
    
    if (!GEI) {
      # Calculate the phenotypic mean
      pheno_mean = mean(data[[trait]], na.rm = T)
      # Extract genetic parameter
      varcomp_lmer <- summary(mod_lmer)$varcor |>
        data.frame() |>
        column_to_rownames(var = "grp")
      
      n_obs <- data |>
        group_by(get(genetic_factor)) |>
        summarise(n = sum(!is.na(get(trait)))) |>
        summarise(mean_n = mean(n)) |>
        pull(mean_n)
      # Extract variance components
      VarComp <- as.data.frame(VarCorr(mod_lmer)) %>%
        dplyr::select(!c("var1", "var2")) %>%
        `colnames<-`(c("FV", "Variance", "Standard Deviation")) %>%
        `rownames<-`(., .$FV)
      
      Plot <- ifelse(is.null(Plot_column), NA, Plot_column)
      nRep <- length(unique(data[[Rep_column]]))
      
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
      blups <- blups_with_se |>
        mutate(
          # Calculate the new mean BV + phenotypic mean
          BLUP_plus_mean = solution + pheno_mean,
          # calculate PEV
          PEV = SE ^ 2,
          # Calculate accuracy
          accuracy = sqrt(1 - PEV / GenPar[which(GenPar[["Components"]] == "Vg"), "Estimates"])
        ) |>
        rename(BV = solution)
      
      if (survival_calc) {
        # get the survival of clones
        survival <- data |>
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
        blups <- blups |>
          left_join(survival, by = genetic_factor)
        
      }
      # MET model
    } else {
      # Calculate the phenotypic mean
      pheno_mean = data |> 
        group_by(get(Trial_column)) |> 
        summarise(
          Mean = mean(get(trait), na.rm = T)
        ) |> 
        `colnames<-`(
          c(
            Trial_column,
            "Mean"
          )
        )
      # Extract genetic parameter
      varcomp_lmer <- summary(mod_lmer)$varcor |>
        data.frame() |>
        column_to_rownames(var = "grp")
      
      n_obs <- data |>
        group_by(get(Trial_column), get(genetic_factor)) |>
        summarise(n = sum(!is.na(get(trait)))) |>
        summarise(mean_n = mean(n)) |>
        column_to_rownames(var = "get(Trial_column)")
      # Extract variance components
      VarComp <- as.data.frame(VarCorr(mod_lmer)) %>%
        dplyr::select(!c("var1", "var2")) %>%
        `colnames<-`(c("FV", "Variance", "Standard Deviation")) %>%
        `rownames<-`(., .$FV)
      
      Plot <- ifelse(is.null(Plot_column), NA, Plot_column)
      nRep <- data |> 
        group_by(get(Trial_column)) |> 
        summarise(
          nRep = length(unique(get(Rep_column)))
        )
      
      GenPar <- tibble(
        Vg = VarComp[genetic_factor, "Variance"],
        Vplot = ifelse(!is.null(Plot_column), VarComp[Plot, "Variance"], Plot),
        Vrandom = ifelse(is.null(random), NA, VarComp[random, "Variance"]),
        Vres = VarComp["Residual", "Variance"],
        Vpheno = sum(Vg, Vplot, Vrandom, Vres, na.rm = T),
        GEI = VarComp[GEI_term, "Variance"] / Vpheno,
        `h2g` = Vg / Vpheno,
        `h2m` = Vg / sum(Vg, (Vplot / mean(nRep$nRep)), (Vres / mean(n_obs$mean_n)), na.rm = T),
        `c2p` = ifelse(is.na(Vplot), NA, Vplot / Vpheno)
      ) %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`("Estimates") %>%
        rownames_to_column(var = "Components")
      
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
      blups <- blups_with_se |>
        mutate(
          # Calculate the new mean BV + phenotypic mean
          BLUP_plus_mean = solution + mean(pheno_mean$Mean),
          # calculate PEV
          PEV = SE ^ 2,
          # Calculate accuracy
          accuracy = sqrt(1 - PEV / GenPar[which(GenPar[["Components"]] == "Vg"), "Estimates"])
        ) |>
        rename(BV = solution)
      
      # Extract gei effects
      
      gei <- ranef(mod_lmer)[[GEI_term]] |> 
        data.frame() |> 
        `colnames<-`("blup_gei") |> 
        rownames_to_column(genetic_factor) |> 
        separate(
          col = genetic_factor,
          into = c(Trial_column, genetic_factor),
          sep = ":"
        ) |> 
        mutate(
          blup_gei = ifelse(blup_gei == 0, NA, blup_gei)
        ) |> 
        left_join(
          blups |> 
            dplyr::select(
              genetic_factor,
              BV
            ),
          by = genetic_factor
        ) |> 
        mutate(blup_gei = blup_gei + BV) |> 
        dplyr::select(!BV) |> 
        pivot_wider(
          names_from = Trial_column,
          values_from = "blup_gei", 
          names_prefix = "BLUP_"
        ) 
      
      
      ## Merge blups and blups_gei
      blups <- blups |> 
        left_join(
          gei,
          by = genetic_factor
        )
      
      if (survival_calc) {
        # get the survival of clones
        survival <- data |>
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
        blups <- blups |>
          left_join(survival, by = genetic_factor)
        
      }
    }
    
    return(list(
      GenPar = GenPar,
      blups = blups,
      adjusted_means = adjusted_means
    ))
  }