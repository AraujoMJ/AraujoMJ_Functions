#--------------------------- Function for extract heritability and correlations ------------------#

Extract_h2a_breedR <- function(Model = NULL,
                               random_effect = NULL,
                               model_type = "std") {
  Var = summary(Model)$var
  # For ar1Comp Model
  if (model_type == "ar1Comp") {
    test_model <- names(Var)
    if (is.null(test_model)) {
      stop("model_type is not 'ar1Comp' model, please, set properly!")
    }
    if (!"spatial" %in% names(Var)) {
      stop("model_type is not 'ar1Comp' model, please, set properly!")
    }
    if (!is.null(random_effect)) {
      Names_Var <- c(names(Var), "competition", "corr_gen_comp")
      Variances_Comp <- data.frame(
        Stats = Names_Var,
        Estimates = c(
          Var[[random_effect]],
          Var$genetic[1, 1],
          Var$pec,
          Var$spatial,
          Var$Residual,
          Var$genetic[2, 2],
          Var$genetic[1, 2]
        ),
        S.E = as.numeric(NA)
      ) |>
        `rownames<-`(Names_Var)
      # Correlation between direct and indirect genetic effect
      Corr_gen_comp <-
        Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
      Variances_Comp["corr_gen_comp", "Estimates"] <- Corr_gen_comp
      # spatial variance
      Spatial <- Variances_Comp["spatial", "Estimates"]
      
      # c2_random_effect
      C2 <- list()
      for (i in random_effect) {
        C2[[i]] <- data.frame(
          Stats = paste0("c2_", i),
          Estimates = Variances_Comp[random_effect, "Estimates"] / (sum(abs(
            Variances_Comp$Estimates
          )) - abs(Corr_gen_comp) - Spatial),
          S.E = NA
        )
        
      }
      c2_random_effect <- bind_rows(C2)
      
      # Variance components
      # h2a
      h2a <-
        c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(abs(
          Variances_Comp$Estimates
        )) - abs(Corr_gen_comp) - Spatial), NA)
      Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
        `rownames<-`(NULL)
      
    } else {
      Names_Var <- c(names(Var), "competition", "corr_gen_comp")
      Variances_Comp <- data.frame(
        Stats = Names_Var,
        Estimates = c(
          Var$genetic[1, 1],
          Var$pec,
          Var$spatial,
          Var$Residual,
          Var$genetic[2, 2],
          Var$genetic[1, 2]
        ),
        S.E = as.numeric(NA)
      ) |>
        `rownames<-`(Names_Var)
      # Correlation between direct and indirect genetic effect
      Corr_gen_comp <-
        Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
      Variances_Comp["corr_gen_comp", "Estimates"] <- Corr_gen_comp
      # spatial variance
      Spatial <- Variances_Comp["spatial", "Estimates"]
      
    }
    # Variance components
    # h2a
    h2a <-
      c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(abs(
        Variances_Comp$Estimates
      )) - abs(Corr_gen_comp) - Spatial), NA)
    
    Data2 <- rbind(Variances_Comp, h2a) |>
      `rownames<-`(NULL)
  } else {
    if (model_type == "Comp") {
      if (is.null(names(Var))) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if ("spatial" %in% names(Var)) {
        stop("model_type is not 'Comp' model, please, set properly!")
      }
      if (!is.null(random_effect)) {
        Names_Var <- c(names(Var), "competition", "corr_gen_comp")
        Variances_Comp <- data.frame(
          Stats = Names_Var,
          Estimates = c(
            Var[[random_effect]],
            Var$genetic[1, 1],
            Var$pec,
            Var$Residual,
            Var$genetic[2, 2],
            Var$genetic[1, 2]
          ),
          S.E = as.numeric(NA)
        ) |>
          `rownames<-`(Names_Var)
        # Correlation between direct and indirect genetic effect
        Corr_gen_comp <-
          Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
        Variances_Comp["corr_gen_comp", "Estimates"] <-
          Corr_gen_comp
        # c2_random_effect
        C2 <- list()
        for (i in random_effect) {
          C2[[i]] <- data.frame(
            Stats = paste0("c2_", i),
            Estimates = Variances_Comp[random_effect, "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp),
            S.E = NA
          )
          
        }
        c2_random_effect <- bind_rows(C2)
        # h2a
        h2a <-
          c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp), NA)
        
        Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
          `rownames<-`(NULL)
        
      } else{
        Names_Var <- c(names(Var), "competition", "corr_gen_comp")
        Variances_Comp <- data.frame(
          Stats = Names_Var,
          Estimates = c(
            Var$genetic[1, 1],
            Var$pec,
            Var$Residual,
            Var$genetic[2, 2],
            Var$genetic[1, 2]
          ),
          S.E = as.numeric(NA)
        ) |>
          `rownames<-`(Names_Var)
        # Correlation between direct and indirect genetic effect
        Corr_gen_comp <-
          Variances_Comp["corr_gen_comp", "Estimates"] / (sqrt(Variances_Comp["genetic", "Estimates"]) * sqrt(Variances_Comp["competition", "Estimates"]))
        Variances_Comp["corr_gen_comp", "Estimates"] <-
          Corr_gen_comp
        # c2_random_effect
        C2 <- list()
        for (i in random_effect) {
          C2[[i]] <- data.frame(
            Stats = paste0("c2_", i),
            Estimates = Variances_Comp[i, "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp),
            S.E = NA
          )
          
        }
        c2_random_effect <- bind_rows(C2)
        
        # h2a
        h2a <-
          c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(Variances_Comp$Estimates) - Corr_gen_comp), NA)
        
        Data2 <- rbind(Variances_Comp, h2a) |>
          `rownames<-`(NULL)
      }
    } else {
      # AR1 Model
      if (model_type == "ar1") {
        if (is.null(dimnames(Var))) {
          stop("model_type is not 'ar1' model, please, set properly!")
        }
        if (!"spatial" %in% dimnames(Var)[[1]]) {
          stop("model_type is not 'ar1' model, please, set properly!")
        }
        if (!is.null(random_effect)) {
          Names_Var <- c(dimnames(Var)[[1]], "rho_r", "rho_c")
          RHO <-
            summary(Model)$rho[which.max(summary(Model)$rho$loglik), ]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var[random_effect, "Estimated variances"],
                          Var["genetic", "Estimated variances"],
                          Var["spatial", "Estimated variances"],
                          Var["Residual", "Estimated variances"],
                          RHO[["rho_r"]],
                          RHO[["rho_c"]]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # spatial variance
          Spatial <- Variances_Comp["spatial", "Estimates"]
          
          # c2_random_effect
          C2 <- list()
          for (i in random_effect) {
            C2[[i]] <- data.frame(
              Stats = paste0("c2_", i),
              Estimates = Variances_Comp[i, "Estimates"] / (sum(
                abs(Variances_Comp$Estimates)
              ) - Spatial),
              S.E = NA
            )
            
          }
          c2_random_effect <- bind_rows(C2)
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            ) - Spatial), NA)
          Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
            `rownames<-`(NULL)
          
        } else {
          Names_Var <- c(dimnames(Var)[[1]], "rho_r", "rho_c")
          RHO <-
            summary(Model)$rho[which.max(summary(Model)$rho$loglik), ]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var["genetic", "Estimated variances"],
                          Var["spatial", "Estimated variances"],
                          Var["Residual", "Estimated variances"],
                          RHO[["rho_r"]],
                          RHO[["rho_c"]]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # spatial variance
          Spatial <- Variances_Comp["spatial", "Estimates"]
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            ) - Spatial), NA)
          Data2 <- rbind(Variances_Comp, h2a) |>
            `rownames<-`(NULL)
        }
      } else{
        if (is.null(dimnames(Var))) {
          stop("model_type is not 'std' model, please, set properly!")
        }
        if ("spatial" %in% dimnames(Var)[[1]]) {
          stop("model_type is not 'std' model, please, set properly!")
        }
        # Standard model
        if (!is.null(random_effect)) {
          Names_Var <- dimnames(Var)[[1]]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var[random_effect, "Estimated variances"],
                          Var["genetic", "Estimated variances"],
                          Var["Residual", "Estimated variances"]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # c2_random_effect
          C2 <- list()
          for (i in random_effect) {
            C2[[i]] <- data.frame(
              Stats = paste0("c2_", i),
              Estimates = Variances_Comp[i, "Estimates"] / (sum(
                abs(Variances_Comp$Estimates)
              )),
              S.E = NA
            )
            
          }
          c2_random_effect <- bind_rows(C2)
          
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            )), NA)
          Data2 <- rbind(Variances_Comp, c2_random_effect, h2a) |>
            `rownames<-`(NULL)
          
        } else {
          Names_Var <- dimnames(Var)[[1]]
          Variances_Comp <- data.frame(
            Stats = Names_Var,
            Estimates = c(Var["genetic", "Estimated variances"],
                          Var["Residual", "Estimated variances"]),
            S.E = as.numeric(NA)
          ) |>
            `rownames<-`(Names_Var)
          
          # Variance components
          # h2a
          h2a <-
            c("h2a", Variances_Comp["genetic", "Estimates"] / (sum(
              abs(Variances_Comp$Estimates)
            )), NA)
          Data2 <- rbind(Variances_Comp, h2a) |>
            `rownames<-`(NULL)
        }
      }
    }
  }
  return(Data2)
}
