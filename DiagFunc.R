#---------------------------- Function for diagnostic analysis -------------------------------#
DiagFunc <-
  function(Rep = "REP",
           Trait = NULL,
           Trat = NULL,
           data1 = data,
           track_feature = "",
           Exp = "Test",
           plot_diag1 = TRUE,
           plot_diag2 = TRUE,
           plotBox1 = TRUE,
           lambda.boxcox = c(0.9, 1.1),
           # Número de análises realizadas simultaneamente
           nDiag = 3,
           verbose = FALSE,
           # Column names to return
           ColumnNames_To_Return =
             c(NULL, NULL)) {
    # Input Exp column if necessary
    if (is.null(Exp)) {
      data1$Exp <- "Exp1"
      Exp <- "Exp"
      Exp2 <- NULL
    } else {
      Exp2 <- "Exp2"
    }
    
    # Set up track_feature
    if (track_feature != "") {
      track_feature_label <- paste0(": ", unique(data1[[track_feature]]))
    } else {
      track_feature_label = track_feature
    }
    
    
    if (!require("pacman")) {
      install.packages("pacman")
    }
    
    pacman::p_load(MASS, lmerTest, nortest, car, moments, Hmisc, tidyverse)
    
    ColumnNames_To_Return <- c(ColumnNames_To_Return, Trait)
    
    if (is.null(Rep)) {
      data1$REP <- 1
      Rep <- "REP"
    }
    colnames(data1)[match(c(Rep, Trat, Exp, Trait), colnames(data1))] <-
      c("Rep", "Trat", "Exp", "Trait")
    
    
    if (length(colnames(data1)[duplicated(colnames(data1))]) != 0) {
      stop(paste("Please, check and rename the following columns in the dataset:",
                 colnames(data1)[duplicated(colnames(data1))]))
    }
    
    # Test and boxcox transformation, if necessary
    FIT <- list()
    BOX <- list()
    LAMBDA <- list()
    DATA <- list()
    MODEL <- list()
    for (i in unique(unlist(data1[["Exp"]]))) {
      if (!is.null(Rep) & (
        data1 %>%
        dplyr::filter(Exp == "exp1") |>
        dplyr::select(Rep, Trait) |>
        filter(!is.na(Trait) & !duplicated(Rep)) |>
        dplyr::select(Trait) |>
        unlist() %>%
        sum(.) %>%
        is.na(.)
      ) == T) {
        Model <-
          as.formula(Trait ~ as.factor(Trat) + as.factor(Rep),
                     env = parent.frame(0))
      } else {
        Model <- as.formula(Trait ~ as.factor(Trat), env = parent.frame(0))
      }
      
      Model2 <-
        as.formula(Trait ~ as.factor(Trat), env = parent.frame(0))
      fit <-
        lm(Model,
           data = subset(data1, Exp == i))
      
      FIT[[i]] <- fit
      MODEL[[i]] <- Model
      DATA[[i]] <- subset(data1, Exp == i)
    }
    for (i in unique(unlist(data1[["Exp"]]))) {
      # boxcox
      box <-
        boxcox(
          MODEL[[i]],
          lambda = seq(-3, 3, 0.05),
          plotit = F,
          data = DATA[[i]]
        )
      lambda <- box$x[which(box$y == max(box$y))]
      
      BOX[[i]] <- box
      LAMBDA[[i]] <- lambda
      
    }
    
    
    for (i in unique(unlist(data1[["Exp"]]))) {
      DATA[[i]]$TraitT <- NA
      if (BOX[[i]]$x[which(BOX[[i]]$y == max(BOX[[i]]$y))] < lambda.boxcox[1] |
          BOX[[i]]$x[which(BOX[[i]]$y == max(BOX[[i]]$y))] > lambda.boxcox[2]) {
        cat(paste0(
          "\nTransformation to ",
          i, " - ", Trait, " ", track_feature, track_feature_label,
          ": ",
          "lambda = ",
          round(LAMBDA[[i]], 3),
          "\n"
        ))
        # Transformando caractere (Trait)
        
        DATA[[i]]$TraitT <- DATA[[i]][["Trait"]] ^ LAMBDA[[i]]
        # data1 <- bind_rows(DATA)
        
        T_BoxCox <- TRUE
      } else {
        cat(paste0("\nNo transformation is required for ", i," - ", Trait, " ", track_feature, 
                   track_feature_label, "\n"))
        DATA[[i]]$TraitT <- DATA[[i]][["Trait"]]
        #data1 <- bind_rows(DATA)
        T_BoxCox <- FALSE
      }
    }
    data1 <- bind_rows(DATA)
    rm(DATA)
    
    # if (plotBox1 == TRUE) {
    #   for (i in unique(unlist(data1[["Exp"]]))) {
    #     boxcox(
    #       FIT[[i]],
    #       seq(-3, 3, 0.1),
    #
    #       xlab = bquote( ~ lambda == .(LAMBDA[[i]]))
    #     )
    #     title(i)
    #   }
    # }
    data1 <- data1
    
    colnames(data1)[match(c("Rep", "Trat", "Exp", "Trait"), colnames(data1))] <-
      c(Rep, Trat, Exp, Trait)
    
    if (is.null(Exp2)) {
      data1 <- data1 |>
        rename(Experiment = Exp)
      
      Exp <- "Experiment"
    }
    
    if (!is.null(Rep) & (
      data1 %>%
      filter(get(Exp) == i) |>
      dplyr::select(Rep, Trait) |>
      filter(!is.na(Trait) & !duplicated(get(Rep))) |>
      dplyr::select(Trait) |>
      unlist() %>%
      sum(.) %>%
      is.na(.)
    ) == F & length(unique(data1[[Rep]])) > 1) {
      Model2 <-
        as.formula(TraitT ~ as.factor(get(Trat)) + as.factor(get(Rep)))
    } else {
      Model2 <- as.formula(TraitT ~ as.factor(get(Trat)))
    }
    
    MOD <- list()
    for (i in unique(unlist(data1[[Exp]]))) {
      mod1 <-
        lm(Model2,
           data = data1 %>%
             filter(get(Exp) == i),
           na.action = na.omit)
      
      # mod2 = sem NA nos resíduos
      mod2 <-
        lm(Model2,
           data = data1 %>%
             filter(get(Exp) == i),
           na.action = na.exclude)
      MOD[[i]] <- list(
        mod1 = mod1,
        mod2 = mod2,
        DATA = data1 %>%
          filter(get(Exp) == i)
      )
    }
    
    
    for (i in unique(unlist(data1[[Exp]]))) {
      # Extraindo o resíduo estudentizado
      MOD[[i]]$DATA[, paste0("residuals.", Trait)] <-
        rstudent(MOD[[i]]$mod2, type = "pearson")
      # Extraindo os valores preditos
      MOD[[i]]$DATA[, paste0("predict.", Trait)] <-
        predict(MOD[[i]]$mod2)
      # Obtendo os resíduos estudentizados
      # Sem NA's
      MOD[[i]]$rs <- rstudent(MOD[[i]]$mod1, type = "pearson")
      
      # Com NA's
      MOD[[i]]$rs.n.na <- rstudent(MOD[[i]]$mod2, type = "pearson")
    }
    
    if (plot_diag1 == TRUE) {
      YP <- list()
      LMAX_RS <- list()
      LMIN_RS <- list()
      LMAX_YP <- list()
      LMIN_YP <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
        if (max(MOD[[i]]$rs, na.rm = T) > 3 |
            min(MOD[[i]]$rs, na.rm = T) < -3) {
          DI1 <- "\nBefore removal of discrepant data"
        } else {
          DI1 <- ""
        }
      }
      # Escrevendo subtítulo nos gráficos
      SUBt <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
        if (T_BoxCox == T) {
          SUBt[[i]] <-
            bquote(`BoxCox Transformation:` ~ lambda == .(LAMBDA[[i]]))
        } else {
          SUBt[[i]] <- ""
        }
      }
      for (i in unique(unlist(data1[[Exp]]))) {
        # Histograma
        hist(
          MOD[[i]]$rs.n.na,
          main = "",
          xlab = "Studentized Residuals",
          ylab = paste("Probability density"),
          freq = F
        )
        title(
          main = paste0("Histogram for ", Trait, " trait: ", i, " ", track_feature, 
                        track_feature_label, " ", DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        xm <-
          seq(min(MOD[[i]]$rs, na.rm = T), max(MOD[[i]]$rs, na.rm = T), length = 40)
        ym <- dnorm(xm)
        lines(xm, ym)
        
        # qqNorm
        ## Plotando os pontos
        qqnorm(MOD[[i]]$rs.n.na,
               main = "")
        # Plotando a reta
        qqline(MOD[[i]]$rs.n.na,
               main = "",
               xlab = "",
               ylab = "")
        title(
          main = paste("Normality Graph - qqNorm:", Trait, "-", i, track_feature,
                       track_feature_label, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        # qqPlot
        qqPlot(MOD[[i]]$rs.n.na,
               main = "",
               ylab = "Studentized Residuals")
        title(
          main = paste("Normality Graph - qqPlot:", Trait, "-", i, track_feature,
                       track_feature_label, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        # Checando presença de valores discrepantes
        # Armazenando yp = Valores preditos
        YP[[i]] <- predict(MOD[[i]]$mod1, na.rm = T)
        
        # Marcando limites para valores normais
        LMAX_RS[[i]] <- max(MOD[[i]]$rs, 3.0, na.rm = T) + 0.1
        LMIN_RS[[i]] <- min(MOD[[i]]$rs, -3.0, na.rm = T) - 0.1
        
        # Armazenando os valores preditos
        LMAX_YP[[i]] <- max(YP[[i]], na.rm = T) + 0.1
        LMIN_YP[[i]] <- min(YP[[i]], na.rm = T) - 0.1
        
        # Plotando o gráfico
        ## Configurando os eixos e títulos
        plot(
          c(LMIN_YP[[i]], LMAX_YP[[i]]),
          c(LMIN_RS[[i]], LMAX_RS[[i]]),
          "n",
          main = "",
          xlab = expression(paste(hat(y)[predito])),
          ylab = "Studentized Residuals"
        )
        title(
          main = paste("Residual Analysis:", Trait, "-", i, track_feature,
                       track_feature_label, DI1),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        ## Plotando as linhas
        abline(h = 0, col = "black")
        abline(h = 3.0, col = "red")
        abline(h = -3.0, col = "red")
        
        ## Plotando os pontos
        points(YP[[i]], MOD[[i]]$rs)
      }
    } else {
      DI1 <- ""
      SUBt <- list()
      for (i in unique(unlist(data1[[Exp]]))) {
        SUBt[[i]] <- ""
        message("plot_diag1 = FALSE: omitting graphs before discrepant analysis")
      }
    }
    
    NORMAL1 <- list()
    ASSIMETRIA1 <- list()
    CURTOSE1 <- list()
    HOMOCED1 <- list()
    # Normalidade e homocedasticidade
    for (i in unique(unlist(data1[[Exp]]))) {
      # Teste de normalidade
      Normal <- try({list(ad.test(MOD[[i]]$rs), lillie.test(MOD[[i]]$rs))})
      NORMAL1[[i]] <- Normal
      # Assimetria
      Assimetria <- skewness(MOD[[i]]$rs.n.na, na.rm = T)
      ASSIMETRIA1[[i]] <- Assimetria
      # Curtose
      Curtose <- kurtosis(MOD[[i]]$rs.n.na, na.rm = T)
      CURTOSE1[[i]] <- Curtose
      
      # Homocedasticidade
      homoced <- try({leveneTest(TraitT ~ as.factor(get(Trat)),
                            data = data1)})
      HOMOCED1[[i]] <- homoced
      
      if (T_BoxCox == TRUE) {
        label(NORMAL1[[i]]) <-
          paste("\nNormality Test for",
                i, track_feature, track_feature_label,
                "after BoxCox transformation",
                DI1)
        label(ASSIMETRIA1[[i]]) <-
          paste("\nSkewness in the data of ",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI1)
        label(CURTOSE1[[i]]) <-
          paste("\nKurtosis in the data of ",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI1)
        label(HOMOCED1[[i]]) <-
          paste("\nHomoscedasticity Test for",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI1)
      } else {
        label(NORMAL1[[i]]) <- paste("\nHomoscedasticity Test for", i, track_feature,
                                     track_feature_label, DI1)
        label(ASSIMETRIA1[[i]]) <-
          paste("\nSkewness in the data of", i, track_feature, track_feature_label, DI1)
        label(CURTOSE1[[i]]) <-
          paste("\nKurtosis in the data of", i, track_feature,
                track_feature_label, DI1)
        label(HOMOCED1[[i]]) <-
          paste("\nHomoscedasticity Test for", i, track_feature,
                track_feature_label, DI1)
      }
    }
    
    # Retirada de dados Discrepantes
    
    disc <- 0
    # rm(disc)
    # disc3 <- 0
    data2 <- data1
    # MOD2 <- list()
    
    repeat {
      data2 <- data2
      # MOD2 <- MOD
      
      if (!is.null(Rep) & (
        data1 %>%
        filter(get(Exp) == i) |>
        dplyr::select(Rep, Trait) |>
        filter(!is.na(Trait) & !duplicated(get(Rep))) |>
        dplyr::select(Trait) |>
        unlist() %>%
        sum(.) %>%
        is.na(.)
      ) == F & length(unique(data1[[Rep]])) > 1) {
        Model2 <-
          as.formula(TraitT ~ as.factor(get(Trat)) + as.factor(get(Rep)))
      } else {
        Model2 <- as.formula(TraitT ~ as.factor(get(Trat)))
      }
      DISCREPANT <- list()
      
      MOD2 <- list()
      for (i in unique(unlist(data2[[Exp]]))) {
        # mod1 = Com NA nos resíduos
        mod1 <-
          lm(Model2,
             data = data2 %>%
               filter(get(Exp) == i),
             na.action = na.omit)
        
        # mod2 = sem NA nos resíduos
        mod2 <-
          lm(Model2,
             data = data2 %>%
               filter(get(Exp) == i),
             na.action = na.exclude)
        MOD2[[i]] <- list(
          mod1 = mod1,
          mod2 = mod2,
          DATA = data2 %>%
            filter(get(Exp) == i)
        )
      }
      
      for (i in unique(unlist(data2[[Exp]]))) {
        # Extraindo o resíduo estudentizado
        MOD2[[i]]$DATA[, paste0("residuals.", Trait)] <-
          rstudent(MOD2[[i]]$mod2, type = "pearson")
        # Extraindo os valores preditos
        MOD2[[i]]$DATA[, paste0("predict.", Trait)] <-
          predict(MOD2[[i]]$mod2)
        # Obtendo os resíduos estudentizados
        # Sem NA's
        MOD2[[i]]$rs <- rstudent(MOD2[[i]]$mod1, type = "pearson")
        
        # Com NA's
        MOD2[[i]]$rs.n.na <-
          rstudent(MOD2[[i]]$mod2, type = "pearson")
      }
      for (i in unique(unlist(data2[[Exp]]))) {
        discrepant <-
          MOD2[[i]]$DATA[which(MOD2[[i]]$DATA[, paste0("residuals.", Trait)] < -3 |
                                 MOD2[[i]]$DATA[, paste0("residuals.", Trait)] > 3), colnames(MOD2[[i]]$DATA)]
        
        DISCREPANT[[i]] <- discrepant
        MOD2[[i]]$DATA[rownames(DISCREPANT[[i]]), "TraitT"] <- NA
        # MOD[[i]]$DATA[rownames(DISCREPANT[[i]]), Trait] <- NA
      }
      
      DATA2 <- list()
      for (i in unique(unlist(data2[[Exp]]))) {
        DATA2[[i]] <- MOD2[[i]]$DATA
      }
      for (i in unique(unlist(data2[[Exp]]))) {
        DATA2[[i]][[Trait]] <- ifelse(is.na(DATA2[[i]]$TraitT), NA,
                                      DATA2[[i]][[Trait]])
        
      }
      data2 <- bind_rows(DATA2)
      
      disc <- disc + 1
      
      if (disc > nDiag) {
        break
      }
    }
    
    
    if (plot_diag2 == TRUE) {
      YP2 <- list()
      LMAX_RS2 <- list()
      LMIN_RS2 <- list()
      LMAX_YP2 <- list()
      LMIN_YP2 <- list()
      
      if (max(MOD[[i]]$rs, na.rm = T) > 3 |
          min(MOD[[i]]$rs, na.rm = T) < -3) {
        DI2 <- "\nAfter removal of discrepant data"
      } else {
        DI2 <- ""
        
      }
      
      for (i in unique(unlist(data2[[Exp]]))) {
        # Histograma
        hist(
          MOD2[[i]]$rs.n.na,
          main = "",
          xlab = "Studentized Residuals",
          ylab = paste("Probability density"),
          freq = F
        )
        title(
          main = paste("Histogram for", Trait, " trait:", i, track_feature,
                       track_feature_label, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        xm <-
          seq(min(MOD[[i]]$rs, na.rm = T), max(MOD[[i]]$rs, na.rm = T), length = 40)
        ym <- dnorm(xm)
        lines(xm, ym)
        
        # qqNorm
        ## Plotando os pontos
        qqnorm(MOD2[[i]]$rs.n.na,
               main = "")
        # Plotando a reta
        qqline(MOD2[[i]]$rs.n.na,
               main = "",
               xlab = "",
               ylab = "")
        title(
          main = paste("Normality Graph - qqNorm:", Trait, "-", i, track_feature,
                       track_feature_label, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        # qqPlot
        qqPlot(MOD2[[i]]$rs.n.na,
               main = "",
               ylab = "Studentized Residuals")
        title(
          main = paste("Normality Graph - qqPlot:", Trait, "-", i, track_feature,
                       track_feature_label, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        # Checando presença de valores discrepantes
        # Armazenando yp = Valores preditos
        YP2[[i]] <- predict(MOD2[[i]]$mod1, na.rm = T)
        
        # Marcando limites para valores normais
        LMAX_RS2[[i]] <- max(MOD2[[i]]$rs, 3.0, na.rm = T) + 0.1
        LMIN_RS2[[i]] <- min(MOD2[[i]]$rs, -3.0, na.rm = T) - 0.1
        
        # Armazenando os valores preditos
        LMAX_YP2[[i]] <- max(YP2[[i]], na.rm = T) + 0.1
        LMIN_YP2[[i]] <- min(YP2[[i]], na.rm = T) - 0.1
        
        # Plotando o gráfico
        ## Configurando os eixos e títulos
        plot(
          c(LMIN_YP2[[i]], LMAX_YP2[[i]]),
          c(LMIN_RS2[[i]], LMAX_RS2[[i]]),
          "n",
          main = "",
          xlab = expression(paste(hat(y)[predito])),
          ylab = "Studentized Residuals"
        )
        title(
          main = paste("Residual Analysis:", Trait, "-", i, track_feature,
                       track_feature_label, DI2),
          cex.main = 1,
          sub = SUBt[[i]]
        )
        
        ## Plotando as linhas
        abline(h = 0, col = "black")
        abline(h = 3.0, col = "red")
        abline(h = -3.0, col = "red")
        
        ## Plotando os pontos
        points(YP2[[i]], MOD2[[i]]$rs)
      }
    } else {
      DI2 <- ""
      message("plot_diag2 = FALSE: omitting graphs after discrepant analysis")
    }
    
    NORMAL2 <- list()
    ASSIMETRIA2 <- list()
    CURTOSE2 <- list()
    HOMOCED2 <- list()
    # Normalidade e homocedasticidade
    for (i in unique(unlist(data2[[Exp]]))) {
      # Teste de normalidade
      Normal <-
        try({list(ad.test(MOD2[[i]]$rs), lillie.test(MOD2[[i]]$rs))})
      NORMAL2[[i]] <- Normal
      # Assimetria
      Assimetria <- skewness(MOD2[[i]]$rs.n.na, na.rm = T)
      ASSIMETRIA2[[i]] <- Assimetria
      # Curtose
      Curtose <- kurtosis(MOD2[[i]]$rs.n.na, na.rm = T)
      CURTOSE2[[i]] <- Curtose
      
      # Homocedasticidade
      homoced <- try({leveneTest(TraitT ~ as.factor(get(Trat)),
                            data = data2)})
      HOMOCED2[[i]] <- homoced
      
      if (T_BoxCox == TRUE) {
        label(NORMAL2[[i]]) <-
          paste("\nNormality test for ",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI2)
        label(ASSIMETRIA2[[i]]) <-
          paste("\nSkewness in the data of ",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI2)
        label(CURTOSE2[[i]]) <-
          paste("\nKurtosis in the data of ",
                i, track_feature, track_feature_label,
                " após transformação BoxCox",
                DI2)
        label(HOMOCED2[[i]]) <-
          paste("\nHomoscedasticity Test for",
                i, track_feature, track_feature_label,
                " after BoxCox transformation",
                DI2)
      } else {
        label(NORMAL2[[i]]) <-
          paste0("\nNormality test for ", i, track_feature,
                 track_feature_label, DI2)
        label(ASSIMETRIA2[[i]]) <-
          paste0("\nSkewness in the data of ", i, track_feature,
                 track_feature_label, DI2)
        label(CURTOSE2[[i]]) <-
          paste0("\nKurtosis in the data of ", i, track_feature,
                 track_feature_label, DI2)
        label(HOMOCED2[[i]]) <-
          paste0("\nHomocedasticity Test for ", i, track_feature,
                 track_feature_label, DI2)
      }
    }
    
    Discrepantes <- suppressMessages(anti_join(data1, data2))
    
    
    if (verbose == TRUE) {
      print(
        list(
          Outliers = Discrepantes,
          `Normality test before discrepant analylis` = NORMAL1,
          `Normality test after discrepant analysis` = NORMAL2,
          `Skewness before discrepant analysis` = ASSIMETRIA1,
          `Skewnwss after discrepant analysis` = ASSIMETRIA2,
          `Kurtosis before discrepant analysis` = CURTOSE1,
          `Kurtosis after discrepant analysis` = CURTOSE2,
          `Homocedasticity test before discrepant analysis` = HOMOCED1,
          `Homocedasticity test after discrepant analysis` = HOMOCED2
        )
      )
    }
    
    Output <-
      list(
        DataOriginals = data1,
        DataAfterDiag = data2[, c(ColumnNames_To_Return)],
        Discrepants = Discrepantes,
        `Normality and Homocedasticity` = list(
          `Normality before diagnostic analysis` = NORMAL1,
          `Normality after diagnostic analysis` = NORMAL2,
          `Skewness before diagnostic analysis` = ASSIMETRIA1,
          `Skewness after diagnostic analysis` = ASSIMETRIA2,
          `Kurtosis before diagnostic analysis` = CURTOSE1,
          `Kurtosis after diagnostic analysis` = CURTOSE2,
          `Homocedasticity before diagnostic analysis` = HOMOCED1,
          `Homocedasticity after diagnostic analysis` = HOMOCED2
        ),
        Lambda = LAMBDA
      )
    return(Output)
  }
