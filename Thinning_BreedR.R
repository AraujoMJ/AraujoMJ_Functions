#' Thinning Strategies Based on Breeding Values
#'
#' This function performs thinning strategies for tree breeding trials using breeding values (BV).
#' It groups families into different categories (up to four groups) based on inflection points in the BV distribution.
#' It offers different strategies for selection between and within families (BW strategy) and selection of the best individuals (BI strategy).
#'
#' @param BV_Column Character. The column name in the data that contains the breeding values (BV). Default is "a_total".
#' @param Trait Character. Name of the trait being analyzed.
#' @param BV_fam Data frame. Contains breeding values of the families.
#' @param Data_Total Data frame. Contains breeding values and experimental information of the trees.
#' @param Family_Data_Total Character. Name of the column in `Data_Total` containing family identifiers.
#' @param Bloc_Column Character. The column in `Data_Total` for blocks or experimental plots. Default is "Block".
#' @param Family_Column Character. The column in `BV_fam` that specifies family identifiers. Default is "Family".
#' @param nGroups Integer. Number of groups for thinning strategies. Must be 2, 3, or 4. Default is 4.
#' @param label.group.y Numeric vector. Positions for group labels on the y-axis. Default is `c(1, 1, 1, 1)`.
#' @param Plot.Rank Logical. If TRUE, plot the ranking of families based on BV. Default is TRUE.
#' @param save_plot_rank Logical. If TRUE, save the plot of families' ranking. Default is TRUE.
#' @param IS Numeric vector. Selection intensities for the BI strategy. Default is NULL.
#' @param id Character. Column name for unique identifiers of trees. Default is "ID".
#' @param nGroups3 Character. Specifies whether the third group should be formed on the "left" or "right" side of the inflection point. Default is "left".
#' @param STP Logical. If TRUE, thinning strategy will be applied Single-Tree-Plot (STP). Default is FALSE.
#' @param n.dodge_plot1 Integer. This argument allows to use different n.dodge values when plotting the graph of BV ~ x-axis factor (e.g., families). Default is 1.
#' @param angle_plot1 Numeric. The angle of the x-axis text labels in the plot. Default is 45 degrees.
#' @param nt_alpha1 Numeric vector. Parameters for inflection point analysis on the best families. Default is `c(2, 5)`.
#' @param nt_alpha2 Numeric vector. Parameters for inflection point analysis on the worst families. Default is `c(2, 5)`.
#' @param finfl_Value Character vector. Specifies whether to use "max" or "min" for finding inflection points in both the best and worst families. Must be a vector of length 2. Default is `c("max", "max")`.
#' @param additional_layer_plot1 ggplot2 layer. Additional ggplot layers to be added to the plot of family rankings. Default is `theme(axis.text.x = element_text(size = 8))`.
#' @param seq_combinations Numeric vector. Custom sequence for generating combinations in thinning strategy. Default is NULL.
#' @param length_seq_combinations Numeric. Length of the sequence for combinations of selected plants per group. Default is 2.
#' @param save_table_xlsx Logical. If TRUE, save the output tables in Excel format. Default is TRUE.
#' @param export_id Character. A string to customize the output file names by appending it to the file name. Default is an empty string.
#' @param format_plot Character. Specifies the file format for the plot output. Default is ".tiff".
#'
#' @return A list containing:
#' \item{Thinning}{Data frame summarizing the thinning strategies and corresponding genetic gains and effective number (NE).}
#' \item{Strategies}{A list of data frames containing the selected families and trees under different thinning strategies.}
#' \item{BV_Fam}{Data frame summarizing breeding values for families, including their assigned group.}
#' \item{ggPlot_families_rank}{A ggplot object showing the ranking of families based on BV.}
#'
#' @details
#' The thinning strategies are based on dividing families into groups based on inflection points in the breeding values. These groups can be used in selection between families and within families (BW strategy) or to select the top-performing individuals regardless of family (BI strategy). The function allows for flexibility in the number of groups (2, 3, or 4) and customizes the selection process based on the input parameters.
#'
#' @note
#' The function relies on the `RootsExtremaInflections` and `gtools` packages for finding inflection points and generating combinations, respectively. The function will automatically install these packages if they are not already installed.
#'
#' @import dplyr ggplot2 openxlsx gtools RootsExtremaInflections
#' @examples
#' \dontrun{
#' BV <- data.frame(
#'   Family = 1:10,
#'   a_total = rnorm(10)
#' )
#' Data_Total <- data.frame(
#'   Family = rep(1:10, each = 10),
#'   Block = rep(1:10, 10),
#'   ID = 1:100,
#'   a = rnorm(100)
#' )
#' Trait <- "Height"
#' Thinning_BreedR(BV_fam = BV, Trait = Trait, Data_Total = Data_Total, Family_Data_Total = "Family")
#' }
#' @export


#---------------------- Function for thinning strategies --------------------------
Thinning_BreedR <- function(BV_Column = "a_total",
                            Trait = Trait,
                            BV_fam = BV_Family,
                            Data_Total = BV$Data_Total,
                            Family_Data_Total = Family,
                            Bloc_Column = "Block",
                            Family_Column = "Family",
                            nGroups = 4,
                            label.group.y = c(1, 1, 1, 1),
                            Plot.Rank = TRUE,
                            save_plot_rank = TRUE,
                            IS = NULL,
                            id = "ID",
                            nGroups3 = "left",
                            STP = FALSE,
                            n.dodge_plot1 = 1,
                            angle_plot1 = 45,
                            nt_alpha1 = c(2, 5),
                            nt_alpha2 = c(2, 5),
                            finfl_Value = c("max", "max"),
                            additional_layer_plot1 = theme(axis.text.x = element_text(size = 8)),
                            seq_combinations = NULL,
                            length_seq_combinations = 2,
                            save_table_xlsx = TRUE,
                            export_id = "",
                            format_plot = ".tiff",
                            verbose = FALSE) {
  if (verbose) {
    message("Starting Thinning_BreedR function...")
  }
  
  # Check if necessary packages are loaded
  if (!require("RootsExtremaInflections")) {
    install.packages("RootsExtremaInflections")
  }
  
  if (!require("gtools")) {
    install.packages("gtools")
  }
  
  # Parameter validation
  if (!nGroups %in% 2:4) {
    stop("nGroups should be: 2, 3 or 4")
  }
  
  if (any(!finfl_Value %in% c("max", "min"))) {
    stop("finfl_value should be a vector of length 2 of: 'max' or 'min'")
  }
  
  # Start of main function logic
  tryCatch({
    # Process BV_fam data frame
    if (verbose) {
      message("Processing BV_fam data frame...")
    }
    
    BV_fam <- BV_fam |>
      rename_with(~ "Family", all_of(Family_Column)) |> 
      arrange(desc(get(BV_Column))) |>
      mutate(Family = factor(Family, levels = Family))
    if (verbose) {
      message("Processed BV_fam data frame.")
    }
    
  }, error = function(e) {
    print(paste("Error processing BV_fam data frame:", e$message))
    stop(e)
  })
    
  # Merge Data_Total and BV_fam
  tryCatch({
    if (verbose) {
      message("Merging Data_Total and BV_fam...")
    }
    
    Data_Total <- Data_Total |> 
      mutate(
        Family = get(Family_Data_Total)
      ) |> 
      left_join(
        BV_fam |> 
          dplyr::select(
            Family,
            BV_Column
          ),
        by = "Family"
      )
    if (verbose) {
      message("Data merged successfully.")
    }
    
  }, error = function(e) {
    print(paste("Error merging Data_Total and BV_fam:", e$message))
    stop(e)
  })
    
  # Identify families with positive BV
  tryCatch({
    if (verbose) {
      message("Identifying families with positive BV...")
    }
    
    MajorBV_fam <- BV_fam |>
      filter(get(BV_Column) >= 0)
    Fam_Zero <- which.min(MajorBV_fam[[BV_Column]])
    MinorBV_fam <- BV_fam |>
      filter(get(BV_Column) < 0)
    
    if (verbose) {
      message("Families identified.")
    }
    
  }, error = function(e) {
    print(paste("Error identifying families with positive BV:", e$message))
    stop(e)
  })

  # Inflection points for best families
  tryCatch({
    if (verbose) {
      message("Finding inflection points for best families...")
    }
    
    Inflex1 <- RootsExtremaInflections::inflexi(
      x = 1:nrow(MajorBV_fam),
      y = MajorBV_fam[[BV_Column]],
      nt = nt_alpha1[1],
      alpha = nt_alpha1[2],
      i1 = 1,
      i2 = nrow(MajorBV_fam),
      plots = F
    )
    
    if (finfl_Value[1] == "max") {
      Inflexi_major <- max(Inflex1$finfl, na.rm = TRUE)
    } else {
      Inflexi_major <- min(Inflex1$finfl, na.rm = TRUE)
    }
    
    message(paste("Inflexi_major:", Inflexi_major))
  }, error = function(e) {
    print(paste("Error finding inflection points for best families:", e$message))
    stop(e)
  })
  
  # Inflection points for worst families
  tryCatch({
    if (verbose) {
      message("Finding inflection points for worst families...")
    }
   
    Inflex2 <- RootsExtremaInflections::inflexi(
      x = 1:nrow(MinorBV_fam),
      y = MinorBV_fam[[BV_Column]],
      nt = nt_alpha2[1],
      alpha = nt_alpha2[2],
      i1 = 1,
      i2 = nrow(MinorBV_fam),
      plots = F
    )
    
    if (finfl_Value[2] == "max") {
      Inflexi_minor <- Fam_Zero + max(Inflex2$finfl, na.rm = TRUE)
    } else {
      Inflexi_minor <- Fam_Zero + min(Inflex2$finfl, na.rm = TRUE)
    }
    message(paste("Inflexi_minor:", Inflexi_minor))
  }, error = function(e) {
    print(paste("Error finding inflection points for worst families:", e$message))
    stop(e)
  })
  
  # Setting up plot groups and annotations
  tryCatch({
    if (verbose) {
      message("Setting up plot groups and annotations...")  
    }
    
    Col1 = "gray70"
    Col2 = "blue"
    Col3 = "red"
    # line in Fam_Zero
    nGroup2 <-
      geom_vline(
        xintercept = Fam_Zero,
        linetype = 4,
        colour = Col1,
        size = 1
      )
    # nGroup = 3
    if (!nGroups3 %in% c("left", "right")) {
      stop("nGroups3 should be: left or right")
    }
    
    if (nGroups3 == "left") {
      nGroup3 <-
        list(
          geom_vline(
            xintercept = Fam_Zero,
            linetype = 4,
            colour = Col1,
            size = 1
          ),
          geom_vline(
            xintercept = Inflexi_major,
            linetype = 4,
            colour = Col2,
            size = 1
          )
        )
    } else {
      nGroup3 <-
        list(
          geom_vline(
            xintercept = Fam_Zero,
            linetype = 4,
            colour = Col1,
            size = 1
          ),
          geom_vline(
            xintercept = Inflexi_minor,
            linetype = 4,
            colour = Col2,
            size = 1
          )
        )
    }
    
    # nGroup = 4
    nGroup4 <-
      list(
        geom_vline(
          xintercept = Fam_Zero,
          linetype = 4,
          colour = Col1,
          size = 1
        ),
        geom_vline(
          xintercept = Inflexi_major,
          linetype = 4,
          colour = Col2,
          size = 1
        ),
        geom_vline(
          xintercept = Inflexi_minor,
          linetype = 4,
          colour = Col3,
          size = 1
        )
      )
    
    # if 2 groups:
    if (nGroups ==  2) {
      nGroupFinal <- nGroup2
      AnnotateGroup <- list(
        annotate(
          "text",
          x = Fam_Zero / 2,
          y = label.group.y[1],
          label = "G1"
        ),
        annotate(
          "text",
          x = (Fam_Zero + (nrow(BV_fam) - Fam_Zero) / 2),
          y = label.group.y[2],
          label = "G2"
        )
      )
    } else {
      # if 3 groups
      if (nGroups == 3) {
        nGroupFinal <- nGroup3
        AnnotateGroup <- list(
          annotate(
            "text",
            x = Inflexi_major / 2,
            y = label.group.y[1],
            label = "G1"
          ),
          annotate(
            "text",
            x = (Fam_Zero - (Fam_Zero - Inflexi_major) / 2),
            y = label.group.y[2],
            label = "G2"
          ),
          annotate(
            "text",
            x = Inflexi_minor + ( nrow(BV_fam) - Inflexi_minor) / 2,
            y = label.group.y[3],
            label = "G3"
          )
        )
      } else {
        # if 4 Groups
        nGroupFinal <- nGroup4
        AnnotateGroup <- list(
          annotate(
            "text",
            x = Inflexi_major / 2,
            y = label.group.y[1],
            label = "G1"
          ),
          annotate(
            "text",
            x = (Fam_Zero - (Fam_Zero - Inflexi_major) / 2),
            y = label.group.y[2],
            label = "G2"
          ),
          annotate(
            "text",
            x = Fam_Zero + (Inflexi_minor - Fam_Zero) / 2,
            y = label.group.y[3],
            label = "G3"
          ),
          annotate(
            "text",
            x = Inflexi_minor + (nrow(BV_fam) - Inflexi_minor) / 2,
            y = label.group.y[4],
            label = "G4"
          )
        )
      }
    }
    if (verbose) {
      message("Plot groups and annotations set.")
    }
    
  }, error = function(e) {
    print(paste("Error in plot group setup:", e$message))
    stop(e)
  })
  
  # Prepare the family ranking plot
  tryCatch({
    if (verbose) {
      message("Preparing the family ranking plot...")
    }
    
    families_rank <-
      ggplot(BV_fam, aes(x = Family, y = get(BV_Column))) +
      geom_point(size = 2) +
      ggtitle(label = paste("Thinning Strategies based on", Trait, "Trait")) +
      ylab(paste0("Additive genetic effect (", BV_Column, ") from ", Trait, " trait")) + xlab("Families") +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "gray60",
        size = 0.5
      ) +
      # Add the groups used
      nGroupFinal +
      # Annotate
      AnnotateGroup +
      # Change x-axis
      scale_x_discrete(
        guide = guide_axis(
          n.dodge = n.dodge_plot1,
          angle = angle_plot1
        )
      )
    if (Plot.Rank == TRUE) {
      if (verbose) {
        message("Displaying plot...")
      }
      
      plot(families_rank + additional_layer_plot1)
    }
    
    if (verbose) {
      message("Family ranking plot prepared successfully.")
    }
    
  }, error = function(e) {
    print(paste("Error in preparing the family ranking plot:", e$message))
    stop(e)
  })
  # Prepare Data_Total groups
  tryCatch({
    if (verbose) {
      message("Preparing Data_Total groups...")
    }
    
    # Map group for trees in dataset that have BV and experiment information
    Data_Total$Group <- NA
    
    # G1
    if (nGroups == 2) {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[1:Fam_Zero, "Family"]), "Group"] <-
        "G1"
    } else {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[1:Inflexi_major, "Family"]), "Group"] <-
        "G1"
    }
    
    # G2
    if (nGroups == 2) {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[(Fam_Zero + 1):nrow(BV_fam), "Family"]), "Group"] <-
        "G2"
    } else {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[(Inflexi_major + 1):Fam_Zero, "Family"]), "Group"] <-
        "G2"
    }
    
    # G3
    if (nGroups == 3) {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[(Fam_Zero + 1):nrow(BV_fam), "Family"]), "Group"] <-
        "G3"
    } else {
      if (nGroups == 4) {
        Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[(Fam_Zero + 1):Inflexi_minor, "Family"]), "Group"] <-
          "G3"
      }
      
    }
    
    # G4
    if (nGroups == 4) {
      Data_Total[Data_Total[, Family_Data_Total] %in% as.vector(BV_fam[(Inflexi_minor + 1):nrow(BV_fam), "Family"]), "Group"] <-
        "G4"
    }
    
    # Map Plot
    if (STP) {
      Data_Total[["Plot"]] <- Data_Total[[Bloc_Column]]
    } else {
      Data_Total[["Plot"]] <- paste0(
        Data_Total[[Family_Data_Total]], sep = ":", Data_Total[[Bloc_Column]]
      )
    }
    
    # Data_Total <- Data_Total |>
    #   mutate(Plot = ifelse(STP == TRUE, get(Bloc_Column), paste0(
    #     get(Family_Data_Total), sep = ":", get(Bloc_Column)
    #   )))
    
    # Combination for thinning strategies
    if (verbose) {
      message("Data_Total groups prepared successfully.")
    }
    
  }, error = function(e) {
    print(paste("Error in preparing Data_Total groups:", e$message))
    stop(e)
  })
  tryCatch({
    if (verbose) {
      message("Configuring thinning strategy combinations...")
    }
    
    if (STP == TRUE) {
      if (is.null(seq_combinations)) {
        nRep <- length(unique(Data_Total[[Bloc_Column]]))
        STP_seq <- seq(1, nRep, length_seq_combinations)
        
        AllComb <- gtools::combinations(
          n = as.numeric(length(STP_seq)) + 1,
          r = as.numeric(nGroups),
          v = as.numeric(c(seq(max(STP_seq), 1, -length_seq_combinations),0)),
          set = F,
          repeats.allowed = T
        )
      } else{
        AllComb <- gtools::combinations(
          n = as.numeric(length(seq_combinations)) + 1,
          r = as.numeric(nGroups),
          v = as.numeric(c(seq(max(seq_combinations), 1, -1),0)),
          set = F,
          repeats.allowed = T
        )
      }
      
      
      AllComb <- AllComb[-nrow(AllComb), ]
      colnames(AllComb) <- paste0("G", 1:nGroups)
      rownames(AllComb) <- seq.int(nrow(AllComb))
    } else {
      
      if (length_seq_combinations > 1) {
        warning(paste0("You choose 'length_seq_combinations' equal ", 
                       length_seq_combinations,
                       ", the combinations of selected trees per groups will be not sequential"))
      }
      
      if (is.null(seq_combinations)) {
        nRep <- length(unique(Data_Total[[Bloc_Column]]))
        NonSTP_seq <- seq(1, nRep, length_seq_combinations)
        
        AllComb <- gtools::combinations(
          n = as.numeric(length(NonSTP_seq)) + 1,
          r = as.numeric(nGroups),
          v = as.numeric(c(seq(max(NonSTP_seq), 1, -length_seq_combinations),0)),
          set = F,
          repeats.allowed = T
        )
      } else{
        AllComb <- gtools::combinations(
          n = as.numeric(length(seq_combinations)) + 1,
          r = as.numeric(nGroups),
          v = as.numeric(c(seq(max(seq_combinations), 1, -1),0)),
          set = F,
          repeats.allowed = T
        )
      }
      
      
      AllComb <- AllComb[-nrow(AllComb), ]
      colnames(AllComb) <- paste0("G", 1:nGroups)
      rownames(AllComb) <- seq.int(nrow(AllComb))
    }
    
    if (verbose) {
      message("Thinning strategy combinations configured.")
    }
    
  }, error = function(e) {
    print(paste("Error in thinning strategy configuration:", e$message))
    stop(e)
  })
  
  tryCatch({
    if (verbose) {
      message("Running BW strategies and collecting results...")
    }
  
    # Mutate Family_Data_Total and eliminate missing data
    Data_Total <- Data_Total |>
      mutate(Family = as.character(get(Family_Data_Total))) |> 
      filter(!is.na(get(Trait)))
    
    
    # BW strategies: Selection between and within families
    gg.BW = list()
    G.BW <- list()
    for (k in sort(unique(Data_Total$Group))) {
      for (i in seq_along(AllComb[, 1])) {
        g.bw <- Data_Total %>%
          filter(Group == k) %>%
          dplyr::group_by(Family, Plot) %>%
          arrange(Family, desc(get(BV_Column)))
        if (AllComb[i, which(colnames(AllComb) == k)] == 0) {
          gg.BW[[k]][[i]] <-
            dplyr::slice_head(g.bw, n = AllComb[i, which(colnames(AllComb) == k)])
        } else {
          gg.BW[[k]][[i]] <-
            dplyr::slice_head(g.bw, n = AllComb[i, which(colnames(AllComb) == k)])
          
        }
        
        if (nGroups == 2) {
          G.BW[[i]] <-
            rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]])
        } else {
          if (nGroups == 3) {
            G.BW[[i]] <-
              rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]], gg.BW[["G3"]][[i]])
          } else {
            G.BW[[i]] <-
              rbind(gg.BW[["G1"]][[i]], gg.BW[["G2"]][[i]], gg.BW[["G3"]][[i]], gg.BW[["G4"]][[i]])
          }
        }
        
        
        G.BW[[i]][, "BV"] <-
          G.BW[[i]][[BV_Column]] + mean(Data_Total[[Trait]], na.rm = T)
        G.BW[[i]][, "Strategy"] <- "BW"
      }
    }
    # Set names in the strategies
    G.BW <- setNames(G.BW, as.character(seq_along(G.BW)))
  }, error = function(e) {
    print(paste("Error in preparing BW strategies:", e$message))
    stop(e)
  }) 
    # BI Strategies: Selection of the best individuals regardless of family that he belongs
    ## Second group strategy
    if (is.null(IS)) {
      IC.TOP <- c(seq(0.01, 0.21, 0.02), 0.3, 0.4, 0.5)
      message(
        "The following IS are being used on the BI strategy:\n0.01 0.03 0.05 0.07 0.09 0.11 0.13 0.15 0.17 0.19 0.21 0.30 0.40 0.50"
      )
    } else {
      message(paste("The following IS are being used on the BI strategy:\n", IS))
      IC.TOP <- IS
    }
    
  tryCatch({
    if (verbose) {
      message("Running BI strategies and collecting results...")
    }
    
    SeleTop.BI <- list()
    for (i in seq_along(IC.TOP)) {
      seletop <- Data_Total %>%
        group_by(Family, Plot) %>%
        arrange(desc(get(BV_Column))) %>%
        head(n = round(nrow(Data_Total) * IC.TOP[i]))
      seletop$BV <- seletop[[BV_Column]] + mean(Data_Total[[Trait]], na.rm = T)
      SeleTop.BI[[i]] <- seletop
      SeleTop.BI[[i]]["Strategy"] <- "BI"
    }
    SeleTop.BI <-
      setNames(SeleTop.BI, as.character(seq(1 + length(G.BW), (length(G.BW) + length(SeleTop.BI)))))
  }, error = function(e) {
    print(paste("Error in preparing BI strategies:", e$message))
    stop(e)
  })
  
  tryCatch({
    if (verbose) {
      message("Combining BW and BI strategies...")
    }
    
    Strategies <- append(G.BW, SeleTop.BI)
    
    AllComb2 <- AllComb
    colnames(AllComb2) <- paste0("nID_fam_plot_G", 1:nGroups)
    
    fill_allcomb <-
      matrix(
        NA,
        nrow = length(Strategies) - (nrow(AllComb2)),
        ncol = nGroups,
        dimnames = list((nrow(AllComb2) + 1):length(Strategies),
                        colnames(AllComb2))
      )
    
    AllComb2 <- rbind(AllComb2, fill_allcomb)
    
    for (i in names(Strategies)) {
      Strategies[[i]][["N_Strategy"]] <- i
    }
  }, error = function(e) {
    print(paste("Error in combining BW and BI strategies:", e$message))
    stop(e)
  })
    # Effective Number (NE) and Genetic Gain with selection (GS)
  
  tryCatch({
    if (verbose) {
      message("Calculating NE and GS...")
    }
    
    GS.NE <- list()
    for (i in seq_along(Strategies)) {
      # GS
      gs <-
        round(mean(Strategies[[i]]$BV, na.rm = T) / mean(Data_Total[[Trait]], na.rm = T) - 1,
              3)
      # NE
      t <- table(as.factor(Strategies[[i]]$Family))
      # Averages of  Trees number for Mother Tree
      n <-
        length(Strategies[[i]][[id]]) / length(levels(as.factor(Strategies[[i]]$Family)))
      # Variance of Trees number for Mother Tree
      VarID <-
        sum((t - (n)) ^ 2) / (length(Strategies[[i]][[id]]) - 1)
      nf <- length(unique(Strategies[[i]]$Family))
      ne <- round((4 * nf * n) / (n + 3 + (VarID / n)), 1)
      gs.ne <-
        data.frame(
          GS = gs,
          `GS(%)` = scales::percent(gs, accuracy = 0.1),
          NE = ne
        ) %>%
        `colnames<-`(., c("GS", "GS(%)", "NE"))
      gs.ne$Trait <- Trait
      gs.ne$N_Strategy <- i
      gs.ne$N_Progeny <- nrow(Strategies[[i]])
      GS.NE[[i]] <-
        gs.ne[, c("N_Strategy", "N_Progeny", "Trait", "GS", "GS(%)", "NE")]
      GS.NE[[i]]$Strategy <- NA
    }
    
    # Joining BW and BI strategies: Pre final table
    for (i in 1:length(G.BW)) {
      GS.NE[[i]]$Strategy <- "BW"
    }
    for (i in (length(G.BW) + 1):(length(Strategies))) {
      GS.NE[[i]]$Strategy <- "BI"
    }
    GS.NE <- bind_rows(GS.NE)
    #GS.NE
  }, error = function(e) {
    print(paste("Error in calculate NE and GS:", e$message))
    stop(e)
  })
  
  tryCatch({
    if (verbose) {
      message("Summarise thinning strategies...")
    }
    
    # Summarise thinning strategies
    nFam_Strategy <- bind_rows(Strategies) |>
        data.frame() |>
        mutate(N_Strategy = as.numeric(N_Strategy)) |>
        group_by(N_Strategy, Group) |>
        summarise(nFam = n_distinct(Family)) |>
        pivot_wider(
          id_cols = N_Strategy,
          names_from = "Group",
          values_from = nFam,
          names_prefix = "nFam.",
          values_fill = 0
    )
  
    BV_fam <- left_join(BV_fam,
                        Data_Total[match(BV_fam$Family, Data_Total$Family), ] |>
                          dplyr::select(Family, Group),
                        by = "Family")
  }, error = function(e) {
    print(paste("Error in summarise thinning strategies:", e$message))
  })
  
  tryCatch({
    if (verbose) {
      message("Building the final table with thinning strategies...")
    }
    
    # Final table with thinning strategies
    Thinning <- left_join(GS.NE, nFam_Strategy, by = "N_Strategy") |>
      cbind(AllComb2)
    
    Thinning <- left_join(
      Thinning[, c("N_Strategy", "N_Progeny")] |>
        mutate(`N_Progeny(%)` = scales::percent(N_Progeny / nrow(Data_Total), accuracy = 0.1)),
      Thinning,
      by = c("N_Strategy", "N_Progeny")
    )
  }, error = function(e){
    print(paste("Error in building the final table with thinning strategies:", e$message))
    stop(e)
  })
  
  tryCatch({
    if (verbose) {
      message("Saving results...")
    }
    
    if (save_table_xlsx == TRUE) {
      require(openxlsx)
      write.xlsx(
        list(
          Thinning = Thinning,
          BV_fam = BV_fam
        ),
        file = paste0(
          "Thinning_Strategies_for_",
          Trait,
          "_trait_",
          export_id,
          Sys.Date(),
          ".xlsx"
        ),
        overwrite = T,
        asTable = T
      )
      write.xlsx(
        Strategies,
        file = paste0(
          "Strategies_",
          Trait,
          "_trait_",
          export_id,
          Sys.Date(),
          ".xlsx"
        ),
        overwrite = T,
        asTable = T
      )
    }
    
    if (save_plot_rank == TRUE) {
      ggsave(filename = paste0("families_rank_", export_id, format_plot), plot = families_rank)
    }
  }, error = function(e) {
    print(paste("Error in saving results:", e$message))
    stop(e)
  })
    return(
      list(
        Thinning = Thinning,
        Strategies = Strategies,
        BV_Fam = BV_fam,
        ggPlot_families_rank = families_rank
      )
    )
  if (verbose) {
    message("Finished running Thinning_BreedR function.")
  }
  
}

