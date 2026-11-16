
#' Thinning_BreedR: Thinning Strategies Based on Genetic Values
#'
#' This function implements thinning strategies between and within families based on additive genetic values (BV) 
#' for a given trait. It identifies family groups with different genetic performance, generates ranking plots, 
#' and calculates genetic gain (GS) and effective number (NE).
#'
#' @param BV_Column Name of the column containing additive genetic values (default = "a_total").
#' @param Trait Name of the trait of interest.
#' @param BV_fam Data frame containing genetic values per family.
#' @param Data_Total Data frame with complete individual data.
#' @param Family_Data_Total Name of the family column in `Data_Total`.
#' @param Bloc_Column Name of the experimental block column.
#' @param Family_Column Name of the family column in `BV_fam`.
#' @param nGroups Number of groups to form (2, 3, or 4).
#' @param label.group.y Vector of Y-axis positions for group annotations in the plot.
#' @param Plot.Rank Logical. If TRUE, plots the family ranking graph.
#' @param save_plot_rank Logical. If TRUE, saves the ranking plot.
#' @param IS Vector of selection proportions for the BI (individual) strategy.
#' @param id Name of the individual ID column.
#' @param nGroups3 Direction for 3-group split ("left" or "right").
#' @param STP Logical. If TRUE, treats blocks as experimental units.
#' @param n.dodge_plot1 Integer. Number of lines for X-axis labels.
#' @param angle_plot1 Angle of X-axis labels.
#' @param nt_alpha1 Parameters for inflection detection in best families.
#' @param nt_alpha2 Parameters for inflection detection in worst families.
#' @param finfl_Value Vector with "max" or "min" to define inflection points.
#' @param additional_layer_plot1 Additional ggplot2 layers for the ranking plot.
#' @param seq_combinations Custom sequence for thinning combinations.
#' @param length_seq_combinations Length of the thinning combination sequence.
#' @param save_table_xlsx Logical. If TRUE, saves output tables as Excel files.
#' @param export_id Suffix for exported file names.
#' @param format_plot Format of the saved plot (e.g., ".tiff").
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param Tree_BV Optional data frame containing individual breeding values (BLUPs).
#' @param TreeID_Column Individual identifier column used to join Tree_BV and Data_Total.
#' @param Tree_BV_Column Individual BLUP column in Tree_BV.
#' @param target_NE Numeric vector of requested effective population sizes; NULL disables the filter.
#' @param target_NE_tol Absolute Ne tolerance, either a scalar or one value per target.
#' @param n_best_strategies Maximum number of retained strategies per Ne target.
#' @param max_candidates Maximum number of BW combinations screened.
#' @param rebuild_selected Logical. Rebuild full individual-level records only for retained strategies.
#' @param seed Random seed used when sampling the candidate-combination space.
#'
#' @details Candidate strategies are screened using compact numerical summaries.
#' Effective population size is calculated from family contribution counts,
#' not from the magnitude or mean of individual BLUPs. Negative BLUPs are valid.
#' GS_absolute is the mean selected BLUP in trait units; GS is its ratio to the
#' phenotypic mean when that mean is nonzero. When the candidate space exceeds
#' max_candidates, a sampled search is used and global optimality is not guaranteed.
#'
#' @return A list containing:
#' \describe{
#'   \item{Thinning}{Final table with thinning strategies, GS, NE, and number of families per group.}
#'   \item{Strategies}{List of selected individuals per strategy.}
#'   \item{BV_Fam}{Data frame with families and their assigned groups.}
#'   \item{ggPlot_families_rank}{Family ranking plot.}
#' }
#'
#' @import dplyr ggplot2 gtools RootsExtremaInflections openxlsx
#' @export

#---------------------- Function for thinning strategies --------------------------
Thinning_BreedR <- function(BV_Column = "a_total",
                            Trait,
                            BV_fam,
                            Data_Total,
                            Family_Data_Total = "Family",
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
                            verbose = FALSE,
                            Tree_BV = NULL,
                            TreeID_Column = "TreeID",
                            Tree_BV_Column = "BV",
                            target_NE = NULL,
                            target_NE_tol = 10,
                            n_best_strategies = 20L,
                            max_candidates = 50000L,
                            rebuild_selected = TRUE,
                            seed = 123L) {
  if (verbose) {
    message("Starting Thinning_BreedR function...")
  }
  
  if (!require("RootsExtremaInflections")) {
    install.packages("RootsExtremaInflections")
  }
  
  if (!require("gtools")) {
    install.packages("gtools")
  }
  
  if (!nGroups %in% 2:4) {
    stop("nGroups should be: 2, 3 or 4")
  }
  
  if (any(!finfl_Value %in% c("max", "min"))) {
    stop("finfl_value should be a vector of lenght 2 of: 'max' or 'min'")
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
      message("BV_fam data frame processed without problem.")
    }
    
  }, error = function(e) {
    print(paste("Error processing BV_fam data frame:", e$message))
    stop(e)
  })
  # # Rename target column to 'Family'
  # Data_Total <- Data_Total |> 
  #   rename_with(~ "Family", all_of(Family_Data_Total))
  
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
      message("Identifying families with positive, zero and negative BV...")
    }
    
    # Families with positive BV
    MajorBV_fam <- BV_fam |>
      filter(get(BV_Column) >= 0)
    
    # Families with BV near of zero
    Fam_Zero <- which.min(MajorBV_fam[[BV_Column]])
    # Families with negative BV
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
    
    # Found inflexion point in the best families
    Inflex1 <- RootsExtremaInflections::inflexi(
      x = 1:nrow(MajorBV_fam),
      y = MajorBV_fam[[BV_Column]],
      nt = nt_alpha1[1],
      alpha = nt_alpha1[2],
      i1 = 1,
      i2 = nrow(MajorBV_fam),
      plots = F,
    )
    
    if (finfl_Value[1] == "max") {
      Inflexi_major <- max(Inflex1$finfl, na.rm = T)
    } else {
      Inflexi_major <- min(Inflex1$finfl, na.rm = T)
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
    
    # Found inflexion point in the worst families
    Inflex2 <- RootsExtremaInflections::inflexi(
      x = 1:nrow(MinorBV_fam),
      y = MinorBV_fam[[BV_Column]],
      nt = nt_alpha2[1],
      i1 = nt_alpha2[2],
      i2 = nrow(MinorBV_fam),
      plots = F,
    )
    
    if (finfl_Value[2] == "max") {
      Inflexi_minor <- Fam_Zero + max(Inflex2$finfl, na.rm = T)
    } else {
      Inflexi_minor <- Fam_Zero + min(Inflex2$finfl, na.rm = T)
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
    
    # Plot BV vs Rank from BV
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
    Data_Total$Group <- NA_character_
    # Family IDs may be factors in BV_fam and character strings in Data_Total.
    # Compare normalized character keys without changing the family ranking order.
    family_group_key <- trimws(as.character(Data_Total[[Family_Data_Total]]))
    family_rank_key <- trimws(as.character(BV_fam$Family))
    if (!any(family_group_key %in% family_rank_key)) {
      stop("No family identifiers overlap between Data_Total and BV_fam. Check Family_Data_Total and Family_Column.")
    }
    
    # G1
    if (nGroups == 2) {
      Data_Total[family_group_key %in% family_rank_key[1:Fam_Zero]), "Group"] <-
        "G1"
    } else {
      Data_Total[family_group_key %in% family_rank_key[1:Inflexi_major]), "Group"] <-
        "G1"
    }
    
    # G2
    if (nGroups == 2) {
      Data_Total[family_group_key %in% family_rank_key[(Fam_Zero + 1):nrow(BV_fam)]), "Group"] <-
        "G2"
    } else {
      Data_Total[family_group_key %in% family_rank_key[(Inflexi_major + 1):Fam_Zero]), "Group"] <-
        "G2"
    }
    
    # G3
    if (nGroups == 3) {
      Data_Total[family_group_key %in% family_rank_key[(Fam_Zero + 1):nrow(BV_fam)]), "Group"] <-
        "G3"
    } else {
      if (nGroups == 4) {
        Data_Total[family_group_key %in% family_rank_key[(Fam_Zero + 1):Inflexi_minor]), "Group"] <-
          "G3"
      }
      
    }
    
    # G4
    if (nGroups == 4) {
      Data_Total[family_group_key %in% family_rank_key[(Inflexi_minor + 1):nrow(BV_fam)]), "Group"] <-
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
  

  # Screen compact records first, reconstruct selected individuals only at the end.
  # Ne is based on family contribution counts, not on average BLUP values.
  # Resolve user-supplied column names before any indexing or block grouping.
  # Match case-insensitively only when the match is unique.
  resolve_column <- function(requested, available, argument) {
    if (requested %in% available) return(requested)
    matches <- available[tolower(available) == tolower(requested)]
    if (length(matches) == 1L) {
      warning(sprintf("%s: using column '%s' instead of '%s'.",
                      argument, matches, requested))
      return(matches)
    }
    stop(sprintf("%s column '%s' not found. Available columns: %s",
                 argument, requested, paste(available, collapse = ", ")))
  }
  Bloc_Column <- resolve_column(Bloc_Column, names(Data_Total), "Bloc_Column")
  Family_Data_Total <- resolve_column(Family_Data_Total, names(Data_Total),
                                      "Family_Data_Total")
  Trait <- resolve_column(Trait, names(Data_Total), "Trait")
  if (!id %in% names(Data_Total)) {
    # The user's experimental dataset stores individual IDs in TreeID, not ID.
    alternative_id <- if (TreeID_Column %in% names(Data_Total)) TreeID_Column else {
      matching <- names(Data_Total)[tolower(names(Data_Total)) ==
                                      tolower(TreeID_Column)]
      if (length(matching) == 1L) matching else character()
    }
    if (length(alternative_id) == 1L) {
      warning(sprintf("Individual identifier '%s' not found; using '%s'.",
                      id, alternative_id))
      id <- alternative_id
    } else {
      stop(sprintf("Individual identifier '%s' not found. Pass id = 'TreeID' or another valid column.", id))
    }
  }
  if (!is.null(Tree_BV)) {
    TreeID_Column <- resolve_column(TreeID_Column, names(Data_Total),
                                    "TreeID_Column (Data_Total)")
    TreeID_Column_BV <- resolve_column(TreeID_Column, names(Tree_BV),
                                       "TreeID_Column (Tree_BV)")
    Tree_BV_Column <- resolve_column(Tree_BV_Column, names(Tree_BV),
                                     "Tree_BV_Column")
  }
  if (!Trait %in% names(Data_Total)) stop("Trait not found in Data_Total.")
  Data_Total$Family <- as.character(Data_Total[[Family_Data_Total]])
  if (!is.null(Tree_BV)) {
    if (!TreeID_Column_BV %in% names(Tree_BV) ||
        !Tree_BV_Column %in% names(Tree_BV) ||
        !TreeID_Column %in% names(Data_Total)) stop("Invalid Tree_BV columns.")
    tb <- Tree_BV[, c(TreeID_Column_BV, Tree_BV_Column), drop = FALSE]
    names(tb) <- c(".tree_id", ".individual_blup")
    tb$.tree_id <- as.character(tb$.tree_id)
    if (anyDuplicated(tb$.tree_id)) stop("Duplicated TreeIDs in Tree_BV.")
    Data_Total$.tree_id <- as.character(Data_Total[[TreeID_Column]])
    Data_Total <- dplyr::left_join(Data_Total, tb, by = ".tree_id")
    blup_col <- ".individual_blup"
  } else {
    warning("Tree_BV missing: falling back to BV_Column (family value if merged).")
    blup_col <- BV_Column
  }
  # Diagnose each exclusion separately before filtering the data.
  # This distinguishes unmatched tree BLUPs, missing trait values, and missing groups.
  # Never replace missing BLUPs with zero: missing values indicate absent predictions.
  if (!is.numeric(Data_Total[[blup_col]])) {
    stop(sprintf("Individual BLUP column '%s' must be numeric; current class: %s",
                 blup_col, paste(class(Data_Total[[blup_col]]), collapse = ", ")))
  }
  if (!is.numeric(Data_Total[[Trait]])) {
    stop(sprintf("Trait column '%s' must be numeric; current class: %s",
                 Trait, paste(class(Data_Total[[Trait]]), collapse = ", ")))
  }
  usable_blup <- is.finite(Data_Total[[blup_col]])
  usable_trait <- is.finite(Data_Total[[Trait]])
  usable_group <- !is.na(Data_Total$Group) & Data_Total$Group %in%
    paste0("G", seq_len(nGroups))
  usable_all <- usable_blup & usable_trait & usable_group
  family_data_keys <- trimws(as.character(Data_Total[[Family_Data_Total]]))
  family_bv_keys <- trimws(as.character(BV_fam$Family))
  family_overlap <- !is.na(family_data_keys) & family_data_keys %in% family_bv_keys
  # Print failure diagnostics even when verbose = FALSE.
  if (isTRUE(verbose) || !any(usable_all)) {
    message(sprintf(
      "Screening counts: total=%d; finite_BLUP=%d; finite_trait=%d; assigned_group=%d; family_match=%d; retained=%d",
      nrow(Data_Total), sum(usable_blup), sum(usable_trait),
      sum(usable_group), sum(family_overlap), sum(usable_all)
    ))
    message(sprintf(
      "Overlaps: BLUP_and_trait=%d; BLUP_and_group=%d; trait_and_group=%d",
      sum(usable_blup & usable_trait), sum(usable_blup & usable_group),
      sum(usable_trait & usable_group)
    ))
    message("Group counts: ", paste(
      capture.output(print(table(Data_Total$Group, useNA = "ifany"))),
      collapse = " "
    ))
    if (!is.null(Tree_BV)) {
      matched_ids <- !is.na(Data_Total$.tree_id) &
        Data_Total$.tree_id %in% tb$.tree_id
      message(sprintf(
        "TreeID join: matched=%d/%d; finite_joined_BLUP=%d",
        sum(matched_ids), nrow(Data_Total), sum(usable_blup)
      ))
      if (any(!matched_ids)) {
        message("Example unmatched TreeIDs: ", paste(
          utils::head(unique(Data_Total$.tree_id[!matched_ids]), 5L),
          collapse = ", "
        ))
      }
    }
    if (any(!family_overlap)) {
      message("Example unmatched families: ", paste(
        utils::head(unique(family_data_keys[!family_overlap]), 5L),
        collapse = ", "
      ))
    }
  }
  if (!any(usable_all)) {
    stop(sprintf(
      paste0("No usable individual records. total=%d; finite_BLUP=%d; ",
             "finite_trait=%d; assigned_group=%d; family_match=%d; ",
             "BLUP_and_trait=%d; BLUP_and_group=%d. ",
             "Review TreeID, family/progenie keys, BLUP and MAI columns."),
      nrow(Data_Total), sum(usable_blup), sum(usable_trait),
      sum(usable_group), sum(family_overlap),
      sum(usable_blup & usable_trait), sum(usable_blup & usable_group)
    ))
  }
  Data_Total <- Data_Total[usable_all, , drop = FALSE]
  if (!is.numeric(Data_Total[[blup_col]])) stop("BLUP must be numeric.")
  if (anyDuplicated(as.character(Data_Total[[id]])))
    stop("Individual IDs must be unique.")
  if (!is.null(target_NE)) {
    if (!is.numeric(target_NE) || any(!is.finite(target_NE)) || any(target_NE <= 0))
      stop("target_NE must contain positive finite values.")
    if (length(target_NE_tol) == 1L) target_NE_tol <- rep(target_NE_tol, length(target_NE))
    if (length(target_NE_tol) != length(target_NE) ||
        any(!is.finite(target_NE_tol)) || any(target_NE_tol < 0))
      stop("target_NE_tol must contain nonnegative tolerances.")
  }
  if (n_best_strategies < 1 || max_candidates < 1)
    stop("n_best_strategies and max_candidates must be positive.")
  n_best_strategies <- as.integer(n_best_strategies)
  max_candidates <- as.integer(max_candidates)
  groups <- paste0("G", seq_len(nGroups))
  cells <- lapply(groups, function(g) {
    ix <- which(Data_Total$Group == g)
    if (!length(ix)) return(list())
    z <- split(ix, interaction(Data_Total$Family[ix], Data_Total$Plot[ix], drop = TRUE))
    lapply(z, function(ii) ii[order(-Data_Total[[blup_col]][ii],
                                       as.character(Data_Total[[id]][ii]))])
  })
  names(cells) <- groups
  nrep <- length(unique(Data_Total[[Bloc_Column]]))
  levels_n <- if (is.null(seq_combinations)) {
    c(0L, seq.int(1L, max(1L, nrep), by = length_seq_combinations))
  } else c(0L, as.integer(seq_combinations))
  levels_n <- sort(unique(levels_n))
  if (anyNA(levels_n) || any(levels_n < 0)) stop("Invalid selection intensities.")
  base_n <- length(levels_n)
  total <- as.double(base_n)^nGroups - 1
  if (!is.finite(total)) stop("Combination space overflow.")
  decode <- function(code) {
    z <- integer(nGroups)
    for (k in seq_len(nGroups)) {
      z[k] <- levels_n[(code %% base_n) + 1L]
      code <- floor(code / base_n)
    }
    z
  }
  select_bw <- function(combo) {
    unlist(lapply(seq_len(nGroups), function(k) {
      if (!combo[k]) return(integer())
      unlist(lapply(cells[[k]], utils::head, n = combo[k]), use.names = FALSE)
    }), use.names = FALSE)
  }
  calc_ne <- function(family) {
    t <- as.numeric(table(family))
    N <- sum(t)
    if (!N) return(NA_real_)
    nf <- length(t)
    mean_n <- N / nf
    var_n <- if (N > 1) sum((t - mean_n)^2) / (N - 1) else 0
    round(4 * nf * mean_n / (mean_n + 3 + var_n / mean_n), 1)
  }
  base_trait <- mean(Data_Total[[Trait]], na.rm = TRUE)
  summary_one <- function(index, kind, key, combo) {
    if (!length(index)) return(NULL)
    mean_blup <- mean(Data_Total[[blup_col]][index])
    ne <- calc_ne(Data_Total$Family[index])
    if (!is.finite(ne)) return(NULL)
    gs <- if (abs(base_trait) > .Machine$double.eps)
      mean_blup / base_trait else NA_real_
    data.frame(Candidate_ID = key, Strategy = kind, NE = ne,
      GS = round(gs, 6), GS_absolute = mean_blup,
      N_Progeny = length(index), nFamily = length(unique(Data_Total$Family[index])),
      as.list(stats::setNames(combo,
        paste0("nID_fam_plot_G", seq_len(nGroups)))),
      check.names = FALSE)
  }
  # Do not allocate a Cartesian-product matrix or store candidate data frames.
  ncheck <- as.integer(min(total, max_candidates))
  if (total > max_candidates) {
    warning("Sampling candidate combinations; global optimum is not guaranteed.")
    set.seed(seed)
    codes <- unique(c(seq_len(min(ncheck, 100L)),
                      floor(stats::runif(ncheck * 2L, 0, total)) + 1))
    codes <- head(codes, ncheck)
  } else codes <- seq_len(ncheck)
  rows <- vector("list", length(codes) + 100L)
  nrows <- 0L
  for (code in codes) {
    one <- summary_one(select_bw(decode(code)), "BW", paste0("BW_", code), decode(code))
    if (!is.null(one)) {nrows <- nrows + 1L; rows[[nrows]] <- one}
  }
  is_vec <- if (is.null(IS)) c(seq(0.01, 0.21, 0.02), 0.3, 0.4, 0.5) else IS
  if (any(!is.finite(is_vec)) || any(is_vec <= 0 | is_vec > 1))
    stop("IS must be proportions in (0,1].")
  order_bi <- order(-Data_Total[[blup_col]], as.character(Data_Total[[id]]))
  for (i in seq_along(is_vec)) {
    n <- max(1L, min(nrow(Data_Total), round(nrow(Data_Total) * is_vec[i])))
    one <- summary_one(order_bi[seq_len(n)], "BI", paste0("BI_",i),
                       rep(NA_integer_, nGroups))
    if (!is.null(one)) {nrows <- nrows + 1L; rows[[nrows]] <- one}
  }
  if (!nrows) stop("No valid strategies.")
  screened <- dplyr::bind_rows(rows[seq_len(nrows)])
  screened$Target_NE <- NA_real_
  screened$NE_Distance <- NA_real_
  if (!is.null(target_NE)) {
    dist <- abs(outer(screened$NE, target_NE, "-"))
    keep <- rowSums(sweep(dist, 2L, target_NE_tol, "<=")) > 0L
    screened <- screened[keep, , drop = FALSE]
    dist <- dist[keep, , drop = FALSE]
    if (nrow(screened)) {
      best <- max.col(-dist, ties.method = "first")
      screened$Target_NE <- target_NE[best]
      screened$NE_Distance <- dist[cbind(seq_len(nrow(screened)), best)]
    }
  }
  if (nrow(screened)) {
    screened <- screened[order(-screened$GS_absolute, screened$NE_Distance), , drop = FALSE]
    if (is.null(target_NE)) {
      screened <- utils::head(screened, n_best_strategies)
    } else {
      screened <- dplyr::bind_rows(lapply(unique(screened$Target_NE), function(x)
        utils::head(screened[screened$Target_NE == x, , drop = FALSE],
                    n_best_strategies)))
    }
  }
  screened$N_Strategy <- seq_len(nrow(screened))
  screened$Trait <- Trait
  screened[["GS(%)"]] <- scales::percent(screened$GS, accuracy = 0.1)
  screened[["N_Progeny(%)"]] <- scales::percent(
    screened$N_Progeny / nrow(Data_Total), accuracy = 0.1)
  Thinning <- screened
  Strategies <- list()
  if (isTRUE(rebuild_selected)) {
    for (i in seq_len(nrow(screened))) {
      candidate <- screened$Candidate_ID[i]
      ix <- if (screened$Strategy[i] == "BW") {
        select_bw(decode(as.double(sub("BW_", "", candidate))))
      } else {
        j <- as.integer(sub("BI_", "", candidate))
        n <- max(1L, min(nrow(Data_Total), round(nrow(Data_Total) * is_vec[j])))
        order_bi[seq_len(n)]
      }
      item <- Data_Total[ix, , drop = FALSE]
      item$BV <- item[[blup_col]] + base_trait
      item$Strategy <- screened$Strategy[i]
      item$N_Strategy <- i
      Strategies[[as.character(i)]] <- item
    }
  }
  BV_fam$Group <- Data_Total$Group[match(as.character(BV_fam$Family),
                                         Data_Total$Family)]
  if (isTRUE(save_table_xlsx)) {
    if (!requireNamespace("openxlsx", quietly = TRUE))
      stop("Install 'openxlsx' for Excel export.")
    openxlsx::write.xlsx(list(Thinning = Thinning, BV_fam = BV_fam),
      file = paste0("Thinning_Strategies_for_",Trait,"_trait_",export_id,
                    Sys.Date(),".xlsx"), overwrite = TRUE)
    if (length(Strategies)) openxlsx::write.xlsx(Strategies,
      file = paste0("Strategies_",Trait,"_trait_",export_id,
                    Sys.Date(),".xlsx"), overwrite = TRUE)
  }
  if (isTRUE(save_plot_rank)) ggplot2::ggsave(
    filename = paste0("families_rank_",export_id,format_plot),
    plot = families_rank)
  if (verbose) message("Finished running Thinning_BreedR function.")
  list(Thinning = Thinning, Strategies = Strategies, BV_Fam = BV_fam,
       ggPlot_families_rank = families_rank)
}
