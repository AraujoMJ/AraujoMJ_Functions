
#' Thinning_BreedR: Thinning Strategies Based on Genetic Values
#'
#' This function implements thinning strategies between and within families based on additive genetic values (BV) 
#' for a given trait. It identifies family groups with different genetic performance, generates ranking plots, 
#' and calculates genetic gain (GS) and effective number (NE).
#'
#' @param BV_Column Name of the column containing additive genetic values (default = "a_total").
#' @param Trait Name of the trait of interest.
#' @param BV_fam Data frame containing genetic values per family.
#' @param BV_fam_full Optional unfiltered family breeding-value data for computing
#'   inflection boundaries and G1..Gn across the entire population. If NULL,
#'   BV_fam is used, and it must itself be unfiltered for global groups.
#'   Filtering Data_Total does not affect the ranking or group boundaries.
#' @param Data_Total Data frame with complete individual data.
#' @param Family_Data_Total Name of the family column in `Data_Total`.
#' @param Bloc_Column Name of the experimental block column.
#' @param Family_Column Name of the family column in `BV_fam`.
#' @param nGroups Number of groups to form (2 through 6).
#' @param positive_groups Number of groups among nonnegative family BLUPs;
#'   NULL uses 1/1 for 2 groups, 2/1 or 1/2 for 3, 2/2 for 4,
#'   3/2 or 2/3 for 5, and 3/3 for 6. Both sides permit 1 to 3 groups.
#' @param label.group.y Vector of Y-axis positions for group annotations in the plot.
#' @param Plot.Rank Logical. If TRUE, plots the family ranking graph.
#' @param save_plot_rank Logical. If TRUE, saves the ranking plot.
#' @param IS Vector of BI proportions; NULL evaluates 1-50% by 1%, 55-100% by 5%.
#' @param id Name of the individual ID column.
#' @param nGroups3 Direction of the extra group for 3 groups ("left" or "right").
#' @param nGroups5 Direction of the extra group for 5 groups ("left":
#'   3 positive / 2 negative; "right": 2 positive / 3 negative).
#' @param STP Logical. If TRUE, treats blocks as experimental units.
#' @param n.dodge_plot1 Integer. Number of lines for X-axis labels.
#' @param angle_plot1 Angle of X-axis labels.
#' @param nt_alpha1 Parameters for inflection detection in best families.
#' @param nt_alpha2 Parameters for inflection detection in worst families.
#' @param finfl_Value Two-element vector specifying which candidate inflection
#'   to prefer on the nonnegative and negative family sides, respectively.
#'   A second split (for three groups on one side) uses another detected
#'   inflection when available; otherwise the largest remaining segment is
#'   divided at its midpoint and reported as a fallback.
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
#' @param all_scenarios Logical. If TRUE, exhaustively evaluate and retain every BW/BI scenario, ignoring target_NE, max_candidates, and n_best_strategies. Use rebuild_selected = FALSE to avoid retaining all individual records.
#' @param rebuild_selected Logical. If TRUE, rebuild full records; if FALSE, return compact TreeID/family/group/BLUP records for retained strategies.
#' @param seed Random seed used when sampling the candidate-combination space.
#' @param min_group_families Minimum represented families in each G1..Gn group.
#'   A scalar applies to all groups; a vector supplies group-specific minima.
#'   Default is one family in G1 and zero in the remaining groups.
#' @param min_group_trees Minimum selected trees in each group (default zero).
#' @param bw_selection_unit One of 'auto', 'plot', or 'family'. Auto switches to
#'   family-level selection for singleton family/plot cells. Plot preserves the
#'   original per-plot selection semantics.
#' @param pareto If TRUE, prioritize nondominated GS_absolute vs NE solutions,
#'   retaining the best remaining feasible solutions only when needed to fill
#'   n_best_strategies. Full results are available with all_scenarios = TRUE.
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
#'   \item{Group_Boundaries}{Global rank cutpoints and their origin: inflection, zero crossing, or midpoint fallback.}
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
                            nGroups5 = "left",
                            positive_groups = NULL,
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
                            all_scenarios = FALSE,
                            rebuild_selected = TRUE,
                            seed = 123L,
                            min_group_families = NULL,
                            min_group_trees = NULL,
                            pareto = TRUE,
                            bw_selection_unit = "auto",
                            BV_fam_full = NULL) {
  if (verbose) {
    message("Starting Thinning_BreedR function...")
  }
  
  if (!require("RootsExtremaInflections")) {
    install.packages("RootsExtremaInflections")
  }
  
  if (!require("gtools")) {
    install.packages("gtools")
  }
  
  if (length(nGroups) != 1L || is.na(nGroups) ||
      !is.numeric(nGroups) || nGroups != floor(nGroups) ||
      !nGroups %in% 2:6) {
    stop("nGroups must be an integer from 2 to 6.")
  }
  if (!nGroups3 %in% c("left", "right") ||
      !nGroups5 %in% c("left", "right"))
    stop("nGroups3 and nGroups5 must be 'left' or 'right'.")
  if (is.null(positive_groups)) {
    positive_groups <- switch(as.character(nGroups),
      "2" = 1L,
      "3" = if (nGroups3 == "left") 2L else 1L,
      "4" = 2L,
      "5" = if (nGroups5 == "left") 3L else 2L,
      "6" = 3L)
  }
  if (length(positive_groups) != 1L ||
      !is.numeric(positive_groups) || !is.finite(positive_groups) ||
      positive_groups != floor(positive_groups) ||
      !positive_groups %in% 1:3 ||
      !(nGroups - positive_groups) %in% 1:3)
    stop("positive_groups must give 1 to 3 groups on each BLUP-sign side.")
  positive_groups <- as.integer(positive_groups)
  negative_groups <- as.integer(nGroups - positive_groups)
  
  if (any(!finfl_Value %in% c("max", "min"))) {
    stop("finfl_value should be a vector of lenght 2 of: 'max' or 'min'")
  }
  
  # Compute global group definitions from complete family predictions before
  # filtering individuals by experiment, type, trait or available tree BLUPs.
  # Data_Total may already be filtered by the calling code: it only defines
  # the eligible population for BW/BI screening, not family-group boundaries.
  if (!is.null(BV_fam_full)) {
    if (!is.data.frame(BV_fam_full)) stop("BV_fam_full must be a data frame.")
    BV_fam <- BV_fam_full
    if (verbose) message("Using unfiltered BV_fam_full for global ranking and groups.")
  } else if (verbose) {
    message("Using BV_fam for ranking: provide unfiltered family values here, ",
            "or pass BV_fam_full, to prevent filter-dependent boundaries.")
  }
  # Start of main function logic
  tryCatch({
    # Process BV_fam data frame
    if (verbose) {
      message("Processing BV_fam data frame...")
    }
    
    if (!all(c(Family_Column, BV_Column) %in% names(BV_fam)))
      stop("The global family table must contain Family_Column and BV_Column.")
    BV_fam <- BV_fam |>
      rename_with(~ "Family", all_of(Family_Column)) |>
      mutate(Family = trimws(as.character(Family))) |>
      arrange(desc(.data[[BV_Column]]))
    if (anyNA(BV_fam$Family) || any(!nzchar(BV_fam$Family)) ||
        anyDuplicated(BV_fam$Family))
      stop("Global family identifiers must be unique, nonmissing values.")
    BV_fam$Family <- factor(BV_fam$Family, levels = BV_fam$Family)
    
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
    
    if (!Family_Data_Total %in% names(Data_Total))
      stop("Family_Data_Total column is missing from Data_Total.")
    Data_Total$Family <- trimws(as.character(Data_Total[[Family_Data_Total]]))
    # Avoid a stale family BV column or a duplicate BV.x/BV.y pair.
    if (BV_Column %in% names(Data_Total))
      Data_Total[[BV_Column]] <- NULL
    family_lookup <- BV_fam[, c("Family", BV_Column), drop = FALSE]
    family_lookup$Family <- as.character(family_lookup$Family)
    Data_Total <- dplyr::left_join(Data_Total, family_lookup, by = "Family")
    
    if (verbose) {
      message("Data merged successfully.")
    }
    
  }, error = function(e) {
    print(paste("Error merging Data_Total and BV_fam:", e$message))
    stop(e)
  })
  
  # Compute every rank boundary on the complete, unfiltered family set.
  # The zero crossing is a mandatory boundary between BLUP-sign sides.
  family_values <- BV_fam[[BV_Column]]
  if (!is.numeric(family_values) || any(!is.finite(family_values)))
    stop("Family breeding values must be finite numeric values.")
  n_family <- nrow(BV_fam)
  n_positive <- sum(family_values >= 0)
  n_negative <- n_family - n_positive
  if (n_positive < positive_groups || n_negative < negative_groups)
    stop(sprintf("Insufficient families for %d positive and %d negative groups.",
                 positive_groups, negative_groups))
  if (length(finfl_Value) != 2L || anyNA(finfl_Value) ||
      any(!finfl_Value %in% c("min", "max")))
    stop("finfl_Value must contain two entries: 'min' or 'max'.")

  inflection_candidates <- function(values, tuning, side) {
    n <- length(values)
    if (n < 3L) return(integer())
    fit <- tryCatch(
      RootsExtremaInflections::inflexi(
        x = seq_len(n), y = values, nt = tuning[1L],
        alpha = tuning[2L], i1 = 1L, i2 = n, plots = FALSE),
      error = function(e) {
        warning(sprintf("%s inflection estimation failed: %s",
                        side, conditionMessage(e)))
        NULL
      })
    if (is.null(fit)) return(integer())
    candidates <- suppressWarnings(as.numeric(unlist(fit$finfl)))
    candidates <- sort(unique(as.integer(round(candidates[is.finite(candidates)]))))
    candidates[candidates >= 1L & candidates < n]
  }

  # Keep the legacy primary min/max inflection choice. For a third group
  # on the same side, prefer a distinct additional detected inflection.
  # If the method provides no second inflection, bisect the longest
  # available rank segment. The fallback is not an inferred inflection.
  split_side <- function(values, n_side_groups, tuning, rule, side) {
    n <- length(values)
    candidates <- if (n_side_groups > 1L)
      inflection_candidates(values, tuning, side) else integer()
    selected <- integer()
    methods <- character()
    while (length(selected) < n_side_groups - 1L) {
      available <- setdiff(candidates, selected)
      if (!length(selected) && length(available)) {
        chosen <- if (rule == "min") min(available) else max(available)
        method <- "inflexion"
      } else if (length(available)) {
        edges <- c(0L, sort(selected), n)
        usable <- available[vapply(available, function(v) {
          j <- findInterval(v - 0.5, edges, rightmost.closed = TRUE)
          (v - edges[j] >= 1L) && (edges[j + 1L] - v >= 1L)
        }, logical(1))]
        if (length(usable)) {
          if (rule == "min") chosen <- min(usable) else chosen <- max(usable)
          method <- "inflexion"
        } else available <- integer()
      }
      if (!length(available) || (!length(selected) && !length(candidates))) {
        edges <- c(0L, sort(selected), n)
        spans <- diff(edges)
        eligible <- which(spans >= 2L)
        if (!length(eligible)) stop("Cannot form nonempty family groups.")
        span <- eligible[which.max(spans[eligible])]
        chosen <- edges[span] + floor(spans[span] / 2L)
        method <- "midpoint_fallback"
        warning(sprintf("%s side: no additional valid inflection; ",
                        side), sprintf("splitting at local rank %d.", chosen))
      }
      selected <- c(selected, as.integer(chosen))
      methods <- c(methods, method)
    }
    ord <- order(selected)
    list(cuts = as.integer(selected[ord]), methods = methods[ord],
         candidates = candidates)
  }
  pos_split <- split_side(family_values[seq_len(n_positive)],
                          positive_groups, nt_alpha1,
                          finfl_Value[1L], "Positive")
  neg_split <- split_side(family_values[n_positive + seq_len(n_negative)],
                          negative_groups, nt_alpha2,
                          finfl_Value[2L], "Negative")
  cuts <- c(pos_split$cuts, n_positive, n_positive + neg_split$cuts)
  cut_methods <- c(pos_split$methods, "zero_crossing", neg_split$methods)
  cuts <- as.integer(cuts)
  if (length(cuts) != nGroups - 1L ||
      any(diff(c(0L, cuts, n_family)) <= 0L))
    stop("Invalid group boundaries: check the family count and inflections.")
  Group_Boundaries <- data.frame(
    Cut_Rank = cuts, Method = cut_methods,
    Side = c(rep("positive", length(pos_split$cuts)), "zero",
             rep("negative", length(neg_split$cuts))),
    stringsAsFactors = FALSE)
  # Use the same boundary vector for plotting and tree group assignment.
  rank_group <- paste0("G", findInterval(
    seq_len(n_family) - 0.5, cuts) + 1L)
  BV_fam$Group <- rank_group
  rank_families <- trimws(as.character(BV_fam$Family))
  tree_families <- trimws(as.character(Data_Total[[Family_Data_Total]]))
  family_index <- match(tree_families, rank_families)
  if (!any(!is.na(family_index)))
    stop("No family identifiers overlap between Data_Total and BV_fam.")
  Data_Total$Group <- rank_group[family_index]
  if (isTRUE(verbose)) {
    message("Global family groups constructed before individual-level eligibility filters.")
    message("Global families: ", n_family, "; eligible input trees: ", nrow(Data_Total))
    message("Positive / negative group allocation: ",
            positive_groups, " / ", negative_groups)
    message("Inflection candidates (positive): ",
            paste(pos_split$candidates, collapse = ", "))
    message("Inflection candidates (negative, local ranks): ",
            paste(neg_split$candidates, collapse = ", "))
    message("Selected boundaries: ", paste(cuts, collapse = ", "),
            "; methods: ", paste(cut_methods, collapse = "/"),
            "; rules: ", paste(finfl_Value, collapse = "/"))
    message("Families per group: ",
            paste(as.integer(table(factor(rank_group,
              levels = paste0("G", seq_len(nGroups))))), collapse = ", "))
  }
  boundary_colors <- ifelse(cut_methods == "zero_crossing", "gray60",
                             ifelse(Group_Boundaries$Side == "positive",
                                    "blue", "red"))
  centers <- (c(0L, cuts) + c(cuts, n_family) + 1) / 2
  families_rank <- ggplot2::ggplot(
    BV_fam, ggplot2::aes(x = Family, y = .data[[BV_Column]])) +
    ggplot2::geom_point(size = 2) +
    ggplot2::ggtitle(paste("Thinning Strategies based on", Trait, "Trait")) +
    ggplot2::ylab(paste0("Additive genetic effect (", BV_Column,
                             ") from ", Trait, " trait")) +
    ggplot2::xlab("Families") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                         colour = "gray60", linewidth = 0.5) +
    ggplot2::scale_x_discrete(
      guide = ggplot2::guide_axis(n.dodge = n.dodge_plot1,
                                  angle = angle_plot1))
  for (j in seq_along(cuts))
    families_rank <- families_rank +
      ggplot2::geom_vline(xintercept = cuts[j] + 0.5,
                          linetype = 4, color = boundary_colors[j],
                          linewidth = 0.7)
  for (j in seq_len(nGroups))
    families_rank <- families_rank +
      ggplot2::annotate("text", x = centers[j],
                        y = label.group.y[min(j, length(label.group.y))],
                        label = paste0("G", j))
  if (isTRUE(Plot.Rank)) plot(families_rank + additional_layer_plot1)
  if (isTRUE(STP)) {
    Data_Total$Plot <- as.character(Data_Total[[Bloc_Column]])
  } else {
    Data_Total$Plot <- paste0(Data_Total[[Family_Data_Total]], ":",
                               Data_Total[[Bloc_Column]])
  }

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
  Data_Total$Family <- trimws(as.character(Data_Total[[Family_Data_Total]]))
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
  if (length(all_scenarios) != 1L || is.na(all_scenarios) || !is.logical(all_scenarios))
    stop("all_scenarios must be TRUE or FALSE.")
  if (n_best_strategies < 1 || max_candidates < 1)
    stop("n_best_strategies and max_candidates must be positive.")
  n_best_strategies <- as.integer(n_best_strategies)
  max_candidates <- as.integer(max_candidates)
  groups <- paste0("G", seq_len(nGroups))
  normalize_minima <- function(x, default, argument) {
    if (is.null(x)) x <- default
    if (length(x) == 1L) x <- rep(x, nGroups)
    if (length(x) != nGroups || !is.numeric(x) ||
        any(!is.finite(x)) || any(x < 0) || any(x != floor(x)))
      stop(sprintf("%s must be a nonnegative integer scalar or a vector of length nGroups.", argument))
    as.integer(x)
  }
  min_group_families <- normalize_minima(
    min_group_families, c(1L, rep(0L, nGroups - 1L)), "min_group_families")
  min_group_trees <- normalize_minima(
    min_group_trees, rep(0L, nGroups), "min_group_trees")
  if (!is.logical(pareto) || length(pareto) != 1L || is.na(pareto))
    stop("pareto must be TRUE or FALSE.")
  if (!is.character(bw_selection_unit) ||
      length(bw_selection_unit) != 1L ||
      !bw_selection_unit %in% c("auto", "plot", "family"))
    stop("bw_selection_unit must be 'auto', 'plot', or 'family'.")
  make_cells <- function(unit) {
    z <- lapply(groups, function(g) {
      ix <- which(Data_Total$Group == g)
      if (!length(ix)) return(list())
      key <- if (unit == "family") Data_Total$Family[ix] else
        interaction(Data_Total$Family[ix], Data_Total$Plot[ix], drop = TRUE)
      partitions <- split(ix, key, drop = TRUE)
      lapply(partitions, function(ii) ii[order(
        -Data_Total[[blup_col]][ii],
        as.character(Data_Total[[id]][ii]))])
    })
    names(z) <- groups
    z
  }
  cells <- make_cells("plot")
  max_cell <- max(0L, unlist(lapply(cells, function(z)
    if (length(z)) max(lengths(z)) else 0L)))
  effective_bw_unit <- if (bw_selection_unit == "auto")
    if (max_cell <= 1L) "family" else "plot" else bw_selection_unit
  if (effective_bw_unit == "family") {
    cells <- make_cells("family")
    max_cell <- max(0L, unlist(lapply(cells, function(z)
      if (length(z)) max(lengths(z)) else 0L)))
    if (bw_selection_unit == "auto")
      message("BW auto mode: singleton family/plot cells detected; ",
              "using per-family BLUP ranking. Set bw_selection_unit = 'plot' ",
              "for strict per-plot selection.")
  }
  if (max_cell < 1L) stop("No usable BW selection cells.")
  levels_n <- if (is.null(seq_combinations)) {
    c(0L, seq.int(1L, max_cell, by = length_seq_combinations))
  } else c(0L, as.integer(seq_combinations))
  levels_n <- sort(unique(pmin(levels_n, max_cell)))
  if (anyNA(levels_n) || any(levels_n < 0L))
    stop("Invalid selection intensities.")
  if (max_cell == 1L)
    warning("Only 0/1 BW levels are possible. Try bw_selection_unit = 'family' ",
            "for finer intensities.")
  if (isTRUE(verbose))
    message("BW selection levels: ", paste(levels_n, collapse = ", "),
            "; maximum trees per family/plot cell: ", max_cell)
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
  # In exhaustive mode no Ne filter is applied, and all scenario summaries are kept.
  if (all_scenarios) {
    if (!is.null(target_NE)) message("all_scenarios = TRUE: target_NE filtering is disabled.")
    target_NE <- NULL
    if (isTRUE(rebuild_selected)) warning("all_scenarios = TRUE with rebuild_selected = TRUE may require considerable memory.")
  }
  # Cache the selected-tree counts and cumulative BLUP sums once.
  # Candidate Ne requires only per-family counts; full data frames are not created.
  family_levels <- sort(unique(Data_Total$Family))
  nf <- length(family_levels)
  fam_index <- match(Data_Total$Family, family_levels)
  blups <- Data_Total[[blup_col]]
  nlevels <- length(levels_n)
  group_count <- vector("list", nGroups)
  group_sum <- vector("list", nGroups)
  group_tree_count <- vector("list", nGroups)
  group_fam_count <- vector("list", nGroups)
  for (k in seq_len(nGroups)) {
    cnt <- matrix(0L, nrow = nf, ncol = nlevels)
    sm <- matrix(0, nrow = nf, ncol = nlevels)
    for (cell in cells[[k]]) {
      if (!length(cell)) next
      fam <- fam_index[cell[1L]]
      running <- c(0, cumsum(blups[cell]))
      take <- pmin(length(cell), levels_n)
      cnt[fam, ] <- cnt[fam, ] + take
      sm[fam, ] <- sm[fam, ] + running[take + 1L]
    }
    group_count[[k]] <- cnt
    group_sum[[k]] <- sm
    group_tree_count[[k]] <- colSums(cnt)
    group_fam_count[[k]] <- colSums(cnt > 0L)
  }
  # Preserve the previous count-based Ne estimator.
  calc_ne_counts <- function(counts) {
    counts <- counts[counts > 0L]
    N <- sum(counts)
    if (!N) return(NA_real_)
    nfamily <- length(counts)
    avg <- N / nfamily
    variance <- if (N > 1L) sum((counts - avg)^2) / (N - 1L) else 0
    round(4 * nfamily * avg / (avg + 3 + variance / avg), 1)
  }
  # The filter is based on Ne, so GS is computed only for eligible candidates.
  matched_target <- function(ne) {
    if (is.null(target_NE)) return(NA_integer_)
    eligible <- which(abs(ne - target_NE) <= target_NE_tol)
    if (!length(eligible)) return(0L)
    eligible[which.min(abs(ne - target_NE[eligible]))]
  }
  base_trait <- mean(Data_Total[[Trait]], na.rm = TRUE)
  new_record <- function(kind, key, combo, ne, ns, nf_selected, sum_blup,
                         matched = NA_integer_, group_trees, group_families) {
    mean_blup <- sum_blup / ns
    gs <- if (is.finite(base_trait) &&
              abs(base_trait) > .Machine$double.eps)
      mean_blup / base_trait else NA_real_
    data.frame(
      Candidate_ID = key, Strategy = kind, NE = ne,
      GS = round(gs, 6), GS_absolute = mean_blup,
      N_Progeny = ns, nFamily = nf_selected,
      as.list(stats::setNames(as.integer(group_trees), paste0("N_Trees_G", seq_len(nGroups)))),
      as.list(stats::setNames(as.integer(group_families), paste0("N_Families_G", seq_len(nGroups)))),
      as.list(stats::setNames(combo,
        paste0("nID_fam_plot_G", seq_len(nGroups)))),
      Target_NE = if (is.na(matched)) NA_real_ else target_NE[matched],
      NE_Distance = if (is.na(matched)) NA_real_ else abs(ne - target_NE[matched]),
      check.names = FALSE
    )
  }
  # Use deterministic exhaustive screening for small search spaces.
  if (all_scenarios && total > .Machine$integer.max)
    stop("Exhaustive search exceeds supported integer indexing; restrict seq_combinations or nGroups.")
  ncheck <- as.integer(if (all_scenarios) total else min(total, max_candidates))
  if (!all_scenarios && total > max_candidates) {
    warning("Bounded candidate search: global Pareto optimality is not guaranteed.")
    # Uniformly sample combination codes, preserving the global RNG state.
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    if (total <= .Machine$integer.max) {
      codes <- sample.int(as.integer(total), ncheck, replace = FALSE)
    } else {
      codes <- unique(floor(stats::runif(ncheck * 3L, 1, total + 1)))
      codes <- utils::head(codes, ncheck)
    }
    # Add a deterministic candidate with G1 represented to protect coverage.
    if (min_group_families[1L] > 0L || min_group_trees[1L] > 0L) {
      positive <- which(levels_n > 0L)
      if (length(positive)) {
        anchor <- as.double(positive[1L] - 1L)
        codes <- unique(c(anchor, codes))
        codes <- utils::head(codes, ncheck)
      }
    }
  } else {
    codes <- seq_len(ncheck)
  }
  rows <- vector("list", length(codes) +
                   length(if (is.null(IS)) c(seq(.01,.50,by=.01),seq(.55,1,by=.05)) else IS))
  nrows <- 0L
  bw_evaluated <- 0L
  bw_matched <- 0L
  for (code in codes) {
    combo <- decode(code)
    lev <- match(combo, levels_n)
    contributions <- integer(nf)
    for (k in seq_len(nGroups)) {
      contributions <- contributions + group_count[[k]][, lev[k]]
    }
    ns <- sum(contributions)
    if (!ns) next
    bw_evaluated <- bw_evaluated + 1L
    gtrees <- vapply(seq_len(nGroups), function(k)
      group_tree_count[[k]][lev[k]], numeric(1))
    gfams <- vapply(seq_len(nGroups), function(k)
      group_fam_count[[k]][lev[k]], numeric(1))
    if (any(gtrees < min_group_trees) || any(gfams < min_group_families))
      next
    ne <- calc_ne_counts(contributions)
    matched <- matched_target(ne)
    if (!is.null(target_NE) && matched == 0L) next
    bw_matched <- bw_matched + 1L
    sum_blup <- 0
    for (k in seq_len(nGroups)) {
      sum_blup <- sum_blup + sum(group_sum[[k]][, lev[k]])
    }
    nrows <- nrows + 1L
    rows[[nrows]] <- new_record(
      "BW", paste0("BW_", code), combo, ne, ns,
      sum(contributions > 0L), sum_blup, matched,
      group_trees = gtrees, group_families = gfams
    )
  }
  # BI selects the highest individual BLUPs across all families.
  is_vec <- if (is.null(IS)) c(seq(0.01, 0.50, by = 0.01), seq(0.55, 1, by = 0.05)) else IS
  if (!is.numeric(is_vec) || any(!is.finite(is_vec)) ||
      any(is_vec <= 0 | is_vec > 1)) {
    stop("IS must be selection proportions in (0, 1].")
  }
  order_bi <- order(-blups, as.character(Data_Total[[id]]))
  blup_cumulative <- c(0, cumsum(blups[order_bi]))
  for (i in seq_along(is_vec)) {
    ns <- max(1L, min(nrow(Data_Total), round(nrow(Data_Total) * is_vec[i])))
    counts <- tabulate(fam_index[order_bi[seq_len(ns)]], nbins = nf)
    gindex <- match(Data_Total$Group[order_bi[seq_len(ns)]], groups)
    gtrees <- tabulate(gindex, nbins = nGroups)
    gfams <- vapply(seq_len(nGroups), function(k) {
      ix <- order_bi[seq_len(ns)][gindex == k]
      length(unique(fam_index[ix]))
    }, integer(1))
    if (any(gtrees < min_group_trees) || any(gfams < min_group_families))
      next
    ne <- calc_ne_counts(counts)
    matched <- matched_target(ne)
    if (!is.null(target_NE) && matched == 0L) next
    nrows <- nrows + 1L
    rows[[nrows]] <- new_record(
      "BI", paste0("BI_", i), rep(NA_integer_, nGroups),
      ne, ns, sum(counts > 0L), blup_cumulative[ns + 1L], matched,
      group_trees = gtrees, group_families = gfams
    )
  }
  if (verbose) {
    message(sprintf(
      "Candidate screening: BW evaluated=%d; BW matched target Ne=%d; retained summaries=%d",
      bw_evaluated, bw_matched, nrows
    ))
  }
  if (!nrows) {
    warning("No candidate strategy matched the requested target Ne range.")
    screened <- data.frame()
  } else {
    screened <- dplyr::bind_rows(rows[seq_len(nrows)])
    screened <- screened[order(-screened$GS_absolute,
                               screened$NE_Distance,
                               screened$Candidate_ID), , drop = FALSE]
    # Identify Pareto-efficient candidates in O(n log n), maximizing both
    # individual BLUP gain and family-based effective population size.
    # Equal objective values are equivalent; keep ties on the same front.
    pareto_order <- order(-screened$GS_absolute, -screened$NE,
                          screened$Candidate_ID)
    front <- logical(nrow(screened))
    best_ne <- -Inf
    best_gs <- Inf
    for (j in pareto_order) {
      gs <- screened$GS_absolute[j]
      ne <- screened$NE[j]
      dominated <- ne < best_ne || (ne == best_ne && gs < best_gs)
      front[j] <- !dominated
      if (ne > best_ne) best_ne <- ne
      if (gs < best_gs) best_gs <- gs
    }
    screened$Pareto <- front
    screened <- screened[order(!screened$Pareto,
                                -screened$GS_absolute,
                                -screened$NE,
                                screened$Candidate_ID), , drop = FALSE]
    if (!all_scenarios) {
      # Preserve both BW and BI when feasible; fill remaining slots by
      # Pareto precedence, then absolute gain, within each Ne target.
      choose <- function(z) {
        if (!pareto) z <- z[order(-z$GS_absolute, -z$NE), , drop = FALSE]
        keep <- integer()
        for (kind in c("BW", "BI")) {
          ix <- which(z$Strategy == kind)
          if (length(ix)) keep <- c(keep, ix[1L])
        }
        keep <- unique(c(keep, seq_len(nrow(z))))
        z[utils::head(keep, max(n_best_strategies, min(2L, nrow(z)))),
          , drop = FALSE]
      }
      if (is.null(target_NE)) screened <- choose(screened)
      else screened <- dplyr::bind_rows(lapply(unique(screened$Target_NE),
        function(x) choose(screened[screened$Target_NE == x, , drop = FALSE])))
    }
  }
  if (!nrow(screened)) {
    warning("No strategies matched the current selection criteria.")
    return(list(Thinning = screened, Strategies = list(), BV_Fam = BV_fam,
                Group_Boundaries = Group_Boundaries,
                ggPlot_families_rank = families_rank))
  }
  screened$N_Strategy <- seq_len(nrow(screened))
  screened$Trait <- Trait
  screened[["GS(%)"]] <- scales::percent(screened$GS, accuracy = 0.1)
  screened[["N_Progeny(%)"]] <- scales::percent(
    screened$N_Progeny / nrow(Data_Total), accuracy = 0.1)
  Thinning <- screened
  Strategies <- list()
  # Always provide selected individual IDs. Full reconstruction is optional.
  # Compact records avoid copying the entire experimental data frame.
  for (i in seq_len(nrow(screened))) {
    candidate <- screened$Candidate_ID[i]
    ix <- if (screened$Strategy[i] == "BW") {
      select_bw(decode(as.double(sub("BW_", "", candidate))))
    } else {
      j <- as.integer(sub("BI_", "", candidate))
      n <- max(1L, min(nrow(Data_Total),
                     round(nrow(Data_Total) * is_vec[j])))
      order_bi[seq_len(n)]
    }
    if (isTRUE(rebuild_selected)) {
      item <- Data_Total[ix, , drop = FALSE]
      item$BV <- item[[blup_col]] + base_trait
    } else {
      item <- data.frame(
        TreeID = as.character(Data_Total[[id]][ix]),
        Family = as.character(Data_Total$Family[ix]),
        Group = as.character(Data_Total$Group[ix]),
        BLUP = as.numeric(Data_Total[[blup_col]][ix]),
        stringsAsFactors = FALSE)
    }
    item$Strategy <- screened$Strategy[i]
    item$N_Strategy <- i
    Strategies[[as.character(i)]] <- item
  }
  # Group assignment already follows the exact ranking-plot boundaries.
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
       Group_Boundaries = Group_Boundaries,
       ggPlot_families_rank = families_rank)
}
