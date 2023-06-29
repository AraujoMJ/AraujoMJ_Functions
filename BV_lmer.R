#----------------- Function for extract BLUPs from lmer ------------#
BV_lmer <- function(Model = Model,
                    Trat = "tratamento",
                    FixEffect = "bloco",
                    Rank = TRUE) {
  if (class(Model)[1] != "lmerModLmerTest") {
    stop("This is not a lmerModLmerTest adjusted model")
  }
  
  # Estimated Marginal Mean
  mu <- emmeans::emmeans(Model, specs = FixEffect) |>
    data.frame() |>
    summarise(mu = mean(emmean)) |>
    unlist()
  
  # BLUPs
  blups <- ranef(Model)[[Trat]] |>
    `colnames<-`("a") |>
    rownames_to_column(Trat) |>
    mutate(`u+a` = mu + a)
  
  # Arrange?
  if (Rank == TRUE) {
    blups <- blups |>
      arrange(desc(a))
  } else {
    blups <- blups
  }
  return(blups)
}