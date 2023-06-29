#--------------------------- Function for extract estimates from lmer model ---------------#
ExtractLmer <- function(Model = Model,
                        Trat = "tratamento",
                        Rep = "bloco",
                        nRep = 6,
                        nPlant = 6,
                        Plot = TRUE,
                        random = NULL) {
  if (class(Model)[1] != "lmerModLmerTest") {
    stop("This is not a lmerModLmerTest adjusted model")
  }
  
  # Extract variance components
  VarComp <- as.data.frame(VarCorr(Model)) %>%
    dplyr::select(!c("var1", "var2")) %>%
    `colnames<-`(c("FV", "Variance", "Standard Deviation")) %>%
    `rownames<-`(., .$FV)
  
  Plot <- ifelse(Plot == FALSE, NA, Plot)
  
  EST <- tibble(
    Vg = VarComp[Trat, "Variance"],
    Vplot = ifelse(Plot == TRUE, VarComp[paste0(Trat, sep = ":", Rep), "Variance"], Plot),
    Vrandom = ifelse(is.null(random), NA, VarComp[random, "Variance"]),
    Vres = VarComp["Residual", "Variance"],
    Vpheno = sum(Vg, Vplot, Vrandom, Vres, na.rm = T),
    `h2g` = Vg / Vpheno,
    `h2m` = Vg / sum(Vg, (Vplot / nRep), (Vres / (nRep * nPlant)), na.rm = T),
    `c2p` = ifelse(is.na(Vplot), NA, Vplot / Vpheno)
  ) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`("Estimates") %>%
    rownames_to_column(var = "Components")
  return(EST)
}