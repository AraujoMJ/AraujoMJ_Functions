---
bibliography: references.bib
---

# README

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons License" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

This repository has been created to preserve some useful functions for use in plant breeding data analysis using the `{lmerTest}` and `{breedR}` R packages. Most of these functions have strong applicability in forest breeding programs, as the datasets used are in most cases from these programs. Please note that these functions are constantly being tested and may be subject to significant changes in the future. Please get in touch for more information.

# DiagFunc

This function has been designed for use in diagnostic analysis to:

-   Check the distribution assumptions of the data (currently normal distribution only )

    -   Check the skewness and kurtosis

    -   Perform the statistic normality test: Anderson-Darling, kolmogorov-Smirnov

-   Test of heterocedasticityt: Levene test

-   Identify discrepant data and outliers

-   Box-Cox test

-   Graphical analysis

# COMP_MODEL

This function allows us to compare different models fitted in the `remlf90` function of the [`{breedR}` R package](https://github.com/famuvie/breedR)[@breedR].

# Extract_h2a_breedR

The `Extract_h2a_breedR` has been developed to extract heritability and correlations from models into the `remlf90` function of the [`{breedR}` R package](https://github.com/famuvie/breedR).

# Deviance_BreedR

The `Deviance_BreedR` function is prompted to be used in the deviance test in models fitted using the `remlf90` function from the [`{breedR}` R package](https://github.com/famuvie/breedR).

# BV_BreedR

The `BV_BreedR` has been developed to extract the breeding values from models into the `remlf90` function of the [`{breedR}` R package](https://github.com/famuvie/breedR).

# ExtractLmer

The `ExtractLmer` has been developed to extract generic parameters from models into the `lmer` function of the `{lmerTest}` R package [@lmerTest].

# BV_lmer

The `BV_lmer` has been developed to extract the breeding values from models into the `lmer` function of the `{lmerTest}` R package.

# Thinning_BreedR

The `Thinning_BreedR` function has been developed for use in tree breeding strategies that involve thinning strategies and help to find a balance between breeding and conservation. It also allows us to test and compare different models included in the `{breedR}` package. More details on the use of this function can be found in the following papers:

1.  [Thinning strategies for *Eucalyptus dunnii* population: balance between breeding and conservation using spatial variation and competition model](https://link.springer.com/article/10.1007/s11295-021-01523-w)[@dearaujo2021].

2.  [Conservative or non-conservative strategy to advance breeding generation? A case study in *Eucalyptus benthamii* using spatial variation and competition model](https://sciendo.com/pdf/10.2478/sg-2023-0001) [@dearaujo2023].

# Reference
