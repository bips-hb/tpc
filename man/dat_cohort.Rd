\name{dat_cohort}

\alias{dat_cohort, true_cohort}

\title{Simulated cohort data}

\description{Simulated data of a European child-and-youth cohort study with
three waves (t0, t1 and t2). See Andrews et al. (2021) for more information
on how the data were generated.}

\usage{
data("dat_cohort")
Rgraphviz::plot(true_cohort)
}

\format{
  A data frame with 5000 observations on 34 variables (10 variables were
  measured at three timepoints each, denoted as "_t0", "_t1" and "_t2").
  \describe{
    \item{\code{sex} (factor)}{Sex. Factor variable with levels "male" and
    "female".}
    \item{\code{country} (factor)}{Country of residence. Factor variable with
    levels "ITA", "EST", "CYP", "BEL", "SWE", "GER", "HUN" and "ESP".}
    \item{\code{fto} (factor)}{Genotype of an FTO gene high-risk allele. Factor
    variable with levels "TT", "AT", "AA".}
    \item{\code{birth_weight} (numeric)}{Birth weight in grams.}
    \item{\code{age_tx} (numeric)}{Age in years.}
    \item{\code{bmi_tx} (numeric)}{Body mass index z-score adjusted for sex and
    age.}
    \item{\code{bodyfat_tx} (numeric)}{Per cent body fat.}
    \item{\code{education_tx} (factor)}{Educational level. Factor variable with
    levels "low education", "medium education" and "high education".}
    \item{\code{fiber_tx} (numeric)}{Fiber intake in log(mg/kcal)}
    \item{\code{media_devices_tx} (numeric)}{Number of audiovisual media in the
    child's bedroom.}
    \item{\code{media_time_tx} (numeric)}{Use of audiovisual media in
    log(h/week+1).}
    \item{\code{mvpa_tx} (numeric)}{Moderate to vigorous physical activity in
    sqrt(min/day).}
    \item{\code{sugar_tx} (numeric)}{Square root of sugar intake score.}
    \item{\code{wellbeing_tx} (numeric)}{Box-Cox-transformed well-being score.}
  }
}

\references{
Andrews RM, Foraita R, Witte J (2021). A practical guide to causal discovery
with cohort data. Working paper.
}

\examples{
data("dat_cohort")
head(dat_cohort)
Rgraphviz::plot(true_cohort)
}
