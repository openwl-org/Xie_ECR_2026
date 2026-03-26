# regression_helpers.R
# Shared regression functions from Xie et al. (2026)
# Extracted from chemcells_rnaseq_stat_1_regression.Rmd

suppressPackageStartupMessages({
  library(MASS)
  library(lme4)
  library(lmerTest)
  library(splines)
  library(tidyverse)
})

# --- Data preparation ---

#' Prepare gene data with lane averaging
#' Groups by SampleName and averages Y across technical replicate lanes
mod.dat <- function(Y, rna, X, pheno) {
  out <- pheno
  out$X <- as.numeric(pheno[, X])
  out$Y <- rna[Y, out$RNASeq_label]
  out <- out %>% group_by(SampleName) %>% mutate(Y = mean(Y)) %>% ungroup()
  out <- out[!duplicated(out$SampleName), ]
  return(data.frame(out))
}

#' Prepare gene data WITHOUT lane averaging (for mixed models)
mod.dat2 <- function(Y, rna, X, pheno) {
  out <- pheno
  out$X <- as.numeric(pheno[, X])
  out$Y <- rna[Y, out$RNASeq_label]
  return(data.frame(out))
}

# --- Linear models ---

#' OLS linear regression: log(Y) ~ X with lane averaging
mod.lm <- function(Y, rna, X, pheno) {
  dat.glm <- mod.dat(Y, rna, X, pheno)
  dat.glm$Y <- log(dat.glm$Y)
  coef.glm <- summary(lm(Y ~ X, dat.glm))$coef

  if ("X" %in% row.names(coef.glm)) {
    out <- data.frame(Gene = Y,
                      Estimate = coef.glm["X", "Estimate"],
                      SE = coef.glm["X", 2],
                      P = coef.glm["X", 4])
  } else {
    out <- data.frame(Gene = Y, Estimate = NA, SE = NA, P = NA)
  }
  return(out)
}

#' Mixed model linear: log(Y) ~ X + (1|SampleName) per-lane data
mod.lmer <- function(Y, rna, X, pheno) {
  dat.glm <- mod.dat2(Y, rna, X, pheno)
  dat.glm$Y <- log(dat.glm$Y)
  coef.glm <- summary(lmer(Y ~ X + (1 | SampleName), dat.glm))$coef

  if ("X" %in% row.names(coef.glm)) {
    out <- data.frame(Gene = Y,
                      Estimate = coef.glm["X", "Estimate"],
                      SE = coef.glm["X", 2],
                      P = coef.glm["X", 5])
  } else {
    out <- data.frame(Gene = Y, Estimate = NA, SE = NA, P = NA)
  }
  return(out)
}

# --- Quadratic models ---

#' OLS quadratic: log(Y) ~ X^2 + X with lane averaging
mod.lm2 <- function(Y, rna, X, pheno) {
  dat.glm <- mod.dat(Y, rna, X, pheno)
  dat.glm$X2 <- dat.glm$X^2
  dat.glm$Y <- log(dat.glm$Y)
  coef.glm <- summary(lm(Y ~ X2 + X, dat.glm))$coef

  if ("X2" %in% row.names(coef.glm)) {
    out <- data.frame(Gene = Y,
                      Estimate = coef.glm["X2", "Estimate"],
                      SE = coef.glm["X2", 2],
                      P = coef.glm["X2", 4])
  } else {
    out <- data.frame(Gene = Y, Estimate = NA, SE = NA, P = NA)
  }
  return(out)
}

#' Mixed model quadratic: log(Y) ~ X^2 + X + (1|SampleName)
mod.lmer2 <- function(Y, rna, X, pheno) {
  dat.glm <- mod.dat2(Y, rna, X, pheno)
  dat.glm$X2 <- dat.glm$X^2
  dat.glm$Y <- log(dat.glm$Y)
  coef.glm <- summary(lmer(Y ~ X2 + X + (1 | SampleName), dat.glm))$coef

  if ("X2" %in% row.names(coef.glm)) {
    out <- data.frame(Gene = Y,
                      Estimate = coef.glm["X2", "Estimate"],
                      SE = coef.glm["X2", 2],
                      P = coef.glm["X2", 5])
  } else {
    out <- data.frame(Gene = Y, Estimate = NA, SE = NA, P = NA)
  }
  return(out)
}

# --- Natural spline models ---

#' Mixed model spline: log(Y) ~ ns(X, df=n) + (1|SampleName)
mod.lmer3 <- function(Y, rna, X, pheno, n) {
  dat.glm <- mod.dat2(Y, rna, X, pheno)
  dat.glm$Y <- log(dat.glm$Y)
  coef.glm <- summary(lmer(Y ~ ns(X, df = n) + (1 | SampleName), dat.glm))$coef

  i <- grep("ns", rownames(coef.glm))

  out <- data.frame(Gene = Y,
                    ns = i - 1,
                    Estimate = coef.glm[i, "Estimate"],
                    SE = coef.glm[i, 2],
                    P = coef.glm[i, 5])
  rownames(out) <- NULL
  return(out)
}

# --- ANOVA model comparison ---

#' AIC comparison: linear vs quadratic vs ns(2) vs ns(3)
mod.anova <- function(Y, rna, X, pheno) {
  dat.glm <- mod.dat2(Y, rna, X, pheno)
  dat.glm$Y <- log(dat.glm$Y)
  dat.glm$X2 <- dat.glm$X^2
  mod1 <- lmer(Y ~ X + (1 | SampleName), dat.glm)
  mod2 <- lmer(Y ~ X + X2 + (1 | SampleName), dat.glm)
  mod3 <- lmer(Y ~ ns(X, df = 2) + (1 | SampleName), dat.glm)
  mod4 <- lmer(Y ~ ns(X, df = 3) + (1 | SampleName), dat.glm)

  out <- data.frame(Gene = Y,
                    AIC1 = anova(mod1, mod2)$AIC[1],
                    AIC2 = anova(mod1, mod2)$AIC[2],
                    AIC3 = anova(mod1, mod3)$AIC[2],
                    AIC4 = anova(mod1, mod4)$AIC[2],
                    P12 = anova(mod1, mod2)$`Pr(>Chisq)`[2],
                    P13 = anova(mod1, mod3)$`Pr(>Chisq)`[2],
                    P14 = anova(mod1, mod4)$`Pr(>Chisq)`[2])
  rownames(out) <- NULL
  return(out)
}

# --- Transcriptome-wide wrappers ---

#' Run OLS regression across all genes
TWAS.lm <- function(Yset, rna, X, pheno) {
  out <- data.frame(Gene = character(), Estimate = numeric(),
                    SE = numeric(), P = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.lm(Yset[i], rna, X, pheno))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}

#' Run mixed model regression across all genes
TWAS.lmer <- function(Yset, rna, X, pheno) {
  out <- data.frame(Gene = character(), Estimate = numeric(),
                    SE = numeric(), P = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.lmer(Yset[i], rna, X, pheno))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}

#' Run OLS quadratic regression across all genes
TWAS.lm2 <- function(Yset, rna, X, pheno) {
  out <- data.frame(Gene = character(), Estimate = numeric(),
                    SE = numeric(), P = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.lm2(Yset[i], rna, X, pheno))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}

#' Run mixed model quadratic regression across all genes
TWAS.lmer2 <- function(Yset, rna, X, pheno) {
  out <- data.frame(Gene = character(), Estimate = numeric(),
                    SE = numeric(), P = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.lmer2(Yset[i], rna, X, pheno))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}

#' Run mixed model spline regression across all genes
TWAS.lmer3 <- function(Yset, rna, X, pheno, n) {
  out <- data.frame(Gene = character(), ns = numeric(),
                    Estimate = numeric(), SE = numeric(), P = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.lmer3(Yset[i], rna, X, pheno, n))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}

#' Run ANOVA model comparison across all genes
TWAS.anova <- function(Yset, rna, X, pheno) {
  out <- data.frame(Gene = character(), AIC1 = numeric(), AIC2 = numeric(),
                    AIC3 = numeric(), AIC4 = numeric(), P12 = numeric(),
                    P13 = numeric(), P14 = numeric())
  n_iter <- length(Yset)
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
  for (i in 1:n_iter) {
    out <- rbind(out, mod.anova(Yset[i], rna, X, pheno))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(out)
}
