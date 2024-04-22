# This file allows you to source in the many functions we have implemented 
# in lab.

# Sourcing is functionally very similar to installing a package with a few less
# perks.

#### varScore ####
varScore = function (Data, Forward, Reverse = NULL, Range = NULL, Prorate = TRUE, 
                     MaxMiss = 0.2) 
{
  d = Data[, c(Forward, Reverse)]
  if (!is.null(Range)) {
    if (min(d, na.rm = TRUE) < Range[1] || max(d, na.rm = TRUE) > 
        Range[2]) {
      stop("Item score(s) out of range")
    }
  }
  if (!is.null(Reverse) && length(Range) != 2) {
    stop("Must specify item range (Range) to reverse score items")
  }
  if (!is.null(Reverse)) {
    for (v in Reverse) {
      d[, v] = (Range[1] + Range[2]) - d[, v]
    }
  }
  if (Prorate) {
    Total = rowMeans(d, na.rm = TRUE) * dim(d)[2]
  }
  else {
    Total = rowSums(d, na.rm = TRUE)
  }
  MissCount = rowSums(is.na(d))
  MissCount = MissCount/dim(d)[2]
  Total[MissCount > MaxMiss] = NA
  return(Total)
}


#### gplotPredict ####
ggplotPredict = function (Model, Data = NULL, Label = NULL, Type = "response") 
{
  if (is.null(Data) & class(Model)[1] == "lm") {
    return(fitted.values(Model))
  }
  else {
    if (is.null(Label)) {
      PredictName = "Predicted"
      CILoName = "CILo"
      CIHiName = "CIHi"
      SEName = "SE"
    }
    else {
      PredictName = paste0("Predicted", Label)
      CILoName = paste0("CILo", Label)
      CIHiName = paste0("CIHi", Label)
      SEName = paste0("SE", Label)
    }
    Predictions = matrix(data = NA, nrow = nrow(Data), ncol = 4, 
                         dimnames = list(1:nrow(Data), c(PredictName, CILoName, 
                                                         CIHiName, SEName)))
    if (class(Model)[1] == "lm") {
      CILevel = 1 - 2 * pt(c(1), df = Model$df.residual, 
                           lower.tail = FALSE)
      Predictions[, 1:3] = predict(Model, newdata = Data, 
                                   interval = "confidence", level = CILevel)
      Predictions[, 4] = Predictions[, 1] - Predictions[, 
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if (class(Model)[1] == "glm") {
      tmpPred = predict(Model, newdata = Data, type = "link", 
                        se.fit = TRUE)
      upr <- tmpPred$fit + tmpPred$se.fit
      lwr <- tmpPred$fit - tmpPred$se.fit
      fit <- tmpPred$fit
      if (Type == "response") {
        fit <- Model$family$linkinv(fit)
        upr <- Model$family$linkinv(upr)
        lwr <- Model$family$linkinv(lwr)
      }
      Predictions[, 1] = fit
      Predictions[, 2] = lwr
      Predictions[, 3] = upr
      Predictions[, 4] = Predictions[, 1] - Predictions[, 
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if ((class(Model)[1] == "lmerMod") || (class(Model)[1] == 
                                           "glmerMod")) {
      Predictions[, c(1, 4)] = predictSE(Model, Data, 
                                         se.fit = TRUE, type = Type, level = 0, print.matrix = TRUE)
      Predictions[, 2] = Predictions[, 1] - Predictions[, 
                                                        4]
      Predictions[, 3] = Predictions[, 1] + Predictions[, 
                                                        4]
    }
    if (any(names(Data) == PredictName) || any(names(Data) == 
                                               CILoName) || any(names(Data) == CIHiName) || any(names(Data) == 
                                                                                                SEName)) {
      warning("Variable names (Predicted, CILo, CIHi, SE with Label PostFix) used in Data.  These variables removed before merging in predicted values")
      Data[, c(PredictName, CILoName, CIHiName, SEName)] = list(NULL)
    }
    Data = data.frame(Predictions, Data)
    return(Data)
  }
}

#### modelCaseAnalysis ####
modelCaseAnalysis = function (Model, Type = "RESIDUALS", Term = NULL, ID = row.names(Model$model)) 
{
  switch(toupper(Type), UNIVARIATE = {
    d = Model$model
    {
      Vars = names(d)
    }
    for (varname in Vars) {
      par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2, ask = TRUE)
      if (is.factor(d[[varname]])) {
        plot(d[varname], xlab = varname, ylab = "Frequency")
      } else {
        hist(d[[varname]], xlab = varname, main = "Red: Mean +- 3SD; Green: Median +- 2.2IQR")
        text(d[[varname]], rep(0, length(d[[varname]])), labels = ID, pos = 3, cex = 0.7)
        abline(v = c(-3, 0, 3) * sd(d[[varname]]) + 
                 mean(d[[varname]]), col = "red", lty = c(1, 
                                                          2, 1))
        abline(v = c(-2.2, 0, 2.2) * IQR(d[[varname]]) + 
                 median(d[[varname]]), col = "green", lty = c(1, 
                                                              2, 1))
      }
    }
  }, HATVALUES = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Model: ", Model$call[2], "\n", "Small sample cut (green) = 3 * mean(Hat)\nLarge sample cut: 2 * mean(Hat)", 
                     sep = "")
    hist(hatvalues(Model), xlab = "Hat Values", main = TheTitle)
    abline(v = c(2, 3) * mean(hatvalues(Model)), col = c("red", 
                                                         "green"))
    text(hatvalues(Model), rep(0, length(hatvalues(Model))), labels = ID, pos = 3, cex = 0.7)
    points(hatvalues(Model), rep(0, length(hatvalues(Model))), 
           pch = "|", col = "blue")
  }, RESIDUALS = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    N = length(rstudent(Model))
    k = length(coef(Model)) - 1
    TCut <- qt(p = 0.025/N, df = N - k - 2, lower.tail = FALSE)
    TheTitle = paste("Model: ", Model$call[2], "\n", "Bonferroni corrected p < .05 cut-off in red", 
                     sep = "")
    hist(rstudent(Model), xlab = "Studentized Residuals", 
         main = TheTitle)
    abline(v = c(-1, 1) * TCut, col = "red")
    text(rstudent(Model), rep(0, length(rstudent(Model))), labels = ID, pos = 3, cex = 0.7)
    points(rstudent(Model), rep(0, length(rstudent(Model))), 
           pch = "|", col = "blue")
  }, COOKSD = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    N = length(cooks.distance(Model))
    k = length(coef(Model)) - 1
    TheTitle = paste("Model: ", Model$call[2], "\n", "4/(N-P) cut-off (red)\nqf(.5,P,N-P) cut-off (green)", 
                     sep = "")
    hist(cooks.distance(Model), xlab = "Cooks d", main = TheTitle)
    abline(v = c((4/(N - k - 1)), qf(0.5, k + 1, N - k - 
                                       1)), col = c("red", "green"))
    text(cooks.distance(Model), rep(0, length(cooks.distance(Model))), labels = ID, pos = 3, cex = 0.7)
    points(cooks.distance(Model), rep(0, length(cooks.distance(Model))), 
           pch = "|", col = "blue")
  }, DFBETAS = {
    if (is.null(Term)) {
      {
        Vars = dimnames(dfbetas(Model))[[2]]
      }
    } else {
      if (!(Term %in% dimnames(dfbetas(Model))[[2]])) {
        stop("Term specified for DFBETAS not valid")
      } else Vars = Term
    }
    for (varname in Vars) {
      par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2, ask = TRUE)
      TheTitle = paste("Model: ", Model$call[2], "\n", 
                       "B= ", coef(Model)[varname], sep = "")
      hist(dfbetas(Model)[, varname], xlab = paste("DFBETAS:", 
                                                   varname), main = TheTitle)
      text(dfbetas(Model)[, varname], rep(0, length(dfbetas(Model)[, 
                                                                   varname])), labels = ID, pos = 3, cex = 0.7)
      abline(v = c(-2, 2), col = "red")
      points(dfbetas(Model)[, varname], rep(0, length(dfbetas(Model)[, varname])), 
             pch = "|", col = "blue")
    }
    
  }, INFLUENCEPLOT = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Influence Bubble plot", "\nModel: ", 
                     Model$call[2], sep = "")
    plot(hatvalues(Model), rstudent(Model), type = "n", 
         xlab = "Hat Values", ylab = "Studentized Residuals", 
         main = TheTitle)
    cooksize = 10 * sqrt(cooks.distance(Model))/max(cooks.distance(Model))
    points(hatvalues(Model), rstudent(Model), cex = cooksize)
    N = length(rstudent(Model))
    k = length(coef(Model)) - 1
    TCut <- qt(p = 0.025/N, df = N - k - 2, lower.tail = FALSE)
    abline(h = c(-1, 0, 1) * TCut, col = "red", lty = c(1, 
                                                        2, 1))
    abline(v = c(1, 2, 3) * mean(hatvalues(Model)), col = "red", 
           lty = c(2, 1, 1))
    text(hatvalues(Model),x = hatvalues(Model), y = rstudent(Model),
         rep(0, length(hatvalues(Model))), labels = ID, pos = 3, cex = 0.7)
  }, COVRATIO = {
    N = length(covratio(Model))
    k = length(coef(Model)) - 1
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Model: ", Model$call[2], "\n", "abs((3*P)/N)-1 cut-off in red", 
                     sep = "")
    hist(covratio(Model), xlab = "CovRatio", main = TheTitle)
    abline(v = abs((3 * (k + 1)/N) - 1), col = "red")
    text(covratio(Model), rep(0, length(covratio(Model))), labels = ID, pos = 3, cex = 0.7)
  }, {
    print("Valid options for type: hatvalues, residuals, cooksd, dfbetas, influenceplot, covratio, univariate")
  })
  Rownames = row.names(Model$model)
  Cases = list(Rownames = Rownames, Values = ID)
  return(Cases)
}

#### modelAssumptions ####
modelAssumptions = function (Model, Type = "NORMAL", ID = row.names(Model$model), 
                             one.page = TRUE) 
{
  switch(toupper(Type), NORMAL = {
    if (one.page) {
      dev.new(width = 14, height = 7)
      par(mfrow = c(1, 2))
    } else {
      dev.new(width = 7, height = 7, record = TRUE)
    }
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    qqPlot(Model, labels = FALSE, sim = TRUE, main = "Quantile-Comparison Plot to Assess Normality", 
           xlab = "t Quantiles", ylab = "Studentized Residuals")
    plot(density(rstudent(Model)), main = "Density Plot to Assess Normality of Residuals", 
         xlab = "Studentized Residual")
    zx <- seq(-4, 4, length.out = 100)
    lines(zx, dnorm(zx, mean = 0, sd = sd(rstudent(Model))), 
          lty = 2, col = "blue")
    cat("Descriptive Statistics for Studentized Residuals\n")
    describe(rstudent(Model))
  }, CONSTANT = {
    if (one.page) {
      dev.new(width = 14, height = 7)
      par(mfrow = c(1, 2))
    } else {
      dev.new(width = 7, height = 7, record = TRUE)
    }
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    plot(rstudent(Model) ~ fitted.values(Model), main = "Studentized Residuals vs. Fitted Values", 
         xlab = "Fitted Values", ylab = "Studentized Residuals")
    abline(h = 0, lty = 2, col = "blue")
    print(spreadLevelPlot(Model))
    cat("\n\n")
    print(ncvTest(Model))
  }, LINEAR = {
    dev.new(width = 7, height = 7, record = TRUE)
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    crPlots(Model, ask = TRUE)
  }, {
    print("Valid options for type: normal, constant, linear")
  })
}

#### modelBoxCox ####
modelBoxCox = function (Model, Lambdas = seq(-2, 2, by = 0.1)) 
{
  LR <- boxCox(Model, lambda = Lambdas)
  Lambda1Index <- sum(LR$x < 1)
  Chi1 <- mean(c(LR$y[Lambda1Index], LR$y[Lambda1Index + 1]))
  ChiLambda <- LR$y[which.max(LR$y)]
  ChiDiff <- 2 * (ChiLambda - Chi1)
  print(paste("Best Lambda=", round(LR$x[which.max(LR$y)], 
                                    2)))
  print(paste("Chi-square (df=1)=", round(ChiDiff, 2)))
  print(paste("p-value=", round(pchisq(ChiDiff, df = 1, lower.tail = FALSE), 
                                5)))
}

#### mdes ####
#' Determine the minimum detectable effect size
#'
#' Given desired power, degrees of freedom, and alpha, this function returns the
#' minimum detectable effect size as an f-squared statistic.
#'
#' @param desired_power Defaults to .8, but this can be changed to test for the
#' minimum detectable effect size for other levels of power.
#' @param num_df Numerator degrees of freedom for the test.
#' @param den_df Denominator degrees of freedom for the test.
#' @param alpha Alpha level for the test.
#' @return Returns an f-squared statistic representing the minimum detectable effect
#' size with the desired power, degrees of freedom, and alpha.
#' @export
mdes <- function(desired_power = .8, num_df, den_df, alpha = .05) {
  
  # overshoot, nearest tenth
  power <- 0
  fsq <- 0
  while (power < desired_power) {
    fsq <- fsq + .01
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  # overshoot, nearest thousandth
  fsq <- fsq - .01
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .0001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  # overshoot, nearest ten thousandth
  fsq <- fsq - .0001
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .00001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  return(fsq)
  
}

#### power_analysis ####

#' Calculate power and for GLM tests
#'
#' Wrapper to calculate power for tests of parameter estimates or full model in
#' GLM based on Cohen's tables and using pwr.f2.test in pwr package. Inspired by
#' modelPower in John Curtin's lmSupport package. Allows the use
#' of partial eta squared or delta R2 rather than just f2 as effect size. If
#' you provide power, it returns N, if you provide N, it returns power. If you
#' provide N and power, it returns the minimum detectable effect size given the
#' specified N and power. You must alwasy specify an effect size as either
#' f2, partial eta2, or delta R2 with model R2. You must also specify the
#' number of parameters in the compact (pc) and  augmented (pa) for the model
#' comparison that will test the effect.
#'
#' @param pc Number of parameters in the compact model; i.e., intercept + all
#' parameters excluding the effect of interest; This is the numerator df of the
#' F test for the effect.
#' @param pa Number of parameters in the augmented model; i.e., the intercept
#' and all parameters including the effect of interest.
#' @param n Sample size.
#' @param alpha Alpha for statistical test.
#' @param power Power for statistical test.
#' @param f2 Effect size.
#' @param peta2 = Partial eta2 effect size.
#' @param dr2 Delta r2 effect; if provided must also specify r2.
#' @param r2 Model r2, only needed if using delta r2.
#' @return Returns a list with n, power, possibly minimum detectable effect size.
#' @examples
#' # return the minimum detectable effect size with 200 participants, power of
#' # .8 (the default), pa of 5, and pc of 4:
#' power_analysis(pa = 5, pc = 4, n = 200)
#'
#' # return the number of participants needed for 70% power given
#' # pa of 3, pc of 2, and peta2 of .01
#' power_analysis(pa = 3, pc = 2, peta2 = .01, power = .7)
#'
#' # return the power of a study with peta2 of .02 and 50 participants, with
#' # pa of 5 and pc of 4.
#' power_analysis(pa = 5, pc = 4, peta2 = .02, n = 50)
#' @export
power_analysis <- function(pc = NULL, pa = NULL, n = NULL, alpha = 0.05, power = NULL,
                           f2 = NULL, peta2 = NULL, dr2 = NULL, r2 = NULL) {
  if (is.null(pa) | is.null(pc)) {
    stop("Must provide pa and pc")
  }
  
  u <- pa - pc
  mdes_peta2 <- NULL
  mdes_f2 <- NULL
  
  nEffs <- 0
  if (!is.null(f2)) {
    nEffs <- nEffs + 1
  }
  if (!is.null(peta2)) {
    f2 <- peta2 / (1 - peta2)
    nEffs <- nEffs + 1
  }
  if (!is.null(dr2) & !is.null(r2)) {
    f2 <- dr2 / (1 - r2)
    nEffs <- nEffs + 1
  }
  if (nEffs > 1) {
    stop("Must not specify more than one of the following: f2, peta2, or both dr2 and r2")
  }
  
  if (!is.null(n)) {
    v <- n - pa
  }
  else {
    v <- NULL
  }
  
  if (!is.null(pa) & !is.null(pc) & !is.null(n)) {
    mdes_power <- power
    if(is.null(mdes_power)) mdes_power <- .8
    mdes_f2 <- mdes(desired_power = mdes_power, u, v, alpha)
    mdes_peta2 <- mdes_f2 / (1 + mdes_f2)
  }
  
  if (nEffs != 0) {
    results <- pwr::pwr.f2.test(
      u = u, v = v, f2 = f2, sig.level = alpha,
      power = power
    )
    
    output <- list(
      n = pa + results$v,
      peta2 = results$f2 / (1 + results$f2),
      f2 = results$f2,
      power = results$power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }
  
  if (nEffs == 0) {
    print("2")
    output <- list(
      n = n,
      peta2 = NULL,
      f2 = NULL,
      power = power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }
  
  return(output)
}