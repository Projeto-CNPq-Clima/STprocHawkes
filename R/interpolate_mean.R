#' Interpolate Mean Accumulated Values for Spatio-Temporal Hawkes Process
#'
#' @description
#' This function interpolates the mean accumulated values for a spatio-temporal Hawkes process
#' using the posterior samples from the `STModelHawkesMCMC` function, spatial and temporal data,
#' and user-defined covariates.
#'
#' @param results The output from the `STModelHawkesMCMC` function, containing posterior samples of model parameters.
#' @param data A matrix of event times for the phenomenon of interest, where rows represent time points and columns represent locations.
#' @param Sites1 A matrix of locations used to estimate the parameters of the spatio-temporal Hawkes process.
#' @param Sites2 A vector representing the location where the mean accumulated value will be interpolated.
#' @param Xw A matrix of covariates for \(W\) at the locations in `Sites1`.
#' @param Xm A matrix of covariates for \(M\) at the locations in `Sites1`.
#' @param Xu A matrix of covariates for \(U\) at the locations in `Sites1`.
#' @param Xmr A matrix of covariates for \(M\) at the interpolation location (`Sites2`).
#' @param Xwr A matrix of covariates for \(W\) at the interpolation location (`Sites2`).
#' @param Xur A matrix of covariates for \(U\) at the interpolation location (`Sites2`).
#' @param R3 A scalar radius defining the neighborhood for interpolation.
#' @param tau A vector of times where the mean accumulated values will be interpolated.
#' @param d A positive real number used in the distance weighting function.
#'
#' @return A matrix (`MatMean`) of interpolated mean accumulated values, where rows correspond to posterior samples, and columns correspond to time points in `tau`.
#'
#' @details
#' This function combines spatial and temporal information to compute mean accumulated values at a specified location and time points. The computation is based on posterior samples of the Hawkes process parameters \eqn{\gamma}, \eqn{\alpha}, \eqn{\eta}, \eqn{\beta}, with spatial interpolation using a distance-based weighting function.
#'
#' @export
interpolate_mean <- function(results, data, Sites1, Sites2, Xw, Xm, Xu, Xmr, Xwr, Xur, R3, tau, d) {
  # Combine current sites and new location for interpolation
  Stotal <- rbind(Sites1, t(as.matrix(Sites2)))
  RIND <- nrow(Stotal)

  # Initialize vectors and matrices
  VecGama <- VecAlpha <- VecEta <- VecBeta <- NULL
  MatMean <- NULL

  # Compute auxiliary information
  m <- nrow(data)
  tempdados <- is.na(data)
  nj <- as.matrix(m - apply(tempdados, 2, sum))
  DistS <- as.matrix(dist(Stotal))[RIND, -RIND]
  amostra <- as.vector(which(DistS <= R3))
  P3 <- FunW(DistS[amostra], d)

  # Iterate over posterior samples
  for (i in 1:nrow(results$MBeta)) {
    # Process W
    SIGMAWtotal <- gSigma(results$Mbw[i], results$Mvw[i], Stotal)
    SIGMAWA12 <- t(as.matrix(SIGMAWtotal[1:ncol(data), (ncol(data) + 1)]))
    SIGMAWA1 <- SIGMAWtotal[1:ncol(data), 1:ncol(data)]
    A2estw <- t(as.matrix(Xwr)) %*% as.matrix(results$MPsiw[i, ]) +
      SIGMAWA12 %*% solve(SIGMAWA1) %*% (as.matrix(results$MW[i, ]) - Xw %*% as.matrix(results$MPsiw[i, ]))
    SIGMAA2estw <- SIGMAWtotal[(ncol(data) + 1), (ncol(data) + 1)] -
      SIGMAWA12 %*% solve(SIGMAWA1) %*% t(SIGMAWA12)
    WNO <- rnorm(1, A2estw, sd = sqrt(SIGMAA2estw))
    Gama <- exp(WNO)
    VecGama <- c(VecGama, Gama)

    # Process M
    SIGMAMtotal <- gSigma(results$Mbm[i], results$Mvm[i], Stotal)
    SIGMAMA12 <- t(as.matrix(SIGMAMtotal[1:ncol(data), (ncol(data) + 1)]))
    SIGMAMA1 <- SIGMAMtotal[1:ncol(data), 1:ncol(data)]
    A2estm <- t(as.matrix(Xmr)) %*% as.matrix(results$MPsim[i, ]) +
      SIGMAMA12 %*% solve(SIGMAMA1) %*% (as.matrix(results$MM[i, ]) - Xm %*% as.matrix(results$MPsim[i, ]))
    SIGMAA2estm <- SIGMAMtotal[(ncol(data) + 1), (ncol(data) + 1)] -
      SIGMAMA12 %*% solve(SIGMAMA1) %*% t(SIGMAMA12)
    MNO <- rnorm(1, A2estm, sd = sqrt(SIGMAA2estm))
    Alphaest <- exp(MNO)
    VecAlpha <- c(VecAlpha, Alphaest)

    # Process U
    SIGMAUtotal <- gSigma(results$Mbu[i], results$Mvu[i], Stotal)
    SIGMAUA12 <- t(as.matrix(SIGMAUtotal[1:ncol(data), (ncol(data) + 1)]))
    SIGMAUA1 <- SIGMAUtotal[1:ncol(data), 1:ncol(data)]
    A2estu <- t(as.matrix(Xur)) %*% as.matrix(results$MPsiu[i, ]) +
      SIGMAUA12 %*% solve(SIGMAUA1) %*% (as.matrix(results$MU[i, ]) - Xu %*% as.matrix(results$MPsiu[i, ]))
    SIGMAA2estu <- SIGMAUtotal[(ncol(data) + 1), (ncol(data) + 1)] -
      SIGMAUA12 %*% solve(SIGMAUA1) %*% t(SIGMAUA12)
    UNO <- rnorm(1, A2estu, sd = sqrt(SIGMAA2estu))
    Eta <- exp(UNO)
    VecEta <- c(VecEta, Eta)

    # Beta
    BetaNO <- mean(results$MBeta[i, ])
    VecBeta <- c(VecBeta, BetaNO)

    # Temporal information
    sind <- sample(amostra, 1, prob = P3)
    temposant <- data[1:nj[sind], sind]

    # Mean interpolation
    VecMean <- NULL
    for (j in 1:length(tau)) {
      VecMean <- c(VecMean, MeanFunction(Gama, Eta, Alphaest, BetaNO, tau[j], Tempos(tau[j], temposant)))
    }
    MatMean <- rbind(MatMean, t(as.matrix(VecMean)))
  }

  return(MatMean)
}
