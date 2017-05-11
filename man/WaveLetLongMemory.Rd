\name{Wavelet}
\alias{WVLM}
\title{Estimating Long Memory using wavelets}
\description{
The function WVLM estimates the long memory prarameter using wavelets as well as using other two methods namely GPH and Semiparametric.
}

\usage{
WVLM(Method,Xt,bandwidth,BetaLagParzen,typeWvtrans,filtertype)
}
\arguments{
  \item{Method}{GPH, SEMIPARAMETRIC, WAVELET}
  \item{Xt}{univariate time series}
  \item{bandwidth}{The bandwidth used in the regression equation}
  \item{BetaLagParzen}{exponent of the bandwidth used in the lag Parzen window}

  \item{typeWvtrans}{type of wavelet transform i.e. dwt or modwt}

  \item{filtertype}{Either a wt.filter object, a character string indicating which   wavelet filter to use in the decomposition, or a numeric vector of wavelet coefficients}
  }
\value{
  \item{Method}{GPH, SEMIPARAMETRIC, WAVELET.}
  \item{xt}{univariate time series.}
  \item{bandwidth}{The bandwidth used in the regression equation.}
  \item{WVLM}{Out Approach.}
  \item{GPH.Estimation}{The GPH estimator is based on the regression equation using the          periodogram function as an estimate of the spectral density.}
  \item{SEM.Estimation}{It is based on the regression equation using the smoothed                periodogram function as an estimate of the spectral density..}
  \item{Wavelet.Estimation}{WAVELET method makes use Jensen (1994) estimator to estimate         the memory parameter d in the ARFIMA(p,d,q) model based on wavelet technique.}
}
\references{

Geweke, J. and Porter-Hudak, S. (1983) The estimation and application of long memory time series models. Journal of Time Series Analysis 4(4), 221 to 238.

Robinson, P. M. (1995). Gaussian Semiparametric Estimation of Long Range Dependence. The Annals of Statistics  23 (5), 1630 to 1661.

Jensen, M.J.(1999). Using wavelets to obtain a consistent ordinary least squares estimator of the long-memory parameter journal of forecasting, Journal of Forecasting 18, 17 to 32.

Paul, R. K., Samanta, S. and Gurung, B. (2015). Monte Carlo simulation for comparison of different estimators of long memory parameter: An application of ARFIMA model for forecasting commodity price. Model Assisted Statistics and Application, 10(2), 116 to 127.

Reisen, V. A. (1994) Estimation of the fractional difference parameter in the ARFIMA(p,d,q) model using the smoothed periodogram. Journal Time Series Analysis, 15(1), 335 to 350.
}
\author{Sandipan Samanta, Ranjit Kumar Paul}
\examples{
## Simulating Long Memory Series
N <- 1000
PHI <- 0.2
THETA <- 0.1
SD <- 1
M <- 0
D <- 0.2
Seed <- 123

set.seed(Seed)
Sim.Series <- fracdiff::fracdiff.sim(n = N, ar = c(PHI), ma = c(THETA),
d = D, rand.gen = rnorm, sd = SD, mu = M)

Yt <- as.ts(Sim.Series$series)

## GPH Estimation
WVLM(Method="GPH",Xt=Yt,bandwidth = 0.5)

## SEMIPARAMETRIC Estimation
WVLM(Method="SEMIPARAMETRIC",Xt=Yt,bandwidth = 0.5,BetaLagParzen = 0.2)

## WAVELET Estimation using different filtertype
WVLM(Method="WAVELET",Xt=Yt,typeWvtrans = "modwt",filtertype = "haar")
WVLM(Method="WAVELET",Xt=Yt,typeWvtrans = "modwt",filtertype = "d6")
WVLM(Method="WAVELET",Xt=Yt,typeWvtrans = "modwt",filtertype = "s8")
}
