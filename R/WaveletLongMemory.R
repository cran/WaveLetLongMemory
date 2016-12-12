#=======================================================================#
# PURPOSE : Long Memory Parameter Estimation                            #
# AUTHOR  : Sandipan Samanta & Dr. Ranjit Kumar Paul                    #
# DATE    : 11 Dec, 2016                                                #
# VERSION : Ver 0.1.0                                                     #
#=======================================================================#

rm(list=ls(all=TRUE))


WVLM<- function(Method,Xt,bandwidth,BetaLagParzen,typeWvtrans,filtertype)
{

  if(Method=="GPH"){
    #===========#
    #GPH METHOD #
    #===========#

    GPH.Estimation = function(X,bandwidth){
      GPH = fracdiff::fdGPH(X, bandw.exp = bandwidth)
      GPH.Est = cbind(GPHEstimates=GPH$d,GPHStandardDev=GPH$sd.as,GPHStandardError=GPH$sd.reg)
      return(GPH.Est)
    }
    Estimates <- GPH.Estimation(Xt,bandwidth)
  }else if(Method=="SEMIPARAMETRIC"){

    #======================#
    #Semiparametric METHOD #
    #======================#

    SEM.Estimation = function(X,bandwidth,BetaLagParzen){
      SPERIO = fracdiff::fdSperio(X, bandw.exp = bandwidth, beta = BetaLagParzen)
      Semi.Para.Est = cbind(SEMEstimates=SPERIO$d,GPHStandardDev=SPERIO$sd.as,GPHStandardError=SPERIO$sd.reg)
      return(Semi.Para.Est)
    }
    Estimates <- SEM.Estimation(Xt,bandwidth,BetaLagParzen)
  }else if(Method=="WAVELET"){
    #===============#
    #Wavelet METHOD #
    #===============#

    Wavelet.Estimation = function(X,typeWvtrans,filtertype){
      WavLet.Transform = wmtsa::wavVar(X, xform=typeWvtrans,wavelet=filtertype)

      WavLet.Var = matrix(t(summary(WavLet.Transform)$vmat),,2)
      WavLet.Var.Unbiased = WavLet.Var[,2]

      Max.j = nrow(as.matrix(WavLet.Var[,2]))
      Ind.Series = (rep.int(2,Max.j))**rep(1:Max.j,1)

      WavLet.Reg = lm(log(WavLet.Var.Unbiased) ~ log(Ind.Series))
      Est.Std = subset( as.data.frame(summary(WavLet.Reg)$coefficients),
                        rownames(summary(WavLet.Reg)$coefficients)==c("log(Ind.Series)"))
      Reg.Result = matrix(Est.Std,1)
      WavLet.D.Est = matrix(0.5*(Reg.Result[[1]] + 1))
      WavLet.D.Std = matrix(0.5*(Reg.Result[[2]]))

      WaveLet.Est = cbind(WaveletEstimates=WavLet.D.Est,WaveletStandardDev=WavLet.D.Std,WaveletStandardError=WavLet.D.Std/sqrt(length(X)))
      return(WaveLet.Est)
    }

    Estimates <- Wavelet.Estimation(Xt,typeWvtrans,filtertype)
  }else("Please choose the correct method !!!")
  return(Estimates)
}
