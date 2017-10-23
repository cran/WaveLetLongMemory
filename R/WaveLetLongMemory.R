#=======================================================================#
# PURPOSE : Long Memory Parameter Estimation                            #
# AUTHOR  : Sandipan Samanta & Dr. Ranjit Kumar Paul                    #
# DATE    : 23 October, 2017                                            #
# VERSION : Ver 0.1.2                                                   #
#=======================================================================#

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
    Estimates <- as.data.frame(GPH.Estimation(Xt,bandwidth))
  }else if(Method=="SEMIPARAMETRIC"){

    #======================#
    #Semiparametric METHOD #
    #======================#

    SEM.Estimation = function(X,bandwidth,BetaLagParzen){
      SPERIO = fracdiff::fdSperio(X, bandw.exp = bandwidth, beta = BetaLagParzen)
      Semi.Para.Est = cbind(SEMEstimates=SPERIO$d,SPERIOStandardDev=SPERIO$sd.as,SPERIOStandardError=SPERIO$sd.reg)
      return(Semi.Para.Est)
    }
    Estimates <- as.data.frame(SEM.Estimation(Xt,bandwidth,BetaLagParzen))
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

    Estimates <- as.data.frame(Wavelet.Estimation(Xt,typeWvtrans,filtertype))
    colnames(Estimates) <- c('WaveletEstimates','WaveletStandardDev','WaveletStandardError')
  }else("Please choose the correct method !!!")
  return(Estimates)
}

