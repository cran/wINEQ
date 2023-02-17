#' @title Weighted inequality measures
#'
#' @description Calculates weighted mean and sum of X, and a set of inequality measures.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param Atkinson.e is a parameter for Atkinson coefficient
#' @param Jenkins.alfa is a parameter for Jenkins coefficient
#' @param Entropy.e is a generalized entropy index parameter
#' @param Kolm.p is a parameter for Kolm index
#' @param Kolm.scale method of data standardization before computing
#' @param Leti.norm (logical). If TRUE (default) then Leti index is divided by a maximum possible value
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#'
#' @return The data frame with weighted mean and sum of X, and all inequality measures.
#'
#' @rdname ineq_weighted
#'
#' @details At this moment ineq.weighted calculates all inequality measures available in wINEQ packages.
#' In future, selection of inequality measures will be available
#'
#' @examples
#' library(dplyr)
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' ineq.weighted(X)
#' ineq.weighted(X,W)
#'
#'
#' data(Tourism)
#' # Results for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' ineq.weighted(X)
#' ineq.weighted(X,W)
#'
#' @export
ineq.weighted=function(X,W=rep(1,length(X)),Atkinson.e=1,Jenkins.alfa=0.8,Entropy.e=0.5,Kolm.p=1,Kolm.scale='Standardization',Leti.norm=T)
{
  n=length(X)
  data=data.frame(X,W)
  data=data[is.na(X)==F,]
  if((min(X,na.rm = T)==0)==T){X=X+0.001*min(X[X>0])}

  result <- data %>%
    summarise(
      Mean = sum(W*X)/sum(W),
      Total = sum(W*X),
      Theil_L= Theil_L(X,W),
      Theil_T=Theil_T(X,W),
      Hoover=Hoover(X,W),
      Gini=Gini(X,W),
      Atkinson=Atkinson(X,W,Atkinson.e),
      Kolm=Kolm(X,W,scale = Kolm.scale,parameter = Kolm.p),
      Entropy=Entropy(X,W,parameter = Entropy.e),
      CoefVar=CoefVar(X,W),
      RicciSchutz=RicciSchutz(X,W),
      Leti=Leti(X,W,norm = Leti.norm),
      Allison_Foster=AF(X,W),
      Prop20_20=Prop20_20(X,W),
      Palma=Palma(X,W),
      Jenkins=Jenkins(X,W,Jenkins.alfa)[,1],
      Cowell_and_Flachaire=Jenkins(X,W,Jenkins.alfa)[,2]
    ) %>% as.data.frame()
  return(result)
}


#' @title Weighted inequality measures with bootstrap
#'
#' @description For weighted mean and weighted total of X as well as for each inequality measure, returns outputs from ineq.weighted and bootstrap outcomes: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param B is a number of bootstrap samples.
#' @param Atkinson.e is a parameter for Atkinson coefficient
#' @param Jenkins.alfa is a parameter for Jenkins coefficient
#' @param Entropy.e is a generalized entropy index parameter
#' @param Kolm.p is a parameter for Kolm index
#' @param Kolm.scale method of data standardization before computing
#' @param Leti.norm (logical). If TRUE (default) then Leti index is divided by a maximum possible value
#' @param keepSamples if TRUE, it returns bootstrap samples of data (Xb) and weights (Wb)
#' @param keepMeasures if TRUE, it returns values of all inequality measures for each bootstrap sample
#' @param conf.alpha significance level for confidence interval
#' @param calib.boot if FALSE, then naive bootstrap is performed, calibrated bootstrap elsewhere
#' @param Xs matrix of calibration variables. By default it is a vector of 1's, applied if calib.boot is TRUE
#' @param total vector of population totals. By default it is a sum of weights, applied if calib.boot is TRUE
#' @param calib.method weights' calibration method for function calib (sampling)
#' @param bounds vector of bounds for the g-weights used in the truncated and logit methods; 'low' is the smallest value and 'upp' is the largest value
#'
#' @importFrom stats sd
#' @importFrom stats quantile
#' @importFrom sampling calib
#'
#' @return By default this functions returns a data frame from ineq.weighted for weighted mean and weighted total of X as well as for each inequality measure extended with bootstrap results: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.  If keepSamples=TRUE or keepMeasures==TRUE then the output becomes a list. If keepSamples=TRUE, the functions returns  Xb and Wb, which are the samples of vector data and the samples of weights, respectively. If keepMeasures==TRUE, the functions returns Mb, which is a set of inequality measures from bootstrapping.
#'
#' @rdname ineq_weighted_boot
#'
#' @details At this moment ineq.weighted.boot calculates all inequality measures available in wINEQ packages.
#' In future, selection of inequality measures will be available. By default, naive bootstrap is performed, that is no weights calibration is conducted.
#' You can choose calibrated bootstrap to calibrate weights with respect to proided variables (Xs) and totals (total).
#' Confidence interval is simply derived with quantile of order \eqn{\alpha} and \eqn{1-\alpha} where \eqn{\alpha} is a significance level for confidence interval.
#'
#' @examples
#' library(sampling);library(dplyr)
#' # Inequality measures with additional statistics
#' X=1:10
#' W=1:10
#' ineq.weighted.boot(X,W)
#'
#'
#'
#' @export
ineq.weighted.boot=function(X,W=rep(1,length(X)),B=100,
                            Atkinson.e=1,Jenkins.alfa=0.8,Entropy.e=0.5,Kolm.p=1,Kolm.scale='Standardization',Leti.norm=T,
                            keepSamples=FALSE,keepMeasures=FALSE,conf.alpha=0.05,
                            calib.boot=FALSE,
                            Xs=rep(1,length(X)),total=sum(W), calib.method='truncated',bounds=c(low=0,upp=10)
)
{
  Z=sample(1:length(X),size = B*length(X),replace = TRUE,prob = 1/W)
  Xb=matrix(X[Z],nrow=length(X),ncol=B)
  Wb=matrix(W[Z],nrow=length(W),ncol=B)
  if(calib.boot)
  {
    for(i in 1:B)
    {
      Wb[,i]=Wb[,i]*calib(Xs = Xs,d = Wb[,i],total = total,method = calib.method,bounds=bounds)
    }
  }

  Mb=matrix(0,nrow=B,ncol=17)
  for(i in 1:B)
  {
    Mb[i,]=unlist(ineq.weighted(Xb[,i],Wb[,i],Atkinson.e,Jenkins.alfa,Entropy.e,Kolm.p))
  }

  M=ineq.weighted(X,W,Atkinson.e,Jenkins.alfa,Entropy.e,Kolm.p)
  Mmean=colMeans(Mb)
  Msd=apply(Mb,2,sd)
  Mcv=Msd/Mmean*100
  Mq1=apply(Mb,2,quantile,probs=conf.alpha/2,na.rm=TRUE)
  Mq2=apply(Mb,2,quantile,probs=1-conf.alpha/2,na.rm=TRUE)

  Results=cbind(Statistics=c('Index (I)','E(I)','Bias(I) [%]','Sd(I)','CV(I) [%]','Lower Q(I)','Higher Q(I)'),rbind(M,Mmean,(M/Mmean-1)*100,Msd,Mcv,Mq1,Mq2))

  if(keepSamples==FALSE & keepMeasures==FALSE){return(Stats=Results)}
  if(keepSamples & keepMeasures){return(list(Stats=Results,Measures=Mb,Variables=Xb,Weights=Wb))}
  if(keepSamples){return(list(Stats=Results,Variables=Xb,Weights=Wb))}
  if(keepMeasures){return(list(Stats=Results,Measures=Mb))}
}
