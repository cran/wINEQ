#' @title ineq.weighted
#'
#' @description Calculates weighted mean and sum of X, and a set of inequality measures.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param Atkinson.e is a parameter for calculating the value of the Atkinson coefficient
#' @param Jenkins.alfa is the Jenkins coefficient parameter
#' @param Entropy.e is a entropy parameter
#' @param Kolm.p is a Kolm parameter
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#'
#' @return The data frame with weighted mean and sum of X, and all inequality measures.
#'
#' @rdname ineq_weighted
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' ineq.weighted(X,W)
#'
#' @export
ineq.weighted=function(X,W=rep(1,length(X)),Atkinson.e=1,Jenkins.alfa=0.8,Entropy.e=0.5,Kolm.p=1)
{
  n=length(X)
  data=data.frame(X,W)
  data=data[is.na(X)==FALSE,]
  if((min(X,na.rm = TRUE)==0)==TRUE){X=X+0.001*min(X[X>0])}

  result <- data %>%
    summarise(
      Mean = sum(W*X)/sum(W),
      Total = sum(W*X),
      Theil_L= Theil_L(X,W),
      Theil_T=Theil_T(X,W),
      Hoover=Hoover(X,W),
      Gini=Gini(X,W),
      Atkinson=Atkinson(X,W,Atkinson.e),
      Kolm=Kolm(X,W),
      Entropy=Entropy(X,W,parameter = Entropy.e),
      CoefVar=CoefVar(X,W),
      RicciSchutz=RicciSchutz(X,W),
      Leti=Leti(X,W),
      Allison_Foster=AF(X,W),
      Prop20_20=Prop20_20(X,W),
      Palma=Palma(X,W),
      Jenkins=Jenkins(X,W,Jenkins.alfa)[,1],
      Cowell_and_Flachaire=Jenkins(X,W,Jenkins.alfa)[,2]
    ) %>% as.data.frame()
  return(result)
}


#' @title ineq.weighted.boot
#'
#' @description For weighted mean and weighted total of X as well as for each inequality measure, returns outputs from ineq.weighted and bootstrap outcomes: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param B numer of bootstrap samples.
#' @param Atkinson.e is a parameter for calculating the value of the Atkinson coefficient
#' @param Jenkins.alfa is the Jenkins coefficient parameter
#' @param Entropy.e is a entropy parameter
#' @param Kolm.p is a Kolm parameter
#' @param keepSamples if TRUE, it returns bootstrap samples of data (Xb) and weights (Wb)
#' @param keepMeasures if TRUE, it returns values of all inequality measures for each bootstrap sample
#' @param conf.alpha significance level for confidence interval
#' @param calib.boot if FALSE, then naive bootstrap is performed, calibrated bootstrap elsewhere
#' @param Xs matrix of calibration variables
#' @param total vector of population totals
#' @param calib.method weights' calibration method for function calib (sampling)
#'
#' @importFrom stats sd
#' @importFrom stats quantile
#' @importFrom sampling calib
#'
#' @return By default this functions returns a data frame from ineq.weighted for weighted mean and weighted total of X as well as for each inequality measure extended with bootstrap results: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.  If keepSamples=TRUE or keepMeasures==TRUE then the output becomes a list. If keepSamples=TRUE, the functions returns  Xb and Wb, which are the samples of vector data and the samples of weights, respectively. If keepMeasures==TRUE, the functions returns Mb, which is a set of inequality measures from bootstrapping.
#'
#' @rdname ineq_weighted_boot
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' ineq.weighted.boot(X,W)
#'
#' @export
ineq.weighted.boot=function(X,W=rep(1,length(X)),B=10,
                            Atkinson.e=1,Jenkins.alfa=0.8,Entropy.e=0.5,Kolm.p=1,
                            keepSamples=FALSE,keepMeasures=FALSE,conf.alpha=0.05,
                            calib.boot=FALSE,
                            Xs=rep(1,length(X)),total=sum(W), calib.method='truncated'
)
{
  Z=sample(1:length(X),size = B*length(X),replace = TRUE,prob = 1/W)
  Xb=matrix(X[Z],nrow=length(X),ncol=B)
  Wb=matrix(W[Z],nrow=length(W),ncol=B)
  if(calib.boot)
  {
    for(i in seq_along(Wb))
    {
      Wb[,i]=Wb[,i]*calib(Xs = Xs,d = Wb[,i],total = total,method = calib.method)
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

  if(keepSamples==FALSE & keepMeasures==FALSE){return(Results)}
  if(keepSamples & keepMeasures){return(list(Results,Mb,Xb,Wb))}
  if(keepSamples){return(list(Results,Xb,Wb))}
  if(keepMeasures){return(list(Results,Mb))}
}
