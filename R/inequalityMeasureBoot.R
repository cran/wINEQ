#' @title Weighted inequality measures
#'
#' @description Calculates weighted mean and sum of X (or median of X), and a set of relevant inequality measures.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param AF.norm (logical). If TRUE (default) then index is divided by its maximum possible value
#' @param Atkinson.e is a parameter for Atkinson coefficient
#' @param Jenkins.alfa is a parameter for Jenkins coefficient
#' @param Entropy.power is a generalized entropy index parameter
#' @param zeroes defines what to do with zeroes in the data vector. Possible options are "remove" and "include". See Entropy function for details.
#' @param Kolm.p is a parameter for Kolm index
#' @param Kolm.scale method of data standardization before computing
#' @param Leti.norm (logical). If TRUE (default) then Leti index is divided by a maximum possible value
#' @param AN_Y.a is a positive parameter for Abul Naga and Yalcin inequality measure
#' @param AN_Y.b is a parameter for Abul Naga and Yalcin inequality measure
#' @param Apouey.a is a parameter for Apouey inequality measure
#' @param Apouey.b is a parameter for Apouey inequality measure
#' @param BL.withsqrt if TRUE function returns index given by BL2, elsewhere by BL (default). See more in details of BL function.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr summarise
#'
#' @return The data frame with weighted mean and sum of X, and all inequality measures relevant for a numeric data.
#' In a case of an ordered factor, the data frame with median of X, and all relevant inequality measures.
#'
#' @rdname ineq_weighted
#'
#' @details Function checks if X is a numeric or an ordered factor. Then it calculates all appropriate inequality measures.
#'
#' @examples
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
ineq.weighted=function(
    X,
    W=rep(1,length(X)),
    AF.norm=TRUE,
    Atkinson.e=1,
    Jenkins.alfa=0.8,
    Entropy.power=0.5,
    zeroes='include',
    Kolm.p=1,
    Kolm.scale='Standardization',
    Leti.norm=T,
    AN_Y.a=1,
    AN_Y.b=1,
    Apouey.a=2/(1-length(W[!is.na(W) & !is.na(X)])),
    Apouey.b=length(W[!is.na(W) & !is.na(X)])/(length(W[!is.na(W) & !is.na(X)])-1),
    BL.withsqrt=FALSE
)
{
  ind=which(!is.na(W) & !is.na(X))
  if(length(ind)==0)return('Input with NAs only')
  W=W[ind];X=X[ind]
  if(!is.numeric(X) & !is.ordered(X))return('X must be numeric or ordered factor')

  n=length(X)
  data=data.frame(X,W)

  if(is.ordered(X))
  {
    result <- data %>%
      summarise(
        Median = medianf(X,W),
        Allison_Foster=AF(X,W,norm = AF.norm),
        Leti=Leti(X,W,norm = Leti.norm),
        Abul_Naga_Yalcin=AN_Y(X,W,a = AN_Y.a,b = AN_Y.b),
        Apouey=Apouey(X,W,a = Apouey.a,b = Apouey.b),
        Blair_Lacy=BL(X,W,withsqrt = BL.withsqrt)
      ) %>% as.data.frame()
  }
  if(is.numeric(X))
  {
    #if(min(X)==0){X=X+0.001*min(X[X>0])}
    result <- data %>%
      summarise(
        Mean = sum(W*X)/sum(W),
        Total = sum(W*X),
        Theil_L= Theil_L(X,W),
        Theil_T=Theil_T(X,W,zeroes),
        Hoover=Hoover(X,W),
        Gini=Gini(X,W),
        Atkinson=Atkinson(X,W,Atkinson.e),
        Kolm=Kolm(X,W,parameter = Kolm.p,scale = Kolm.scale),
        Entropy=Entropy(X,W,power = Entropy.power,zeroes),
        CoefVar=CoefVar(X,W),
        RicciSchutz=RicciSchutz(X,W),
        Prop20_20=Prop20_20(X,W),
        Palma=Palma(X,W),
        Jenkins=Jenkins(X,W,Jenkins.alfa)[,1],
        Cowell_and_Flachaire=Jenkins(X,W,Jenkins.alfa)[,2]
      ) %>% as.data.frame()
  }
  return(result)
}


#' @title Weighted inequality measures with bootstrap
#'
#' @description For weighted mean and weighted total of X (or median of X) as well as for each relevant inequality measure, returns outputs from ineq.weighted and bootstrap outcomes: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param B is a number of bootstrap samples.
#' @param AF.norm (logical). If TRUE (default) then index is divided by its maximum possible value
#' @param Atkinson.e is a parameter for Atkinson coefficient
#' @param Jenkins.alfa is a parameter for Jenkins coefficient
#' @param Entropy.power is a generalized entropy index parameter
#' @param zeroes defines what to do with zeroes in the data vector. Possible options are "remove" and "include". See Entropy function for details.
#' @param Kolm.p is a parameter for Kolm index
#' @param Kolm.scale method of data standardization before computing
#' @param Leti.norm (logical). If TRUE (default) then Leti index is divided by a maximum possible value
#' @param AN_Y.a is a positive parameter for Abul Naga and Yalcin inequality measure
#' @param AN_Y.b is a parameter for Abul Naga and Yalcin inequality measure
#' @param Apouey.a is a parameter for Apouey inequality measure
#' @param Apouey.b is a parameter for Apouey inequality measure
#' @param BL.withsqrt if TRUE function returns index given by BL2, elsewhere by BL (default). See more in details of BL function.
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
#' @return This functions returns a data frame from ineq.weighted extended with bootstrap results: expected value, bias (in %), standard deviation, coefficient of variation, lower and upper bound of confidence interval.
#' If keepSamples=TRUE or keepMeasures==TRUE then the output becomes a list. If keepSamples=TRUE, the functions returns  Xb and Wb, which are the samples of vector data and the samples of weights, respectively.
#' If keepMeasures==TRUE, the functions returns Mb, which is a set of inequality measures from bootstrapping.
#'
#' @rdname ineq_weighted_boot
#'
#' @details By default, naive bootstrap is performed, that is no weights calibration is conducted.
#' You can choose calibrated bootstrap to calibrate weights with respect to provided variables (Xs) and totals (total).
#' Confidence interval is simply derived with quantile of order \eqn{\alpha} and \eqn{1-\alpha} where \eqn{\alpha} is a significance level for confidence interval.
#'
#' @examples
#' # Inequality measures with additional statistics for numeric variable
#' X=1:10
#' W=1:10
#' ineq.weighted.boot(X,W,B=10)
#'
#' # Inequality measures with additional statistics for ordered factor variable
#' X=factor(c('H','H','M','M','L','L'),levels = c('L','M','H'),ordered = TRUE)
#' W=c(2,2,3,3,8,8)
#' ineq.weighted.boot(X,W,B=10)
#'
#' @export
ineq.weighted.boot=function(X,
                            W=rep(1,length(X)),
                            B=100,
                            AF.norm=TRUE,
                            Atkinson.e=1,
                            Jenkins.alfa=0.8,
                            Entropy.power=0.5,
                            zeroes='include',
                            Kolm.p=1,
                            Kolm.scale='Standardization',
                            Leti.norm=T,
                            AN_Y.a=1,
                            AN_Y.b=1,
                            Apouey.a=2/(1-length(W[!is.na(W) & !is.na(X)])),
                            Apouey.b=length(W[!is.na(W) & !is.na(X)])/(length(W[!is.na(W) & !is.na(X)])-1),
                            BL.withsqrt=FALSE,
                            keepSamples=FALSE,
                            keepMeasures=FALSE,
                            conf.alpha=0.05,
                            calib.boot=FALSE,
                            Xs=rep(1,length(X)),
                            total=sum(W),
                            calib.method='truncated',
                            bounds=c(low=0,upp=10)
)
{
  ind=which(!is.na(W) & !is.na(X))
  if(length(ind)==0)return('Input with NAs only')
  W=W[ind];X=X[ind]
  if(!is.numeric(X) & !is.ordered(X))return('X must be numeric or ordered factor')

  n=length(X)

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

  cols=ifelse(is.ordered(X),6,15)
  Mb=matrix(0,nrow=B,ncol=cols) %>% as.data.frame()

  medians=vector('character',B)
  for(i in 1:B)
  {
    if(is.ordered(X)){XX=factor(Xb[,i],levels = levels(X),ordered = T)}else{XX=Xb[,i]}
    m=unlist(ineq.weighted(XX,
                           Wb[,i],
                           AF.norm,
                           Atkinson.e,
                           Jenkins.alfa,
                           Entropy.power,
                           zeroes,
                           Kolm.p,
                           Kolm.scale,
                           Leti.norm,
                           AN_Y.a,
                           AN_Y.b,
                           Apouey.a,
                           Apouey.b,
                           BL.withsqrt))
    if(!is.ordered(X)){Mb[i,]=m}else{Mb[i,-1]=m[-1];medians[i]=as.character(medianf(XX,Wb[,i]))}
  }

  M=ineq.weighted(X,W,
                  AF.norm,
                  Atkinson.e,
                  Jenkins.alfa,
                  Entropy.power,
                  zeroes,
                  Kolm.p,
                  Kolm.scale,
                  Leti.norm,
                  AN_Y.a,
                  AN_Y.b,
                  Apouey.a,
                  Apouey.b,
                  BL.withsqrt)
  if(is.ordered(X)){medianM=M$Median;M$Median=NA}

  Mmean=colMeans(Mb)
  Msd=apply(Mb,2,sd)
  Mcv=Msd/Mmean*100
  Mq1=apply(Mb,2,quantile,probs=conf.alpha/2,na.rm=TRUE)
  Mq2=apply(Mb,2,quantile,probs=1-conf.alpha/2,na.rm=TRUE)
  Statistics.names=c('Index (I)','E(I)','Bias(I) [%]','Sd(I)','CV(I) [%]','Lower Q(I)','Higher Q(I)')
  Statistics.values=rbind(M,Mmean,(M/Mmean-1)*100,Msd,Mcv,Mq1,Mq2)
  Results=cbind(Statistics.names,Statistics.values)
  colnames(Mb)=colnames(M)
  if(is.ordered(X))
  {
    Mb=as.data.frame(Mb)
    Mb$Median=medians
    Results$Median=factor(NA,levels = levels(X),ordered = TRUE)
    Results$Median[1]=medianM
  }

  if(keepSamples==FALSE & keepMeasures==FALSE){return(Stats=Results)}
  if(keepSamples & keepMeasures){return(list(Stats=Results,Measures=Mb,Variables=Xb,Weights=Wb))}
  if(keepSamples){return(list(Stats=Results,Variables=Xb,Weights=Wb))}
  if(keepMeasures){return(list(Stats=Results,Measures=Mb))}
}
