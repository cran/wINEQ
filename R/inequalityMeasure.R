#' @title Allison and Foster index
#'
#' @description Computes Allison and Foster inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @importFrom stats aggregate
#'
#' @return The value of Allison and Foster coefficient.
#'
#' @rdname Allison_and_Foster
#'
#' @details Let \eqn{m} be the median category and \eqn{p_i} be a share of \eqn{i}-th category. The following index was proposed by Allison and Foster (2004):
#'
#' @details
#' \deqn{AF =  \frac{\sum_{i=m}^n c_{i} p_{i} }{\sum_{i=m}^n p_{i}} - \frac{\sum_{i=1}^{m-1} c_{i} p_{i}}{\sum_{i=1}^{m-1} p_{i}}}
#'
#' @references Atkinson A. B.: (1970) On the measurement of inequality
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Cowell F. A.: (2011) Measuring Inequality (Third Edition) ISBN 978-0-19-959403-0
#' @references Allison R. A., Foster J E.: (2004) Measuring health inequality using qualitative data, Journal of Health Economics
#' @references Kobus M., Miłoś P.: (2011) Inequality decomposition by population subgroups for ordinal data, Working Papers, No. 24
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' AF(X)
#' AF(X,W)
#'
#' data(Well_being)
#' # Allison and Foster index for health assessment with sample weights:
#' X=Well_being$V11
#' W=Well_being$Weight
#' AF(X,W)
#'
#'
#' @export
AF=function(X,W=rep(1,length(X)))
{
  W=W[is.na(W)==FALSE]
  W=aggregate(W,by=list(X),FUN=sum)
  Fx=(cumsum(W$x)/sum(W$x))
  SW=cumsum(W$x)
  min=min(X)
  max=max(X)
  a=which(Fx<0.5)
  if(length(a)!=0)
  {A=max(a)
  b=floor(sum(W$x)/2)-SW[max(a)]
  c=W$x[max(a)+1]-b }else{A=0
  b=floor(sum(W$x)/2)
  c=W$x[1]-b}
  AF=((sum(W$Group.1[-c(a,A+1)]*W$x[-c(a,A+1)])+W$Group.1[A+1]*c)/(sum(W$x[-c(a,A+1)])+c))-((sum(W$Group.1[a]*W$x[a])+W$Group.1[A+1]*b)/(sum(W$x[a])+b))
  AFNorm=AF/(max-min)
  return(AFNorm)
}




#' @title Atkinson index
#'
#' @description Computes Atkinson inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param e is a coefficient of aversion to inequality, by default 1
#'
#' @return The value of Atkinson coefficient.
#'
#' @rdname Atkinson
#'
#' @details Atkinson coefficient with respect to parameter \eqn{\epsilon} is given by
#' \deqn{1-\frac{1}{\mu}{(\frac{1}{n}\sum_{i=1}^{n} x_{i}^{1-\epsilon} )}^{\frac{1}{1-\epsilon}}}
#' for \eqn{\epsilon \neq 1} and
#' \deqn{1-\frac{1}{\mu}{(\prod_{i=1}^{n} x_i)}^{\frac{1}{n}}}
#' for \eqn{\epsilon=1}.
#'
#' @references Atkinson A. B.: (1970) On the measurement of inequality, Journal of Economic Theory
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Allison P.D.: (1978) Measures of Inequality, American Sociological Review
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Atkinson(X)
#' Atkinson(X,W)
#'
#' data(Tourism)
#' # Atkinson index for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Atkinson(X,W)
#'
#' @export
Atkinson=function(X,W=rep(1,length(X)),e=1)
{
  if(e==1){A=1-prod(X^(W/sum(W)))/ (sum(W*X)/sum(W))}else{
    A=1-((1/sum(W)*sum(W*(X^(1-e))))^(1/(1-e)))/(sum(W*X)/sum(W))}
  return(A)
}


#' @title Generalized entropy index
#'
#' @description Computes generalized entropy index of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param parameter is a entropy parameter
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of generalized entropy index
#'
#' @rdname Entropy
#'
#' @details Entropy coefficient with respect to parameter \eqn{\alpha} is equal to Theil_L(X,W)  whenever \eqn{\alpha=0},
#' is equal to Theil_T(X,W) whenever \eqn{\alpha=1}, and whenever \eqn{\alpha \in (0,1)} we have
#' \deqn{GE(\alpha) = \frac{1}{\alpha(\alpha-1)W}\sum_{i=1}^{n}w_{i}\left(\left(\frac{x_{i}}{\mu}\right)^\alpha-1\right)}
#' where \eqn{W} is a sum of weights and \eqn{\mu} is the arithmetic mean of \eqn{x_{1},...,x_{n}}.
#'
#' @references Shorrocks A. F.: (1980) The Class of Additively Decomposable Inequality Measures. Econometrica
#' @references Pielou E.C.: (1966) The measurement of diversity in different types of biological collections. Journal of Theoretical Biology
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Entropy(X)
#' Entropy(X,W)
#'
#' data(Tourism)
#' # Generalized entropy index for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Entropy(X,W)
#'
#'
#' @export
Entropy=function (X,W=rep(1,length(X)), parameter = 0.5, na.rm = TRUE)
{
  if (!na.rm && any(is.na(X)))
    return(NA_real_)
  X <- as.numeric(na.omit(X))
  X <- as.numeric(X)
  if (is.null(parameter))
    parameter <- 0.5
  if (parameter == 0)
    e <- Theil_L(X,W)
  else if (parameter == 1)
    e <- Theil_T(X,W)
  else {
    k <- parameter
    e <- W*(X/(sum(W*X)/sum(W)))^k
    e <- sum(e - 1)/(k * (k - 1))/sum(W)
  }
  e
}


#' @title Kolm index
#'
#' @description Computes Kolm inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param parameter is a Kolm parameter
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#' @param scale method of data scaling (None, Normalization, Unitarization, Standardization)
#'
#' @importFrom stats na.omit
#'
#' @return The value of Kolm coefficient.
#'
#' @rdname Kolm
#'
#' @details Kolm index with parameter \eqn{\alpha}  is defined as:
#'          \deqn{K = \frac{1}{ \alpha} (log( \sum_{i=1}^n \exp(\alpha (w_{i} - \mu)) - log(n)))}
#'
#'  Kolm index is scale-dependent. Basic normalization methods can be applied before final computation.
#'
#' @references Kolm S. C.: (1976) Unequal inequalities I and II
#' @references Kolm S. C.: (1996) Intermediate measures of inequality
#' @references Chakravarty S. R.: (2009) Inequality, Polarization and Poverty  e-ISBN 978-0-387-79253-8
#'
#'
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Kolm(X)
#' Kolm(X,W)
#'
#' # Compare raw and standardized data.
#' Kolm(X,W)
#' Kolm(X,W, scale ="Standardization")
#'
#' # Changing units has an impact on the final result
#' Kolm(X)
#' Kolm(10*X)
#'
#' # Changing units has no impact on the final result with standardized data
#' Kolm(X,scale ="Standardization")
#' Kolm(10*X,scale ="Standardization")
#'
#' @export
Kolm=function (X,W=rep(1,length(X)), parameter = 1, na.rm = TRUE, scale = "None")
{
  if (!na.rm && any(is.na(X)))
    return(NA_real_)
  X <- as.numeric(na.omit(X))
  X <- as.numeric(X)

  if(scale=="Standardization"){X=(X-mean(X))/sd(X)}
  if(scale=="Unitarization"){X=(X-min(X))/(max(X)-min(X))}
  if(scale=="Normalization"){X=X/sqrt(sum(X^2))}

  if (is.null(parameter))
    parameter <- 1
  W <- W/sum(W)
  KM <- parameter * (sum(W*X) - X)
  KM <- sum(W*(exp(KM)))
  KM <- (1/parameter) * log(KM)
  KM
}


#' @title Ricci and Schutz index
#'
#' @description Computes Ricci and Schutz inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of Ricci and Schutz coefficient.
#'
#' @rdname RicciSchutz
#'
#' @details In the case of an empirical distribution with n elements where y_{i} denotes the wealth of household i and \eqn{\overline{y}} the sample average, the Ricci and Schutz coefficient can be expressed as:
#'          \deqn{RS =  \frac{1}{2n} \sum_{i=1}^{n} \frac{\mid y_{i} - \overline{y} \mid}{\overline{y}}}
#'
#'
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Eliazar I. I., Sokolov I. M.: (2010) Measuring statistical heterogeneity: The Pietra index
#' @references Costa R. N., Pérez-Duarte S.: (2019) Not all inequality measures were created equal, Statistics Paper Series, No 31
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' RicciSchutz(X)
#' RicciSchutz(X,W)
#'
#' data(Tourism)
#' #Ricci and Schutz index for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' RicciSchutz(X,W)
#'
#'
#' @export
RicciSchutz=function (X,W=rep(1,length(X)), na.rm = TRUE)
{
  if (!na.rm && any(is.na(X)))
    return(NA_real_)
  X <- as.numeric(na.omit(X))
  d <- abs(X - (sum(W*X)/sum(W)))
  d <- (sum(W*d)/sum(W))/(2 * (sum(W*X)/sum(W)))
  d
}




#' @title Coefficient of Variation
#'
#' @description Computes Coefficient of Variation inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param square logical, argument of the function CoefVar, for details see below
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of CoefVar coefficient.
#'
#' @rdname CoefVar
#'
#' @details Coefficient of variation is given by:
#' \deqn{CV= \frac{\sigma}{\mu}\times 100}
#' where \eqn{\sigma} is a standard deviation and \eqn{\mu} is arithmetic mean.
#'
#'
#' @references Sheret M.: (1984) Social Indicators Research, An International and Interdisciplinary Journal for Quality-of-Life Measurement, Vol. 15, No. 3, Oct. ISSN 03038300
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' CoefVar(X)
#' CoefVar(X,W)
#'
#' data(Tourism)
#' #Coefficient of variation for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' CoefVar(X,W)
#'
#'
#' @export
CoefVar=function (X,W=rep(1,length(X)), square = FALSE, na.rm = TRUE)
{
  if (!na.rm && any(is.na(X)))
    return(NA_real_)
  x <- as.numeric(na.omit(X))
  w.m <- sum(W*X)/sum(W)
  V <- sqrt(sum(W*(X-w.m)^2)/sum(W))/w.m
  if (square)
    V <- V^2
  V
}


#' @title Gini coefficient
#'
#' @description Computes Gini coefficient of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of Gini coefficient.
#'
#' @rdname Gini
#'
#' @details Gini coefficient is given by:
#'          \deqn{G = \frac{ \sum_{i=1}^n \sum_{j=1}^n \mid x_{i} - x_{j} \mid}{2n^{2} \overline{x}}}
#'
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Dixon P. M., Weiner, J., Mitchell-Olds, T., and Woodley, R.: (1987) Bootstrapping the Gini Coefficient of Inequality. Ecology , Volume 68 (5)
#' @references Firebaugh G.: (1999) Empirics of World Income Inequality, American Journal of Sociology
#' @references Deininger K.; Squire L.: (1996) A New Data Set Measuring Income Inequality, The World Bank Economic Review, Vol. 10, No. 3
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Gini(X)
#' Gini(X,W)
#'
#' data(Tourism)
#' #Gini coefficient for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Gini(X,W)
#'
#'
#' @export
Gini=function(X,W=rep(1,length(X))){sum(abs(matrix(X,length(X),length(X),TRUE)-matrix(X,length(X),length(X),FALSE))*matrix(W,length(W),length(W),TRUE)*matrix(W,length(W),length(W),FALSE)/(2*sum(W)^2*(sum(W*X)/sum(W))))}



#' @title Hoover index
#'
#' @description Computes Hoover inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of Hoover coefficient.
#'
#' @rdname Hoover
#'
#' @details Let x_{i} be the income of the i-th person and \eqn{\overline{x}} be the mean income. Then the Hoover index H is:
#' \deqn{H={\frac {1}{2}}{\frac {\sum_{i}|x_{i}-{\overline{x}}|}{\sum_{i}x_{i}}}}
#'
#'
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Hoover E. M. Jr.: (1936) The Measurement of Industrial Localization, The Review of Economics and Statistics, 18
#' @references Hoover E. M. Jr.: (1984) An Introduction to Regional Economics, ISBN 0-07-554440-7
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Hoover(X)
#' Hoover(X,W)
#'
#' data(Tourism)
#' #Hoover index for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Hoover(X,W)
#'
#'
#' @export
Hoover=function(X,W=rep(1,length(X))){(1/2)*(sum(W*abs(X - (sum(W*X)/sum(W))))/sum(X*W))}


#' @title Leti index
#'
#' @description Computes Leti inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector (ordered factor or numeric)
#' @param W is a vector of weights
#' @param norm (logical). If TRUE (default) then Leti index is divided by a maximum possible value which is \eqn{(k-1)/2} where \eqn{k} in a number of categories.
#'
#' @importFrom stats aggregate
#'
#' @return The value of Leti coefficient.
#'
#' @rdname Leti
#'
#' @details Let \eqn{n_{i}} be the number of individuals in category \eqn{i} and let \eqn{N} be the total sample size.
#' Cumulative distribution is given by \eqn{F_{i} = \frac{\sum_{j=1}^{i} n_{j}}{N}}. Leti index is defined as:
#'          \deqn{L =2 \sum_{i=1}^{k-1} F_{i}(1-F_{i})}
#'
#' @references Leti G.: (1983). Statistica descrittiva, il Mulino, Bologna. ISBN: 8-8150-0278-2
#'
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Leti(X)
#' Leti(X,W)
#'
#' data(Tourism)
#' #Leti index for Total expenditure with sample weights:
#' X= Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Leti(X,W)
#'
#'
#' @export
Leti=function(X,W=rep(1,length(X)),norm=T)
{
  W=W[is.na(W)==FALSE]
  W=aggregate(W,by=list(sort(X,decreasing = F)),FUN=sum)
  W=W$x
  Fx=(cumsum(W)/sum(W))
  Leti=2*sum(Fx*(1-Fx))
  if(norm){Leti=Leti*2/(length(Fx)-1)}
  return(Leti)
}


#' @title Jenkins, Cowell and Flachaire
#'
#' @description Computes Jenkins as well as Cowell and Flachaire inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param alfa is the Jenkins coefficient parameter
#'
#' @importFrom stats aggregate
#'
#' @return The value of Jenkins, Cowell and Flachaire coefficient.
#'
#' @rdname Jenkins
#'
#' @details Jenkins coefficient is given by:
#' \deqn{J=1-\sum_{j=0}^{K-1} (p_{j+1}-p_{j})(GL_{j}+GL_{j+1})}
#'
#' @details where GL is Generalized Lorenz curve.
#'
#' @details Cowell and Flachaire coefficient with alpha parameter is given by:
#' \deqn{I(\alpha)=\frac{1}{\alpha(\alpha-1)}(\frac{1}{N}\sum_{i=1}^{N}s_{i}^{\alpha}-1)}
#'
#' @details for \eqn{\alpha \in (0,1)}, and
#' \deqn{I(0)=-\frac{1}{N}\sum_{i=1}^{N} log(s_{i})}
#'
#' @details for \eqn{\alpha = 0}.
#'
#'
#'
#' @references Coulter P. B.: (1989) Measuring Inequality ISBN 0-8133-7726-9
#' @references Jenkins S. P. and P. J. Lambert: (1997) Three ‘I’s of Poverty Curves, with an Analysis of U.K. Poverty Trends
#' @references Cowell F. A.: (2000) Measurement of Inequality, Handbook of Income Distribution
#' @references Cowell F. A., Flachaire E.: (2017) Inequality with Ordinal Data
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Jenkins(X)
#' Jenkins(X,W)
#'
#' data(Tourism)
#' #Jenkins, Cowell and Flachaire coefficients for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Jenkins(X,W)
#'
#'
#' @export
Jenkins=function(X,W=rep(1,length(X)), alfa=0.8) #dane to wektor (konkretne pytanie), 0<= alfa < 1
{
  dane=as.vector(X)
  tab=aggregate(W,by=list(X),FUN=sum)

  yy=as.data.frame(tab)
  yy$fk=yy$x/sum(yy$x) #odsetek
  yy$Fk=cumsum(yy$x)/sum(yy$x) #dystrybuanta

  dejta=as.data.frame(tab$Group.1)
  dejta = merge(dejta,yy[,c("Group.1","fk","Fk")],by.y="Group.1",by.x = "tab$Group.1",all.x = TRUE)


  if(alfa!=0)
  {Cowell_and_Flachaire=1/(alfa*(alfa-1))*(1/nrow(dejta)*sum((dejta$Fk)^alfa)-1)
  }else{
    Cowell_and_Flachaire=-1/nrow(dejta)*sum(log(dejta$Fk))}

  P_I=seq(1:nrow(dejta))/nrow(dejta)

  GL=1/nrow(dejta)*cumsum(dejta$Fk)

  dejta$P_I=P_I
  dejta$GL=GL

  GL_J=aggregate(dejta$GL,by=list(dejta$`tab$Group.1`),FUN=max)
  P_I_J=aggregate(dejta$P_I,by=list(dejta$`tab$Group.1`),FUN=max)

  JJ=rbind(c(0,0,0),cbind(GL_J,P_I_J$x))
  colnames(JJ)=c("Group","GL_J","PI_J")

  J=1-sum((JJ$PI_J[2:nrow(JJ)]-JJ$PI_J[1:(nrow(JJ)-1)])*(JJ$GL_J[1:(nrow(JJ)-1)]+JJ$GL_J[2:nrow(JJ)]))

  wynik=matrix(c(J,Cowell_and_Flachaire),1,2)
  colnames(wynik)=c("Jenkins","Cowell_and_Flachaire")

  return(wynik)
}




#' @title Palma index
#'
#' @description Palma proportion - originally the ratio of the total income of the 10% richest people to the 40% poorest people.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @importFrom stats aggregate
#'
#' @return The value of Palma coefficient.
#'
#' @rdname Palma
#'
#' @details Palma index is calculated by the following formula:
#'          \deqn{Palma =\frac{H}{L}}
#'          where \eqn{H} is share of 10% of the highest values,
#'          \eqn{L} is share of 40% of the lowest values.
#'
#' @references Cobham A., Sumner A.: (2013) Putting the Gini Back in the Bottle? 'The Palma' as a Policy-Relevant Measure of Inequality
#' @references Palma J. G.: (2011) Homogeneous middles vs. heterogeneous tails, and the end of the ‘Inverted-U’: the share of the rich is what it’s all about
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Palma(X)
#' Palma(X,W)
#'
#' data(Tourism)
#' #Palma index for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Palma(X,W)
#'
#'
#'
#' @export
Palma=function(X,W=rep(1,length(X)))
{
  W=W[is.na(W)==F]
  W=aggregate(W,by=list(X),FUN=sum)
  Fx=(cumsum(W$x)/sum(W$x))
  SW=cumsum(W$x)
  a=which(Fx<0.4)
  if(length(a)!=0){A=max(a)
  b=floor(sum(W$x)*0.4)-SW[max(a)]}else{A=0
  b=floor(sum(W$x)*0.4)}
  c=which(Fx<0.9)
  if(length(c)!=0){C=max(c)
  d=W$x[max(c)+1]-(floor(sum(W$x)*0.9)-SW[max(c)])}else{C=0
  d=W$x[1]-(floor(sum(W$x)*0.9))}
  Palma=(sum(W$Group.1[-c(c,C+1)]*W$x[-c(c,C+1)])+W$Group.1[C+1]*d)/(sum(W$Group.1[a]*W$x[a])+W$Group.1[A+1]*b)

  return(Palma)
}

#' @title Proportion 20:20
#'
#' @description 20:20 ratio - originally the ratio of the total income of the 20% richest people to the 20% poorest people.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @importFrom stats aggregate
#'
#' @return The value of 20:20 ratio coefficient.
#'
#' @rdname Prop20_20
#'
#' @details 20:20 ratio is calculated as follows:
#'          \deqn{Prop =\frac{H}{L}}
#'          where \eqn{H} is share of 20% of the highest values,
#'          \eqn{L} is share of 20% of the lowest values.
#'
#'
#' @references Panel Data Econometrics: Theoretical Contributions And Empirical Applications edited by Badi Hani Baltag
#' @references Notes on Statistical Sources and Methods - The Equality Trust.
#'
#'
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Prop20_20(X)
#' Prop20_20(X,W)
#'
#' data(Tourism)
#' #Prop20_20 proportion for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Prop20_20(X,W)
#'
#'
#' @export
Prop20_20=function(X,W=rep(1,length(X)))
{
  W=W[is.na(W)==FALSE]
  W=aggregate(W,by=list(X),FUN=sum)
  Fx=(cumsum(W$x)/sum(W$x))
  SW=cumsum(W$x)
  a=which(Fx<0.2)
  a=which(Fx<0.2)
  if(length(a)!=0){A=max(a)
  b=floor(sum(W$x)*0.2)-SW[max(a)]}else{A=0
  b=floor(sum(W$x)*0.2)}
  c=which(Fx<0.8)
  if(length(c)!=0){C=max(c)
  d=W$x[max(c)+1]-(floor(sum(W$x)*0.8)-SW[max(c)])}else{C=0
  d=W$x[1]-(floor(sum(W$x)*0.8))}
  Prop20_20=(sum(W$Group.1[-c(c,C+1)]*W$x[-c(c,C+1)])+W$Group.1[C+1]*d)/(sum(W$Group.1[a]*W$x[a])+W$Group.1[A+1]*b)

  return(Prop20_20)
}

#' @title Theil L
#'
#' @description Computes Theil_L inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of Theil_L coefficient.
#'
#' @rdname Theil_L
#'
#' @details Theil L index is defined as:
#'   \deqn{T_{L} =  T_{\alpha=0} =  \frac{1}{N} \sum_{i=1}^N ln \big(\frac{\mu }{x_{i}} \big)}
#'   where \deqn{\mu = \frac{1}{N} \sum_{i=1}^N x_{i}}
#'
#' @references Serebrenik A., van den Brand M.: Theil index for aggregation of software metrics values. 26th IEEE International Conference on Software Maintenance. IEEE Computer Society.
#' @references Conceição P., Ferreira P.: (2000) The Young Person’s Guide to the Theil Index: Suggesting Intuitive Interpretations and Exploring Analytical Applications
#' @references OECD: (2020) Regions and Cities at a Glance 2020, Chapter: Indexes and estimation techniques
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Theil_L(X)
#' Theil_L(X,W)
#'
#' data(Tourism)
#' # Theil L coefficient for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Theil_L(X,W)
#'
#'
#'
#' @export
Theil_L=function(X,W=rep(1,length(X))){1/sum(W)*sum(W*log((sum(W*X)/sum(W))/X))}



#' @title Theil T
#'
#' @description Computes `Theil_T` inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of `Theil_T` coefficient.
#'
#' @rdname Theil_T
#'
#' @details Theil T index is defined as:
#'          \deqn{T_{T} =  T_{\alpha=1} =  \frac{1}{N}  \sum_{i=1}^N  \frac{ x_{i} }{\mu} ln \big( \frac{ x_{i} }{\mu} \big)}
#'          where \deqn{\mu = \frac{1}{N} \sum_{i=1}^N x_{i}}
#'
#' @references Serebrenik A., van den Brand M.: Theil index for aggregation of software metrics values. 26th IEEE International Conference on Software Maintenance. IEEE Computer Society.
#' @references Conceição P., Ferreira P.: (2000) The Young Person’s Guide to the Theil Index: Suggesting Intuitive Interpretations and Exploring Analytical Applications
#' @references OECD: (2020) Regions and Cities at a Glance 2020, Chapter: Indexes and estimation techniques
#'
#' @examples
#' # Compare weighted and unweighted result.
#' X=1:10
#' W=1:10
#' Theil_T(X)
#' Theil_T(X,W)
#'
#' data(Tourism)
#' # Theil T coefficient for Total expenditure with sample weights:
#' X=Tourism$`Total expenditure`
#' W=Tourism$`Sample weight`
#' Theil_T(X,W)
#'
#'
#'
#' @export
Theil_T=function(X,W=rep(1,length(X))){1/sum(W)*sum(W*(X/(sum(W*X)/sum(W))*log(X/(sum(W*X)/sum(W)))))}
#' @title Sample survey on trips
#' @description Data from sample survey on trips conducted in Polish households.
#' @docType data
#' @keywords datasets
#' @name Tourism
#' @usage data(Tourism)
#' @format A data frame with 384 observations of 14 variables
#'\itemize{
#' \item Year
#' \item Group of countries
#' \item Purpose of trip
#' \item Accommodation type
#' \item Total expenditure
#' \item Expenditure for organiser
#' \item Private expenditure
#' \item Expenditure on accommodation
#' \item Expenditure on restaurants & café
#' \item Expenditure on transport
#' \item Expenditure on commodities
#' \item Number of trip's participants
#' \item Nights spent
#' \item Sample weight
#'}
#' @details Answers were modified due to disclosure control. Data presents only part of full database.
NULL
#' @title Sample survey on quality of life
#' @description Data from sample survey on quality of life conducted on Polish-Ukrainian border in 2015 and 2019.
#' @docType data
#' @keywords datasets
#' @name Well_being
#' @usage data(Well_being)
#' @format A data frame with 1197 observations of 27 variables
#'\itemize{
#' \item Area. Rural and urban
#' \item Gender. Male and female
#' \item Year. Year of survey (2015 and 2019)
#' \item V1. I have good opportunities to use my talents and skills at work
#' \item V2. I am treated with respect by others at work
#' \item V3. I have adequate opportunities for vacations or leisure activities
#' \item V4. The quality of local services where (I) live is good
#' \item V5. There is very little pollution from cars or other sources where I spend most of my time
#' \item V6. There are parks and green areas near my residence
#' \item V7. I have the freedom to plan my life the way I want to
#' \item V8. I feel safe walking around my neighborhood during the day
#' \item V9. Overall, to what extent are you currently satisfied with your life
#' \item V10. Overall, to what extent do you feel that the things you do in life are worthwhile
#' \item V11. How do you rate your health
#' \item V12. How do you rate your work
#' \item V13. How do you rate your sleep
#' \item V14. How do you rate your leisure time
#' \item V15. How do you rate your family life
#' \item V16. How do you rate your community and public affairs life
#' \item V17. How do you rate your personal plans
#' \item V18. How do you rate your housing conditions
#' \item V19. How do you rate your personal income
#' \item V20. How do you rate your personal prospects
#' \item V21. Does being part of the local community make you feel good about yourself
#' \item V22. Do you have a say in what the local community is like
#' \item V23. Is your neighborhood a good place for you to live
#' \item Weight. Sample weight for each household
#'}
#' @details Questions are on Likert scale: 1 - the worst assessment, 5 - the best assessment.
#' Only 23 question were selected out of over 100 questions.
#' Answers were modified due to disclosure control.
NULL

