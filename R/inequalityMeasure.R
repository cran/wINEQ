#' @title Allison and Foster
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' AF(X,W)
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




#' @title Atkinson
#'
#' @description Computes Atkinson inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param e is a parameter for calculating the value of the Atkinson coefficient
#'
#' @return The value of Atkinson coefficient.
#'
#' @rdname Atkinson
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Atkinson(X,W)
#'
#' @export
Atkinson=function(X,W=rep(1,length(X)),e=1)
{
  if(e==1){A=1-prod(X^(W/sum(W)))/ (sum(W*X)/sum(W))}else{
    A=1-((1/sum(W)*sum(W*(X^(1-e))))^(1/(1-e)))/(sum(W*X)/sum(W))}
  return(A)
}


#' @title Entropy
#'
#' @description Computes Entropy inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param parameter is a entropy parameter
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of Entropy coefficient.
#'
#' @rdname Entropy
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Entropy(X,W)
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





#' @title Kolm
#'
#' @description Computes Kolm inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param parameter is a Kolm parameter
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of Kolm coefficient.
#'
#' @rdname Kolm
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Kolm(X,W)
#'
#' @export
Kolm=function (X,W=rep(1,length(X)), parameter = 1, na.rm = TRUE)
{
  if (!na.rm && any(is.na(X)))
    return(NA_real_)
  X <- as.numeric(na.omit(X))
  X <- as.numeric(X)
  if (is.null(parameter))
    parameter <- 1
  KM <- parameter * ((sum(W*X)/sum(W)) - X)
  KM <- sum(W*(exp(KM)))/sum(W)
  KM <- (1/parameter) * log(KM)
  KM
}




#' @title RicciSchutz
#'
#' @description Computes RicciSchutz inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param na.rm logical, should missing values (NAs) be removed prior to computations? If set to FALSE the computations yield NA
#'
#' @importFrom stats na.omit
#'
#' @return The value of RicciSchutz coefficient.
#'
#' @rdname RicciSchutz
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' RicciSchutz(X,W)
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




#' @title CoefVar
#'
#' @description Computes CoefVar inequality measure of a given variable taking into account weights.
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' CoefVar(X,W)
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


#' @title Gini
#'
#' @description Computes Gini inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of Gini coefficient.
#'
#' @rdname Gini
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Gini(X,W)
#'
#' @export
Gini=function(X,W=rep(1,length(X))){sum(abs(matrix(X,length(X),length(X),TRUE)-matrix(X,length(X),length(X),FALSE))*matrix(W,length(W),length(W),TRUE)*matrix(W,length(W),length(W),FALSE)/(2*sum(W)^2*(sum(W*X)/sum(W))))}


#' @title Hoover
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Hoover(X,W)
#'
#' @export
Hoover=function(X,W=rep(1,length(X))){(1/2)*(sum(W*abs(X - (sum(W*X)/sum(W))))/sum(X*W))}


#' @title Leti
#'
#' @description Computes Leti inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @importFrom stats aggregate
#'
#' @return The value of Leti coefficient.
#'
#' @rdname Leti
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Leti(X,W)
#'
#' @export
Leti=function(X,W=rep(1,length(X)))
{
  W=W[is.na(W)==FALSE]
  W=aggregate(W,by=list(X),FUN=sum)
  W=W$x
  Fx=(cumsum(W)/sum(W))
  Leti=2*sum(Fx*(1-Fx))
  LetiNorm=Leti*2/(length(Fx)-1)
  return(LetiNorm)
}


#' @title Jenkins and Cowell_and_Flachaire
#'
#' @description Computes Jenkins and Cowell_and_Flachaire inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#' @param alfa is the Jenkins coefficient parameter
#'
#' @importFrom stats aggregate
#'
#' @return The value of Jenkins and Cowell_and_Flachaire coefficient.
#'
#' @rdname Jenkins
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'             F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Jenkins(X,W)
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




#' @title Palma
#'
#' @description Palma proportion - the ratio of the total income of the 10% richest people to the 40% poorest people.
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'             Putting the Gini Back in the Bottle? 'The Palma' as a Policy-Relevant Measure of Inequality
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Palma(X,W)
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

#' @title Prop20:20
#'
#' @description 20:20 ratio - the ratio of thr total income of the 20% richest people to the 20% poorest people.
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'             Panel Data Econometrics: Theoretical Contributions And Empirical Applications edited by Badi Hani Baltag
#'             Notes on Statistical Sources and Methods - The Equality Trust.
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Prop20_20(X,W)
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
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'             A. Serebrenik, M. van den Brand. Theil index for aggregation of software metrics values. 26th IEEE International Conference on Software Maintenance. IEEE Computer Society.
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Theil_L(X,W)
#'
#' @export
Theil_L=function(X,W=rep(1,length(X))){1/sum(W)*sum(W*log((sum(W*X)/sum(W))/X))}



#' @title Theil T
#'
#' @description Computes Theil_T inequality measure of a given variable taking into account weights.
#'
#' @param X is a data vector
#' @param W is a vector of weights
#'
#' @return The value of Theil_T coefficient.
#'
#' @rdname Theil_T
#'
#' @references Philip B. Coulter: (1989) Measuring Inequality ISBN 0-8133-7726-9
#'             A. Serebrenik, M. van den Brand. Theil index for aggregation of software metrics values. 26th IEEE International Conference on Software Maintenance. IEEE Computer Society.
#'
#' @examples
#' X=c(1,2,3,4,5,6,7,8,9)
#' W=c(2,5,6,7,3,4,5,2,5)
#' Theil_T(X,W)
#'
#' @export
Theil_T=function(X,W=rep(1,length(X))){1/sum(W)*sum(W*(X/(sum(W*X)/sum(W))*log(X/(sum(W*X)/sum(W)))))}
