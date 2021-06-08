#'  T-test for comparing the element in two log ratio list
#'
#'
#' @param i ith species/genes which are used for comparison
#' @param z the first samples
#' @param w the second samples
#'
#'
#' @return Vector of t-test p-value
#'
#' @examples
#' ecc.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==0),],version=1,epsilon=1)})
#' free.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==1),],version=1,epsilon=1)})
#' ttest(1,ecc.ratio,free.ratio)
#'
#'
#'
#'
#' @export
#'
ttest <- function(i,z,w)
{
  x = z[[i]]
  y = w[[i]]
  x1 = x[x!= -Inf & x!= Inf]
  y1 = y[y!= -Inf & y!= Inf]
  if(length(x1) <= 1 | length(y1) <= 1)
  {
    res <- NA
  }
  else if(all(x1==x1[1]))
  {
    res <- t.test(y1-x1[1])$p.value
  }
  else if(all(y1==y1[1]))
  {
    res <- t.test(x1-y1[1])$p.value
  }
  else
  {
    res <- t.test(x1,y1)$p.value
  }
  return(res)
}


#' Pre-test for DNA and RNA
#'
#' @param i ith genera or species specific taxonomic units which are used for testing
#' @param otu Operational taxonomic unit table
#' @param method the test we want to used for pre-analysis, \code{method = "LB"} represents logistics Beta test, while \code{method = "LN"} represents log Normal test
#' @param val DNA & RNA that needs to be tested
#' @param xm the main predictor we use (consider age and batch as covariates)
#' @param ranef whether we should put batch as a random effect
#'
#'
#' @return Coefficient estimate and p-value of the main predictor
#'
#' @examples
#'
#' pure.pretest(1,otu)
#'
#'
#'
#' @importFrom lmerTest lmer
#' @importFrom gamlss gamlss
#' @export
#'
pure.pretest <- function(i,otu,method = "LN",val = "DNA",xm="caries_free",ranef = FALSE)
{
  df1 <- data.frame(df,DNA = otu[i,,1],RNA = otu[i,,2])
  if(method == "LN")
  {
    if(ranef)
    {
      modf <- lmerTest::lmer(log(1+get(val))~get(xm)+age+(1|batch),data=df1)
      x= summary(modf)$coefficients[2,c(1,5)]
    }
    else{
      modf <- lm(log(1+get(val))~get(xm)+age+batch,data=df1)
      x = summary(modf)$coefficients[2,c(1,4)]
    }
    names(x) = c("coeffcient of main predictor","p-value of LN test")
  }
  else if (method=="LB")
  {
    scals <- function(x){x/sum(x)}
    if(sum(df1[,paste0(val)]!=0) <= 1){
      x = rep(NA,4);
    }
    else
    {
      if(ranef)
      {
        modlb <- gamlss::gamlss(scals(get(val))~get(xm)+age+random(batch),data = df1,
                        nu.formula = ~ get(xm)+age+random(batch), family = BEZI(sigma.link = "log"), control = gamlss.control(n.cyc = 100, trace = FALSE))

      }
      else{
        modlb <- gamlss(scals(get(val))~get(xm)+age+batch,data = df1,
                        nu.formula = ~ get(xm)+age+batch, family = BEZI(sigma.link = "log"), control = gamlss.control(n.cyc = 100, trace = FALSE))

      }
        lb1<-summary(modlb, save=TRUE);
        x = as.vector(lb1$coef.table[which(grepl("get",rownames(lb1$coef.table))==TRUE),c(1,4)]);
        names(x) = c("coeffcient of main predictor on non-zero part","coeffcient of main predictor on zero part","p-value of LB test on nonzero part","p-value of LB test on zero part")
    }
  }

  return(x)

}


#' Correlation-test between log ratio and predictor
#'
#' @param i ith genera or species specific taxonomic units which are used for testing
#' @param simu Operational taxonomic unit table
#' @param all Log ratio list containing every genera or species
#' @param xm the main predictor we use
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#'
#' @return Coefficient estimate and p-value of the main predictor
#'
#' @examples
#' all.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,,],version=1,epsilon=1)})
#' cor.otu(1,otu,all.ratio)
#'
#'
#'
#' @export
#'
cor.otu <- function(i,simu,all,xm = "t3c", method = "pearson")
{
  df1 <-  newdf(i,simu,all)
  if(nrow(df1)==0)
  {
    x = rep(NA,2)
  }
  else
  {
    p <- with(df1,cor.test(get(xm),ratio, method = method))
    x = c(p$estimate,p$p.value)
  }
  names(x) = c("estimated corr","p-value")
  return(x)
}

#' Ratio level test
#'
#' @param i ith genera or species specific taxonomic units which are used for testing
#' @param simu Operational taxonomic unit table
#' @param all Log ratio list containing every genera or species
#' @param xm the main predictor we use
#' @param reg kind of regression that we want to see the p-value of it, it has \code{linear} and \code{logistics} two choices. For \code{logistics} the outcome would become whether RNA equals to zero.
#' @param ranef whether we should regard batch as random effect
#' @param df meta data we use
#' @param firth_correction whether we need Firth Logistics Regression for \code{reg="logistics"}
#'
#' @return Coefficient estimate and p-value of the main predictor
#'
#' @examples
#' all.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,,],version=1,epsilon=1)})
#' ratio.test(1,otu,all.ratio,reg="linear",df=df)
#'
#' @importFrom lmerTest lmer
#' @importFrom lme4 glmer
#' @importFrom logistf logistf
#'
#' @export
#'
ratio.test <- function(i,simu,all,xm="t3c",reg = "linear",ranef = TRUE,df = df, firth_correction = FALSE)
{
  df1 <- newdf(i,simu,all,df = df)
  df1$ratioindex <- ifelse(simu[i,which(simu[i,,1]>0),2] == 0, 0, 1)
  if(nrow(df1)<=4 |(sum(table(df1$caries_free)==0) > 0 & xm == "caries_free"))
  {
    x= rep(NA,2)

  }else {
    if(sum(table(df1$batch)==0) > 0)
    {
      if(reg == "linear")
      {
        if(ranef == TRUE)
        {
          stop("No random batch effect")
        }
        else
        {
          mod <- lm(ratio ~get(xm)+age, data = df1)
        }
        x = summary(mod)$coefficients[2,c("Estimate","Pr(>|t|)")]

      }
      else if(reg == "logistics")
      {
        if(ranef == TRUE)
        {
          stop("No random batch effect")
        }
        else
        {
          if(xm == "caries_free" & firth_correction ==TRUE)
          {
            mod <-  logistf::logistf(ratioindex ~ get(xm)+batch+age,data = df1)
            x = c(mod$coefficients[2],mod$prob[2])
          }
          else{
            mod <- glm(ratioindex ~get(xm)+age,family="binomial", data = df1)
            x = summary(mod)$coefficients[2,c("Estimate","Pr(>|z|)")]

          }

        }
      }
    }
    else{
      if(reg == "linear")
      {
        if(ranef == TRUE)
        {
          mod <- lmerTest::lmer(ratio ~get(xm)+age+(1|batch), data = df1)
        }
        else
        {
          mod <- lm(ratio ~get(xm)+age+batch, data = df1)
        }
        x = summary(mod)$coefficients[2,c("Estimate","Pr(>|t|)")]
      }
      else if(reg == "logistics")
      {
        if(ranef == TRUE)
        {
          mod <- lme4::glmer(ratioindex ~get(xm)+age+(1|batch),family=binomial(), data = df1)
          x = summary(mod)$coefficients[2,c("Estimate","Pr(>|z|)")]
        }
        else
        {
          if(xm == "caries_free" & firth_correction ==TRUE)
          {
            mod <-  logistf::logistf(ratioindex ~ get(xm)+batch+age,data = df1)
            x = c(mod$coefficients[2],mod$prob[2])
          }
          else{
            mod <- glm(ratioindex ~get(xm)+age+batch,family="binomial", data = df1)
            x = summary(mod)$coefficients[2,c("Estimate","Pr(>|z|)")]
          }
        }


      }
    }
  }
  names(x) = c("coefficient estimates","p-value")
  return(x)
}


#' Gaussian Mixture Regression for otu data
#'
#' @param i ith genera or species specific taxonomic units which are used for testing
#' @param simu Operational taxonomic unit table
#' @param all Log ratio list containing every genera or species
#' @param formula0 a two-sided linear formula object describing all parts of the model. The main predictor should be written in the first place after tilde.
#' @param bindex the peak status of ith genera or species of otu, it has three values: \code{left},\code{right} and \code{both}, which refers to left-only peak, right-only peak and two peaks.
#'
#' @return Coefficient estimate and p-value of each peak
#'
#' @examples
#' all.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,,],version=1,epsilon=1)})
#' bindex = sapply(1:nrow(otu),function(i){ifelse(mean(otu[i,otu[i,,1]!=0,2] > 0) < 0.15,"left",ifelse(mean(otu[i,otu[i,,1]!=0,2] == 0)<0.15,"right","both"))})
#' GMR.peaktest(3,otu,all.ratio,formula0 = "ratio~t3c+age",bindex,df=df)
#'
#' @importFrom flexmix stepFlexmix FLXMRglm
#'
#' @export
#'

GMR.peaktest <- function(i,simu,all,formula0 = "ratio ~ t3c + age",bindex,df = df)
{
  df1 <- newdf(i,simu,all,df = df)
  fitted = flexmix::stepFlexmix(formula(formula0), model = flexmix::FLXMRglm(family = "gaussian"), nrep = 3, k = 1:2,
                         data = df1)
  fitted1= getModel(fitted,which="BIC")
  fitted.lab <- relabel(fitted1,"model",which = "(Intercept)")
  up <- refit(fitted.lab)
  coef.matrix <- up@components[[1]]
  if(is.null(coef.matrix$Comp.2))
  {
    if(bindex[i] == "left")
    {
      return(c(coef.matrix$Comp.1[2,c(1,4)],rep(NA,2)))
    }
    else
    {
      return(c(rep(NA,2),coef.matrix$Comp.1[2,c(1,4)]))
    }
  }
  else
  {
    return(c(coef.matrix$Comp.1[2,c(1,4)],coef.matrix$Comp.2[2,c(1,4)]))
  }
}

#' AUC probability and empirical p-value calculation
#'
#' @param x numeric vector represents 1st sample
#' @param y numeric vector represents 2nd sample
#' @param N sampling times with replacement
#' @param rdm random seed
#' @param p.output indicator whether to return p-value
#'
#' @return AUC probability estimates (or with p-value)
#'
#' @examples
#' ecc.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==0),],version=1,epsilon=1)})
#' free.ratio =  sapply(1:nrow(otu),function(x){log_ratio(otu[x,which(meta$cariesfree==1),],version=1,epsilon=1)})
#' auc_prob_pvalue(ecc.ratio[[3]],free.ratio[[3]],rdm = i)
#'
#' @export
#'

auc_calculate <- function(x,y,N=5000,rdm,p.output = FALSE)
{
  set.seed(rdm)
  pos <- sample(x, N, replace=TRUE)
  neg <- sample(y, N, replace=TRUE)
  xhat = (sum(pos > neg) + sum(pos == neg)/2) / N
  x = xhat
  if(p.output == TRUE)
  {
    g1hat <- sapply(1:N,function(i){(sum(pos > neg[i]) + sum(pos == neg[i])/2) /N - xhat})
    g2hat <- sapply(1:N,function(i){(sum(pos[i] > neg) + sum(pos[i] == neg)/2) /N - xhat})
    sigm <- mean(g1hat^2+g2hat^2)/(N-1)
    pval <- 2-2*pnorm(abs(xhat-0.5)/sqrt(sigm))
    x = c(xhat, pval)
  }
   return(x)

}
