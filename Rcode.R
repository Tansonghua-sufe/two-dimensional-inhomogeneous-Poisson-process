pot.external <- function(data, threshold,external.regressors){
  # fit a two-dimensional inhomogeneous Poisson process model with Explanatory Variable Model 
  # please refer to Chapter 7.7.6 in Tsay, R. S. (2005). Analysis of financial time series. John wiley & sons.
  # This function is modified based on evir::pot()
  # ===========================================
  # data(Vector)==>numeric vector of data
  # threshold(double)==>a threshold value
  # external.regressors(matrix)==>Explanatory Variable
  # ===========================================
  p <- dim(external.regressors)[2] # dimension of Explanatory Variable
  if(length(data)!=dim(external.regressors)[1])
    stop("The number of rows in external regressors is not equal to the length of data") # 检查数据长度是否相同
  n <- length(as.numeric(data))
  times <- 1:n
  attributes(data)$times <- times
  start <- 1
  end <- n
  span <- end - start # D
  exceedances.its <- structure(data[data > threshold], times =times[data > threshold])
  exceedances.external <- as.matrix(external.regressors[data > threshold,]) 
  n.exceed <- length(as.numeric(exceedances.its)) 
  time.exceed <- times[data > threshold] 
  p.less.thresh <- 1 - n.exceed/n 
  exceedances <- as.numeric(exceedances.its) 
  
  # ==================
  # init
  xbar <- mean(exceedances) - threshold
  s2 <- var(exceedances)
  shape0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  extra <- ((length(exceedances)/span)^( - shape0) - 1)/shape0
  betahat <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  scale0 <- log(betahat/(1 + shape0 * extra))
  loc0 <- 0
  theta <- c(shape0,rep(0,p), scale0,rep(0,p), loc0,rep(0,p)) #初值
  
  # ==================
  # Negative log-likelihood function
  negloglik <- function(theta, exceedances,exceedances.external,external.regressors, threshold, span)
  {
    p <- (length(theta)/3-1)
    n.exceed <- length(exceedances)
    
    loc_t <- theta[2*(p+1)+1]+external.regressors%*%theta[(2*(p+1)+2):(3*(p+1))]# 位置参数,教材中的beta_t
    shape_t <- theta[1]+external.regressors%*%theta[2:(p+1)]# 形状参数,教材中的xi_t
    scale_t <- exp(theta[(p+1)+1]+external.regressors%*%theta[((p+1)+2):(2*(p+1))])# 尺度参数,教材中的alpha_t
    
    loc_exceed_t <- theta[2*(p+1)+1]+exceedances.external%*%theta[(2*(p+1)+2):(3*(p+1))] 
    shape_exceed_t <- theta[1]+exceedances.external%*%theta[2:(p+1)]
    scale_exceed_t <- exp(theta[(p+1)+1]+exceedances.external%*%theta[((p+1)+2):(2*(p+1))])
    
    if(scale_exceed_t<0||min(1 + (shape_exceed_t * (exceedances -loc_exceed_t)) / scale_exceed_t) <= 0){
      f <- 1e+6
    }else{
      # eq(7.41) in Tsay(2005)
      y <- logb(1 + (shape_exceed_t * (exceedances - loc_exceed_t)) / scale_exceed_t)
      term3 <- (1/shape_exceed_t + 1)*(y)
      term1 <- (1 + (shape_t * (threshold - loc_t)) /
                  scale_t)^(-1/shape_t)/span
      term2 <- logb(scale_exceed_t)
      term4 <- logb(span)*n.exceed
      f <- sum(term2 + term3) + sum(term1)+term4
    }
    f
  }
  fit <- optim(theta, negloglik, hessian = TRUE, exceedances =
                 exceedances, threshold = threshold,exceedances.external=exceedances.external,external.regressors=external.regressors, span = span,control = list(maxit=10000)) # 使用optim()最小化负对数似然函数
  if(fit$convergence)
    warning("optimization may not have succeeded")
  par.ests <- fit$par
  # varcov <- solve(fit$hessian)
  # par.ses <- sqrt(diag(varcov))
  
  loc.est = par.ests[(2*(p+1)+1):(3*(p+1))]
  scale.est = par.ests[(1*(p+1)+1):(2*(p+1))]
  shape.est = par.ests[(0*(p+1)+1):(1*(p+1))]
  
  # ==============================
  # model check
  
  # Exceedance Rate
  loc_t <- loc.est[1]+external.regressors%*%loc.est[2:(p+1)]
  shape_t <- shape.est[1]+external.regressors%*%shape.est[2:(p+1)]
  scale_t <- exp(scale.est[1]+external.regressors%*%scale.est[2:(p+1)])
  temp1 <- (sapply(1+shape_t*(threshold-loc_t)/scale_t,FUN=function(x) max(x,0)))^(-1/shape_t)
  z <- NULL
  z[1] <- sum(temp1[1:time.exceed[1]])/span
  for(i in 2:n.exceed){
    z[i] <- sum(temp1[(time.exceed[i-1]+1):time.exceed[i]])/span #教材中(7.42)
  }
  par(mfrow=c(2,2))
  qplot(z)
  
  ## Distribution of Excesses
  loc_t <- loc.est[1]+exceedances.external%*%loc.est[2:(p+1)]
  shape_t <- shape.est[1]+exceedances.external%*%shape.est[2:(p+1)]
  scale_t <- exp(scale.est[1]+exceedances.external%*%scale.est[2:(p+1)])
  beta_t <- scale_t + shape_t * (threshold - loc_t)
  
  w <- 1/shape_t*log(sapply((1+shape_t*(exceedances-threshold)/beta_t),FUN = function(x) max(x,0)))
  qplot(w)
  
  
  # Independence
  acf(z) 
  acf(w)
  par(mfrow=c(1,1))
  
  # ==============================
  # output result
  out <- list(n = length(data), period = c(start, end), data = exceedances.its, 
              span = span, threshold = threshold, p.less.thresh = p.less.thresh, 
              n.exceed = n.exceed,
              loc.est=loc.est,
              shape.est=shape.est,
              scale.est=scale.est,
              converged = fit$convergence
  )
  class(out) <- "potd"
  return(out)
}