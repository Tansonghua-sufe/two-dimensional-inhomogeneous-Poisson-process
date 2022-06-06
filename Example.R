# This is an illustration, i.e. the exercises 7.2(f) in Tsay(2005)
Cisco <- read.table("d-csco9808.txt",header = T)
Cisco$date <- as.Date(as.character(Cisco$date),format = "%Y%m%d")

span <- 1:length(Cisco$rtn)
index_year <- (year(Cisco$date)-min(year(Cisco$date)))/(max(year(Cisco$date))-min(year(Cisco$date))) # an annual time trend
index_month <- as.numeric(month(Cisco$date)%in%c(10,11,12)) # a dummy variable for October, November, and December
# a fitted volatility based on a Gaussian GARCH(1,1) model.
spec.GARCH.norm=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                           mean.model = list(armaOrder = c(0, 0),include.mean = FALSE),
                           distribution.model = "norm")
Cisco.GARCH.norm <- ugarchfit(spec = spec.GARCH.norm, data = -Cisco$rtn, solver = 'hybrid', fit.control = list(stationarity = 1))

# fit a two-dimensional inhomogeneous Poisson process model
Cisco.pot.external <- pot.external(data = -Cisco$rtn,threshold = 0.02,external.regressors = cbind(index_year,index_month,sigma(Cisco.GARCH.norm)))

# forecast based on equation (7.36) in Tsay(2005)
n <- length(index_year)
beta <- c(1,index_year[n],index_month[n],sigma(Cisco.GARCH.norm)[n])%*%Cisco.pot.external$loc.est
xi <- c(1,index_year[n],index_month[n],sigma(Cisco.GARCH.norm)[n])%*%Cisco.pot.external$shape.est
alpha <- exp(c(1,index_year[n],index_month[n],sigma(Cisco.GARCH.norm)[n])%*%Cisco.pot.external$scale.est)
beta-alpha/xi*(1-(-span*log(1-tau))^(-xi))