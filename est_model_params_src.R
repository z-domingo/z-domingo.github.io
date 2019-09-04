install.packages("nortest")
install.packages("pracma")
install.packages("HDInterval")
library(HDInterval)
library(pracma)
library(nortest)
set.seed(2)

data = read.delim("Y6.dat", header = T, sep = "\n")

data = data[,1]
n = length(data)

summary(data)


den_dat = density(data)
x.modes = den_dat$x[findpeaks(den_dat$y)[,2]]
y.modes = findpeaks(den_dat$y)[,1] * 8000

# appears to be bi-modal suggesting 2 sub populations 
his = hist(data, freq = F, breaks = 80, xlim = c(60, 200), ylim = c(0, .03),
           xlab = "Data",
           main = "Histogram and Kernel Density of Data")

lines(den_dat)

abline(v= x.modes, col= 'red', lty = 2)

legend(180, .03, c("Kernel", "Mode"), lty = c(1, 2), col = c("black", "red"),
       box.lty = 0)

# normality tests
ad.test(data)
shapiro.test(data)


#scatterplot

plot(1:n, data, type = "p", ylim = c(68, 230), 
     xlab = "Index", ylab = "Observations",
     main = "Scatter Plot of Data")
abline(h = 123, col = "red", lty = 2)

legend(760, 230 , c("Sub-pop. Sample Mean", "Partition"), lty = c(2, 2), col = c("black", "red"),
       box.lty = 0)

# pop 1
samp.1 = data[data <= 123]
n1 = length(samp.1)
ybar1 = mean(samp.1)
s21 = var(samp.1)

prop1.est = n1 / n
abline(h = ybar1, lty = 2)

# pop 2
samp.2 = data[data > 123]
n2 = length(samp.2)
ybar2 = mean(samp.2)
s22 = var(samp.2)

prop2.est = 1 - prop1.est
abline(h = ybar2, lty = 2)

################  Gibbs sampler ################

S = 10000

### prior selections

a = 1
b = 1
p.n = 0.2
mu.0 = mean(ybar2)
t2.0 = 1500
sigma2.0 = var(samp.2)
nu.0 = 1

#initial values
theta.1 = ybar1
sigma2.1 = var(samp.1)
theta.2 = ybar2
sigma2.2 = var(samp.2)

#storage

PHI = data.frame(thetas.1 = numeric(S), thetas.2 = numeric(S), sigmas.1 = numeric(S),
                 sigmas.2 = numeric(S), ps = numeric(S), y.pred = numeric(S))

#store initial values
PHI$thetas.1[1] = theta.1
PHI$thetas.2[1] = theta.2
PHI$sigmas.1[1] = sigma2.1
PHI$sigmas.2[1] = sigma2.2
PHI$ps[1] = p.n


for (i in 2:S) {

  p.1 = p.n * dnorm(data, theta.1, sqrt(sigma2.1))
  p.2 = (1 - p.n) * dnorm(data, theta.2, sqrt(sigma2.2))
  p.bern = p.1 / (p.1 + p.2)
  x.is = rbinom(n, 1, p.bern)
  
  #partition data
  
  y1 = data[which(x.is == 1)] 
  y2 = data[which(x.is == 0)]
  
  n1 = length(y1)
  n2 = length(y2)
  
  ybar.1 = mean(y1)
  ybar.2 = mean(y2)
  sigma2.1 = var(y1)
  sigma2.2 = var(y2)
  
  p.n = rbeta(1, a + n1, b + n2)
  
  
  ####### thetas #######
  #calc params
  mun.1 = ((mu.0 / t2.0) + (n1 * ybar.1 / sigma2.1)) / (t2.0^-1 + (n1 / sigma2.1))
  t2n.1 = ((t2.0)^-1 + (n1 / sigma2.1))^-1
  
  #sample from posterior
  theta.1 = rnorm(1, mun.1, sqrt(t2n.1))
  
  #calc params
  mun.2 = ((mu.0 / t2.0) + (n2 * ybar.2 / sigma2.2)) / (t2.0^-1 + (n2 / sigma2.2))
  t2n.2 = ((t2.0)^-1 + (n2 / sigma2.2))^-1
  
  #sample from posterior
  theta.2 = rnorm(1, mun.2, sqrt(t2n.2))
  
  
  ####### sigmas #######
  nun.1 = nu.0 + n1
  sigma2n.1 = ((nu.0 * sigma2.0) + ((n1 - 1) * sigma2.1) + (n1 * (ybar.1 - theta.1)^2)) / nun.1
  
  sigma2.1 = 1 / rgamma(1, (nun.1 / 2), (sigma2n.1 * nun.1 / 2))
  
  nun.2 = nu.0 + n2
  sigma2n.2 = ((nu.0 * sigma2.0) + ((n2 - 1) * sigma2.2) + (n2 * (ybar.2 - theta.2)^2)) / nun.2
  
  sigma2.2 = 1 / rgamma(1, (nun.2 / 2), (sigma2n.2 * nun.2 / 2))
  
  
  ####### ppd #######
  
  # randomly sample 
  x.pred = runif(1)
  if (x.pred < p.n) {
    y.pred = rnorm(1, theta.1, sqrt(sigma2.1))
  } else {
    y.pred = rnorm(1, theta.2, sqrt(sigma2.2))
  }
  
  #store values
  PHI$thetas.1[i] = theta.1
  PHI$thetas.2[i] = theta.2
  PHI$sigmas.1[i] = sigma2.1
  PHI$sigmas.2[i] = sigma2.2
  PHI$ps[i] = p.n
  #will have 999 predictions
  PHI$y.pred[i-1] = y.pred
}

mean(PHI$thetas.1)
mean(PHI$thetas.2)
sqrt(mean(PHI$sigmas.1))
sqrt(mean(PHI$sigmas.2))
mean(PHI$ps)

### observations vs preditction
plot(density(data), lty = 2, ylim = c(0, .025),
     xlab = expression(y[i]),
     ylab = "Density",
     main = "Predictive Distribution vs Observations Distrubtion")
lines(density(PHI$y.pred)) 

grid()
legend(175, .025 , c("Predictions", "Observations"), lty = c(1, 2),
       box.lty = 0)


#credible intervals
hdi.thetas.1 = hdi(PHI$thetas.1, prob = c(0.025, .975))
hdi.thetas.2 = hdi(PHI$thetas.2, prob = c(0.025, .975))

hdi.sigmas.1 = hdi(sqrt(PHI$sigmas.1), prob = c(0.025, .975))
hdi.sigmas.2 = hdi(sqrt(PHI$sigmas.2), prob = c(0.025, .975))

hdi.ps = hdi(PHI$ps, prob = c(0.025, .975))

### traceplots
plot(PHI$thetas.1, type = "p", ylab = expression(theta[1]), xlab ="Iteration",
     main= "Traceplot of 10,000 Gibbs Samples")

plot(PHI$sigmas.1, type = "p", ylab = "", xlab ="Iteration",
     main= "Traceplot of 10,000 Gibbs Samples")
title(mgp = c(2.5,1,0), ylab = expression(paste(sigma[1]^2)))

#TODO: compute this and prior sensitivity 

# finding p(y_i from group 1 | y_i = 120)

p = mean(PHI$ps)
theta_1 = mean(PHI$thetas.1)
theta_2 = mean(PHI$thetas.2)
sigma2_1 = mean(PHI$sigmas.1)
sigma2_2 = mean(PHI$sigmas.2)
y.i = 120

d1 = dnorm(y.i, theta_1, sqrt(sigma2_1))
d2 = dnorm(y.i, theta_2, sqrt(sigma2_2))

p.1 = (p * d1) / ((p * d1) + ((1-p) * d2))
p.2 = ((1-p) * d2) / ((p * d1) + ((1-p) * d2))
p.n = p.1 + p.2

prob_y.i = dbinom(0, 1, (p.1 / (p.n)))
dbinom(1, 1, (p.2 / (p.n)))

dbinom(0, 1, (p.1 / (p.n)))
dbinom(0, 1, (p.2 / (p.n)))

dbinom(1, 1, (p.1 / (p.n)))
dbinom(1, 1, (p.2 / (p.n)))




### single iteration with y_i = 120
a = 1
b = 2
p.n = 0.5
mu.0 = mean(data)
t2.0 = 1500
sigma2.0 = var(data)
nu.0 = 1

#initial values
theta.1 = ybar1
sigma2.1 = var(samp.1)
theta.2 = ybar2
sigma2.2 = var(samp.2)

#storage

PHI = data.frame(thetas.1 = numeric(S), thetas.2 = numeric(S), sigmas.1 = numeric(S),
                 sigmas.2 = numeric(S), ps = numeric(S), y.pred = numeric(S))

#store initial values
PHI$thetas.1[1] = theta.1
PHI$thetas.2[1] = theta.2
PHI$sigmas.1[1] = sigma2.1
PHI$sigmas.2[1] = sigma2.2
PHI$ps[1] = p.n


  
  p.1 = p.n * dnorm(y.i, theta.1, sqrt(sigma2.1))
  p.2 = (1 - p.n) * dnorm(y.i, theta.2, sqrt(sigma2.2))
  p.bern = p.1 / (p.1 + p.2)
  x.is = rbinom(1, 1, p.bern)
  
  #partition data
  
  y1 = data[which(x.is == 1)] 
  y2 = data[which(x.is == 0)]
  
  n1 = length(y1)
  n2 = length(y2)
  
  ybar.1 = mean(y1)
  ybar.2 = mean(y2)
  sigma2.1 = var(y1)
  sigma2.2 = var(y2)
  
  p.n = rbeta(1, a + n1, b + n2)
  
  
  ####### thetas #######
  #calc params
  mun.1 = ((mu.0 / t2.0) + (n1 * ybar.1 / sigma2.1)) / (t2.0^-1 + (n1 / sigma2.1))
  t2n.1 = ((t2.0)^-1 + (n1 / sigma2.1))^-1
  
  #sample from posterior
  theta.1 = rnorm(1, mun.1, sqrt(t2n.1))
  
  #calc params
  mun.2 = ((mu.0 / t2.0) + (n2 * ybar.2 / sigma2.2)) / (t2.0^-1 + (n2 / sigma2.2))
  t2n.2 = ((t2.0)^-1 + (n2 / sigma2.2))^-1
  
  #sample from posterior
  theta.2 = rnorm(1, mun.2, sqrt(t2n.2))
  
  
  ####### sigmas #######
  nun.1 = nu.0 + n1
  sigma2n.1 = ((nu.0 * sigma2.0) + ((n1 - 1) * sigma2.1) + (n1 * (ybar.1 - theta.1)^2)) / nun.1
  
  sigma2.1 = 1 / rgamma(1, (nun.1 / 2), (sigma2n.1 * nun.1 / 2))
  
  nun.2 = nu.0 + n2
  sigma2n.2 = ((nu.0 * sigma2.0) + ((n2 - 1) * sigma2.2) + (n2 * (ybar.2 - theta.2)^2)) / nun.2
  
  sigma2.2 = 1 / rgamma(1, (nun.2 / 2), (sigma2n.2 * nun.2 / 2))
  
  
  ####### ppd #######
  
  # randomly sample 
  x.pred = runif(1)
  if (x.pred < p.n) {
    y.pred = rnorm(1, theta.1, sqrt(sigma2.1))
  } else {
    y.pred = rnorm(1, theta.2, sqrt(sigma2.2))
  }
  
  #store values
  PHI$thetas.1[i] = theta.1
  PHI$thetas.2[i] = theta.2
  PHI$sigmas.1[i] = sigma2.1
  PHI$sigmas.2[i] = sigma2.2
  PHI$ps[i] = p.n
  #will have 999 predictions
  PHI$y.pred[i-1] = y.pred


mean(PHI$thetas.1)
mean(PHI$thetas.2)

