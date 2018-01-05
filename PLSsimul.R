rm(list = ls())
if (!require(clusterGeneration)) {install.packages('clusterGeneration'); library(clusterGeneration)} 
if (!require(CCA)) {install.packages('CCA'); library(CCA)}
if (!require(MASS)) {install.packages('MASS'); library(MASS)}
if (!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if (!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
if (!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if (!require(pls)) {install.packages('pls'); library(pls)}
if (!require(spls)) {install.packages('spls'); library(spls)}
if (!require(spls)) {install.packages('spls'); library(spls)}
if (!require(spls)) {install.packages('spls'); library(spls)}

# Monte Carlo simulation PLS
# Generate random covariance matrices with eig
# Sigma =
# rotate Y?

r.XY <- 0.2 # correlation between X and Y

n.X <- 5 # dimension of X
n.Y <- 5  # dimension of Y
n = n.X + n.Y

nfactors.X <- 3 # number of factors for X (nfactors.X < n.X)
nfactors.Y <- 3 # number of factors for Y (nfactors.Y < n.Y)

generate_sigma = function(n.X, n.Y, r.XY, wn){ 
eigenvalues.X <- c(1, runif(n.X - 1, 0, 0.4))
eigenvalues.Y <- c(1, runif(n.Y - 1, 0, 0.4))
n = n.X + n.Y

# eigenvalue matrix, ensure a correlation between factors
A.X = diag(eigenvalues.X)
A.Y = diag(eigenvalues.Y)
A = bdiag(diag(c(eigenvalues.X, eigenvalues.Y)))
A[1, n.X + 1] = A[n.X + 1, 1] = r.XY

# generate random orthogonal matrix and take nfactors columns
U.X <- qr.Q(qr(matrix(rnorm(n.X^2), n.X)))
U.Y <- diag(n.Y)#qr.Q(qr(matrix(rnorm(n.Y^2), n.Y)))
U <- bdiag(U.X, U.Y)

#U <- qr.Q(qr(matrix(rnorm(n^2), n)))

# rotate eigenvalue matrix
Sigma <- t(U) %*% A %*% U + wn*diag(n)
#Sigma <- A

}

ITERS = 5
STEPS = 5
STEP = 100
etas = seq(0.1,0.9,0.4)

wns = 0.1*seq(5) # white noise variance
mycor.test = array(0, dim = c(STEPS, ITERS, length(wns)))
mycor.train = array(0, dim = c(STEPS, ITERS, length(wns)))
PROPORTION = 0.8 # proportion of train data
  
for(j in seq(length(wns))){
for (iter in seq(ITERS)){
  wn = wns[j]
  Sigma = generate_sigma(n.X, n.Y, r.XY, wn)
  print(iter)
  for (i in seq(STEPS)){
  NSAMP = i*STEP
  S = mvrnorm(n = NSAMP, rep(0, n.X + n.Y), Sigma)
  X = S[, 1:n.X] 
  Y = S[, (n.X + 1):n] 

  MYSAMPLE = seq(round(NSAMP*PROPORTION))
  X.train = X[MYSAMPLE, ]
  X.test = X[-MYSAMPLE, ]
  Y.train = Y[MYSAMPLE, ]
  Y.test = Y[-MYSAMPLE, ]
  
  cv <- cv.spls( X.train, Y.train, eta = etas, K = c(1:n.X) )
  mypls <- spls( X.train, Y.train, eta = cv$eta.opt, K = cv$K.opt )
  

  y.pred.train = predict(mypls, X.train)
  y.pred.test = predict(mypls, X.test)
  mycor.train[i, iter, j] = cor.test(Y.train[, 1], y.pred.train[, 1])$estimate
  mycor.test[i, iter, j] = cor.test(Y.test[, 1], y.pred.test[, 1])$estimate
}
}
}
mycor.melt.test = melt(mycor.test) 
mycor.melt.train = melt(mycor.train) 
mycor.melt.test$Type = 'test'
mycor.melt.train$Type = 'training'

mycor.melt = rbind(mycor.melt.test, mycor.melt.train)

colnames(mycor.melt) = c('Step', 'Iteration', 'Noise', 'r', 'Type')
mycor.melt = mycor.melt %>% mutate( Samples = Step*STEP)


mycor.sum = mycor.melt %>%  
  group_by(Samples, Noise, Type) %>% 
  summarise(lwr = quantile(r, .05), 
            upr = quantile(r, .95), 
            r = mean(r))

print(
ggplot(mycor.sum, aes(x = Samples, y = r)) + geom_line() + geom_point() + 
   geom_ribbon(data = mycor.sum, aes(ymin = lwr, ymax = upr), alpha = 0.6) +
   theme_classic() + geom_hline(yintercept = r.XY, col = 'red') + 
   geom_hline(yintercept = 0) + ylim(-1, 1) +
   facet_grid( Noise ~ Type) 
)
