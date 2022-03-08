rm(list=ls())

library(sgPLS)
library(abind)

# Loading the data
data(linnerud)
dataX <- linnerud$exercise
dataY <- linnerud$physiological

# Creating traning/test sets
set.seed(1)
ids_train=sample(1:nrow(dataX),size=0.8*nrow(dataX))
ids_test=c(1:nrow(dataX))[!c(1:nrow(dataX))%in%ids_train]
dataX_train=dataX[ids_train,]
dataY_train=dataY[ids_train,]
dataX_test=dataX[ids_test,]

# Fitting pls model on training set
mypls <- pls(dataX_train, dataY_train)


### Prediction

# Arguments
newdata=dataX_test # dataset to use for predictions
ncomp=2 # number of components to use for predictions

# Extract scaled X and Y datasets used for model fitting
X=mypls$X
Y=mypls$Y

# Average and standard deviation of un-scaled X dataset used for model fitting
means.X=attr(mypls$X,"scaled:center")
sigma.X=attr(mypls$X,"scaled:scale")

# Scaling new X dataset with average and standard deviation of un-scaled X dataset used for model fitting
newdata=t(apply(newdata,1,FUN=function(x){(x-means.X)/sigma.X}))

# Compute P, C and W matrices that will be used to get hat(Y)=XW(P'W)^-1C'
Pmat = crossprod(as.matrix(X), mypls$variates$X)
Cmat = crossprod(as.matrix(Y), mypls$variates$X)
Wmat = mypls$loadings$X

# Average and standard deviation of un-scaled Y dataset used for model fitting
means.Y=attr(mypls$Y,"scaled:center")
sigma.Y=attr(mypls$Y,"scaled:scale")

# Compute Y prediction
Yhat=NULL
for (i in 1:ncomp){
  Ypred=newdata%*%Wmat[,1:i]%*%solve(t(Pmat[,1:i])%*%Wmat[,1:i])%*%t(Cmat)[1:i,]
  Ypred=apply(Ypred, 1, function(x){x*sigma.Y+means.Y})
  Ypred=t(Ypred)
  Yhat=abind(Yhat,Ypred,along=3)
}
Yhat

# Compare to predictions from the predict() function
mypredictions=predict(mypls, newdata=dataX_test)
mypredictions$predict[,,1]

par(mar=c(5,5,1,1))
plot(Yhat[,,1], mypredictions$predict[,,1], pch=19, col="navy", las=1, cex.lab=1.5,
     xlab=expression(hat(Y)[1]~"(manual)"), ylab=expression(hat(Y)[1]~"(implemented function)"))
abline(0,1,lty=3)



