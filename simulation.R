
source("functions.R")

library(MASS)
library(mvtnorm)
library(OpenMx)
library(hydroGOF)
library(statmod)
library(lme4)


N <- 300

Rho <-  matrix( 
  c( 1, .3, .4 , .1,
     .3, 1, .3, .3,
     .4, .3, 1, .3,
     .1,  .3,  .3, 1
  ),  
  nrow=4,              
  ncol=4,              
  byrow = TRUE)        

sigma <- c(1,1,4,5)
#Sigma <- sigma *  Rho 
Sigma <- Rho * sigma * rep(sigma, each = nrow(Rho))

mu = c(0,0,0,0)
e <- rmvnorm(N,mu,Sigma)

x1 <- rnorm(N,0,1)
x2 <- rnorm(N,0,1)
x3 <- rnorm(N,0,1)
x4 <- rnorm(N,0,1)

# add the new columns  #

x1.1 <- rnorm(N,0,1)
x2.1 <- rnorm(N,0,1)
x3.1 <- rnorm(N,0,1)
x4.1 <- rnorm(N,0,1)

 
X1<-cbind(1,x1, x1.1)
X2<-cbind(1,x2, x2.1)
X3<-cbind(1,x3, x3.1)
X4<-cbind(1,x4, x4.1)



Beta = matrix( 
  c( .2, .1, .2, .2,  0.3 , .1,
     .5, 8, 10, .4,  5, 8  ),  
  nrow=4,              
  ncol=3,              
  byrow = TRUE)        

para_true <-c( t (Beta[1,]) , t(Beta[2,]) , t(Beta[3,]), t(Beta[4,]) , sigma[3] ,sigma[4] , Rho[1,2], Rho[1,3], Rho[1,4], Rho[2,3], Rho[2,4], Rho[3,4]   )


y1 = X1 %*% Beta[1,] 
y2 = X2 %*% Beta[2,] 
y3 = X3 %*% Beta[3,] 
y4 = X4 %*% Beta[4,] 



Y<-cbind (y1,y2,y3,y4) + e


nrows <- c(1:N)

for (val in nrows) {
  Y[val,1] <- if( Y[val,1] > 0) 1 else 0
  Y[val,2] <- if( Y[val,2] > 0) 1 else 0
}


colnames(Y) <- c("Y1","Y2","Y3","Y4")
 
 data_1 <- data.frame( cbind(Y[,1] , X1) )
 data_2 <- data.frame( cbind(Y[,2] , X2) )
 data_3 <- data.frame( cbind(Y[,3] , X3) )
 data_4 <- data.frame( cbind(Y[,4] , X4) )
 
colnames(data_1) <- c("Y1", "X10", "X11", "X12")
colnames(data_2) <- c("Y2", "X20", "X21", "X22")
colnames(data_3) <- c("Y3", "X30", "X31", "X32")
colnames(data_4) <- c("Y4", "X40", "X41", "X42")

### init values  ##

glm1 <- glm(formula = Y1 ~  X11 + X12 , family = binomial(link = logit),     data = data_1)
glm2 <- glm(formula = Y2 ~  X21 + X22 , family = binomial(link = logit),     data = data_2) 
lm3 <- lm(formula = Y3 ~  X31 + X32  ,   data = data_3)
lm4 <- lm(formula = Y4 ~  X41 + X42 ,    data = data_4)



beta_1<-glm1$coefficients
beta_2<-glm2$coefficients
beta_3<-lm3$coefficients
beta_4<-lm4$coefficients

sig_1 <- 1
sig_2 <- 1
sig_3 <- sd(residuals.lm(lm3))
sig_4 <- sd(residuals.lm(lm4))


rho_12 <-cor( residuals.glm(glm1),  residuals.glm(glm2) )
rho_13 <-cor( residuals.glm(glm1),  residuals.lm(lm3) )
rho_14 <-cor( residuals.glm(glm1),  residuals.lm(lm4) )
rho_23 <-cor( residuals.glm(glm2),  residuals.lm(lm3) )
rho_24 <-cor( residuals.glm(glm2),  residuals.lm(lm4) )
rho_34 <-cor( residuals.lm(lm3) ,   residuals.lm(lm4) )


eta_0 <- c(beta_1,beta_2,beta_3,beta_4, sig_3, sig_4, rho_12, rho_13,rho_14, rho_23, rho_24, rho_34 )

names(eta_0) <- c("beta_10", "beta_11","beta_12",  "beta_20", "beta_21" , "beta_22", "beta_30", "beta_31" , "beta_32" ,"beta_40", "beta_41", "beta_42", "sig_3", "sig_4", "rho_12", "rho_13",
                  "rho_14", "rho_23", "rho_24", "rho_34")



l_12<- Case_2_discrete  (data_1,  data_2, beta_1 , beta_2, rho_12 )
l_13<- Case_1_discrete  (data_3,  data_1, beta_3 , beta_1, sig_3, rho_13 )
l_14<- Case_1_discrete  (data_4,  data_1, beta_4 , beta_1, sig_4, rho_14 )
l_23<- Case_1_discrete  (data_3,  data_2, beta_3 , beta_2, sig_3, rho_23 )
l_24<- Case_1_discrete  (data_4,  data_2, beta_4 , beta_2, sig_4, rho_24 )
l_34<- Case_2_continuous(data_3,  data_4, beta_3 , beta_4, sig_3, sig_4, rho_34 )


l_12_grad <- attr(l_12,"gradient")
l_13_grad <- attr(l_13,"gradient")
l_14_grad <- attr(l_14,"gradient")
l_23_grad <- attr(l_23,"gradient")
l_24_grad <- attr(l_24,"gradient")
l_34_grad <- attr(l_34,"gradient")

l_12_hess <- attr(l_12,"hessian")
l_13_hess <- attr(l_13,"hessian")
l_14_hess <- attr(l_14,"hessian")
l_23_hess <- attr(l_23,"hessian")
l_24_hess <- attr(l_24,"hessian")
l_34_hess <- attr(l_34,"hessian")


parameter_0 <- eta_0 


L_return <- L_return1(data_1, data_2, data_3, data_4, beta_1, beta_2, beta_3, beta_4, sig_3, sig_4, rho_12, rho_13, rho_14, rho_23, rho_24, rho_34)

N_1 <- length(beta_1)
N_2 <- length(beta_2)
N_3 <- length(beta_3)
N_4 <- length(beta_4)

for (i in (1:30))
{
  eta_new <-   eta_0 - solve (attr(L_return,"hessian")) %*% attr(L_return,"gradient")
  diff <- sum (abs( eta_new - eta_0 ))
  print (diff ) 
  eta_0 <- eta_new
  
  # for (i_rho in (11:16))
  #  {
  #   eta_0[i_rho] <- if ( eta_0[i_rho] < 1 &&  eta_0[i_rho] > -1  ) eta_0[i_rho] else  eta_10 [i_rho]
  #  }
  
  beta_1 <- eta_0[1:N_1]
  beta_2 <- eta_0[(N_1+1):(N_1+N_2)]
  beta_3 <- eta_0[(N_1+N_2+1):(N_1+N_2+N_3)]
  beta_4 <- eta_0[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)]
  sig_3 <- eta_0[(N_1+N_2+N_3+N_4+1)]
  sig_4 <- eta_0[(N_1+N_2+N_3+N_4+2)]
  rho_12 <- eta_0[(N_1+N_2+N_3+N_4+3)]
  rho_13 <- eta_0[(N_1+N_2+N_3+N_4+4)]
  rho_14 <- eta_0[(N_1+N_2+N_3+N_4+5)]
  rho_23 <- eta_0[(N_1+N_2+N_3+N_4+6)]
  rho_24 <- eta_0[(N_1+N_2+N_3+N_4+7)]
  rho_34 <- eta_0[(N_1+N_2+N_3+N_4+8)]
  
  L_return <- L_return1(data_1, data_2, data_3, data_4, beta_1, beta_2, beta_3, beta_4, sig_3, sig_4, rho_12, rho_13, rho_14, rho_23, rho_24, rho_34)
}

parameters <- cbind(para_true , parameter_0, eta_new )
colnames(parameters) <- c("True" , "Marginal", "Method")
# parameters



J_theta <-    t ( attr(L_return,"gradient_point") /N)  %*%  attr(L_return,"gradient_point")/N 
H_theta <- - attr(L_return,"hessian")/N

I_Godambe <- solve(H_theta) %*%  J_theta  %*% solve( H_theta )


sd<-sqrt( diag2vec( I_Godambe))
colnames (sd) <-"SD of Method"
#sd
results <- cbind(parameters,sd)
results <- data.frame(results)

results$diff_Marginal <- abs ((results[,1] - results[,2])/results[,1])
results$diff_Method <- abs((results[,1] - results[,3]) / results[,1])
 
#results

beta_1_hat <- eta_new[1:3,1]
beta_2_hat <- eta_new[4:6,1]
beta_3_hat <- eta_new[7:9,1]
beta_4_hat <- eta_new[10:12,1]

y1_hat = X1 %*% beta_1_hat 
y2_hat = X2 %*% beta_2_hat 
y3_hat = X3 %*% beta_3_hat
y4_hat = X4 %*% beta_4_hat

test_1 <- data.frame(cbind( y1,y1_hat ) )
test_2 <- data.frame(cbind( y2,y2_hat ) )
test_3 <- data.frame(cbind( y3,y3_hat ) )
test_4 <- data.frame(cbind( y4,y4_hat ) )

perf_1<-performance(test_1)
perf_2<-performance(test_2)
perf_3<-performance(test_3)
perf_4<-performance(test_4)

perf_method <- c(perf_1,perf_2,perf_3,perf_4)


beta_1_hat <- parameter_0[1:3]
beta_2_hat <- parameter_0[4:6]
beta_3_hat <- parameter_0[7:9]
beta_4_hat <- parameter_0[10:12]

y1_hat = X1 %*% beta_1_hat 
y2_hat = X2 %*% beta_2_hat 
y3_hat = X3 %*% beta_3_hat
y4_hat = X4 %*% beta_4_hat

test_1 <- data.frame(cbind( y1,y1_hat ) )
test_2 <- data.frame(cbind( y2,y2_hat ) )
test_3 <- data.frame(cbind( y3,y3_hat ) )
test_4 <- data.frame(cbind( y4,y4_hat ) )

perf_1<-performance(test_1)
perf_2<-performance(test_2)
perf_3<-performance(test_3)
perf_4<-performance(test_4)

perf_marginal <- c(perf_1,perf_2,perf_3,perf_4)

results
rbind(perf_marginal, perf_method)

