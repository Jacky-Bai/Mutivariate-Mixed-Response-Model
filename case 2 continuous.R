
Case_2_continuous<- function(data_3, data_4, beta_3, beta_4, sig_3, sig_4, rho_34) {
  

  rho <- rho_34
  X_3 <- as.matrix (data_3 [,-1] )
  X_4 <- as.matrix (data_4 [,-1] )
  Y3<-data_3 [,1]
  Y4<-data_4 [,1]

  N_3 <- length(beta_3)
  N_4 <- length(beta_4)
    
  mu_3 = X_3 %*% beta_3
  mu_4 = X_4 %*% beta_4
  
  
  eta = 1/(1-rho^2) *( (Y3 - mu_3 )^2/sig_3^2 + ( Y4 - mu_4 )^2/sig_4^2  - 2 * rho*(Y3 - mu_3) * (Y4 - mu_4)/(sig_3 *sig_4)  )
  
  
  l_34.beta3<- t(X3) %*% (1/(1-rho^2)*( ( Y3 - mu_3 )/sig_3^2 - rho * (Y4 - mu_4)/(sig_3 *sig_4) ))
  l_34.beta4<- t(X4) %*% (1/(1-rho^2)*( ( Y4 - mu_4 )/sig_4^2 - rho * (Y3 - mu_3)/(sig_3 *sig_4) ))
  l_34.sig_3 = sum(-1/(sig_3)*( 1- 1/(1-rho^2)*(( Y3 - mu_3)^2/sig_3^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )) )
  l_34.sig_4 = sum(-1/(sig_4)*( 1- 1/(1-rho^2)*(( Y4 - mu_4)^2/sig_4^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )) )
  l_34.rho = sum( rho/(1-rho^2) * (1 - eta + (Y3 - mu_3)*(Y4 - mu_4)/(rho*sig_3 *sig_4) ) )
  
  l_34.beta3_point<- X3 * rep ( (1/(1-rho^2)*( ( Y3 - mu_3 )/sig_3^2 - rho * (Y4 - mu_4)/(sig_3 *sig_4) )), N_3)
  l_34.beta4_point<- X4 * rep ( (1/(1-rho^2)*( ( Y4 - mu_4 )/sig_4^2 - rho * (Y3 - mu_3)/(sig_3 *sig_4) )), N_4)
  l_34.sig_3_point = -1/(sig_3)*( 1- 1/(1-rho^2)*(( Y3 - mu_3)^2/sig_3^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )) 
  l_34.sig_4_point = -1/(sig_4)*( 1- 1/(1-rho^2)*(( Y4 - mu_4)^2/sig_4^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )) 
  l_34.rho_point =  rho/(1-rho^2) * (1 - eta + (Y3 - mu_3)*(Y4 - mu_4)/(rho*sig_3 *sig_4) ) 
  
  
  # Hessian  ;
  
  
  l_34.beta_3beta_3<--1/((1-rho^2)*sig_3^2) * t(X3) %*% X3
  l_34.beta_4beta_4<--1/((1-rho^2)*sig_4^2) * t(X4) %*% X4
  l_34.beta_3beta_4<- rho/((1-rho^2)*sig_3* sig_4) * t(X3) %*% X4
  
  l_34.sig_3sig_3 <-  sum( 1/(sig_3^2) * ( ( 1- 1/(1-rho^2)*(( Y3-mu_3)^2/sig_3^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) ) ) +
                                             1/(1-rho^2)*( - 2*(Y3-mu_3)^2 /sig_3^2 + rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )  )   ) 
  
  
  l_34.sig_4sig_4 <- sum( 1/(sig_4^2)* ( ( 1- 1/(1-rho^2)*(( Y4 - mu_4)^2/sig_4^2 - rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )) +
                                           1/(1-rho^2)*( - 2*(Y4-mu_4)^2 /sig_4^2 + rho*(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) )  )   )
  
  
  l_34.rho_rho <- sum ( (1/(1-rho^2) + 2*rho^2/(1-rho^2)^2) * (1 - eta + (Y3 - mu_3)*(Y4 - mu_4)/(rho*sig_3 *sig_4) )  -
                          rho/(1-rho^2)* (  2/(1-rho^2)* ( rho * eta - (Y3 - mu_3)*(Y4 - mu_4)/(sig_3 *sig_4) ) +
                                              (Y3 - mu_3)*(Y4 - mu_4)/(rho^2 * sig_3 *sig_4) ) )
  
  
  
  l_34.beta_3sig_3<- t(X3) %*% (-1/(1-rho^2) * ( 2*(Y3-mu_3)/sig_3^3 - rho * (Y4-mu_4)/(sig_3^2 *sig_4) ) )
  l_34.beta_3sig_4<- t(X3) %*% ( 1/(1-rho^2) * (  rho * (Y4-mu_4)/(sig_3 *sig_4^2) ) )
  l_34.beta_3_rho<-  t(X3) %*% ( ( 2*rho/(1-rho^2)^2 * ( (Y3-mu_3)/sig_3^2 - rho * (Y4-mu_4)/(sig_3 *sig_4) ) ) - 1/(1-rho^2) * (Y4-mu_4)/(sig_3 *sig_4))
  
  l_34.beta_4sig_3<- t(X4) %*% ( 1/(1-rho^2) * (  rho * (Y3-mu_3)/(sig_3^2 *sig_4) ) )
  l_34.beta_4sig_4<- t(X4) %*% (-1/(1-rho^2) * ( 2*(Y4-mu_4)/sig_4^3 - rho * (Y3-mu_3)/(sig_3 *sig_4^2) ) )
  l_34.beta_4_rho<-  t(X4) %*% ( ( 2*rho/(1-rho^2)^2 * ( (Y4-mu_4)/sig_4^2 - rho * (Y3-mu_3)/(sig_3 *sig_4) ) ) - 1/(1-rho^2) * (Y3-mu_3)/(sig_3 *sig_4))
  
  
  l_34.sig_3sig_4 <- sum ( 1/(1-rho^2) * rho *(Y3 - mu_3)*(Y4 - mu_4)/(sig_3^2 * sig_4^2) )
  
  
  l_34.sig_3_rho <- sum(  2 *rho /( (1-rho^2)^2 * sig_3) * ( (Y3-mu_3)^2/sig_3^2 - rho * (Y3-mu_3) * (Y4-mu_4)/(sig_3 *sig_4) )-
                            1 /( (1-rho^2) * sig_3) *(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 * sig_4)  )
  
  
  l_34.sig_4_rho <- sum(  2 *rho /( (1-rho^2)^2 * sig_4) * ( (Y4-mu_4)^2/sig_4^2 - rho * (Y3-mu_3) * (Y4-mu_4)/(sig_3 *sig_4) )-
                            1 /( (1-rho^2) * sig_4) *(Y3 - mu_3)*(Y4 - mu_4)/(sig_3 * sig_4)  )
  
  
  hessian<-matrix(0, ( N_3 + N_4 +3) ,(N_3 + N_4 +3) )
  
  
  hessian[1:N_3,1:N_3] <-l_34.beta_3beta_3
  hessian[(N_3 + 1) : (N_3 + N_4) , (N_3 + 1) : (N_3 + N_4)  ] <-l_34.beta_4beta_4
  
 
  hessian[(N_3 + N_4 +1) , (N_3 + N_4 +1)] <-l_34.sig_3sig_3
  hessian[(N_3 + N_4 +2) , (N_3 + N_4 +2)] <-l_34.sig_4sig_4
  hessian[(N_3 + N_4 +3), (N_3 + N_4 +3)] <-l_34.rho_rho
  
  hessian[1:N_3, (N_3 +1):(N_3 + N_4)] <-l_34.beta_3beta_4
  hessian[(N_3 +1):(N_3 + N_4) , 1:N_3] <-t(l_34.beta_3beta_4)
  
  hessian[1:N_3,(N_3 + N_4 +1)] <- l_34.beta_3sig_3
  hessian[(N_3 + N_4 +1),1:N_3] <- t(l_34.beta_3sig_3 )
  hessian[1:N_3,(N_3 + N_4 +2)] <- l_34.beta_3sig_4
  hessian[(N_3 + N_4 +2),1:N_3] <- t(l_34.beta_3sig_4)
  hessian[1:N_3,(N_3 + N_4 +3)] <- l_34.beta_3_rho
  hessian[(N_3 + N_4 +3),1:N_3] <- t(l_34.beta_3_rho)
  
  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +1)] <- l_34.beta_4sig_3
  hessian[(N_3 + N_4 +1),(N_3 +1):(N_3 + N_4 )] <- t(l_34.beta_4sig_3 )
  
  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +2)] <- l_34.beta_4sig_4
  hessian[(N_3 + N_4 +2),(N_3 +1):(N_3 + N_4) ] <- t(l_34.beta_4sig_4)
  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +3)] <- l_34.beta_4_rho
  hessian[(N_3 + N_4 +3),(N_3 +1):(N_3 + N_4) ] <- t(l_34.beta_4_rho)
  
  hessian[(N_3 + N_4 +1),(N_3 + N_4 +2)] <- l_34.sig_3sig_4
  hessian[(N_3 + N_4 +2),(N_3 + N_4 +1)] <- l_34.sig_3sig_4
  hessian[(N_3 + N_4 +1),(N_3 + N_4 +3)] <- l_34.sig_3_rho
  hessian[(N_3 + N_4 +3),(N_3 + N_4 +1)] <- l_34.sig_3_rho
  hessian[(N_3 + N_4 +2),(N_3 + N_4 +3)] <- l_34.sig_4_rho
  hessian[(N_3 + N_4 +3),(N_3 + N_4 +2)] <- l_34.sig_4_rho
  
  
  grad<- c(t(l_34.beta3), t(l_34.beta4), l_34.sig_3, l_34.sig_4 , l_34.rho)
  grad_point <- cbind(l_34.beta3_point, l_34.beta4_point, l_34.sig_3_point, l_34.sig_4_point , l_34.rho_point)
  
  l_34<-0
  attr(l_34,"gradient")<- grad
  attr(l_34,"hessian")<- hessian
  attr(l_34,"grad_point")<-grad_point
  
  return(l_34)
  }



lm3 <- lm(formula = Y3 ~  X31 + X32  ,   data = data_3)
lm4 <- lm(formula = Y4 ~  X41 + X42 ,    data = data_4)

beta_3<-lm3$coefficients
beta_4<-lm4$coefficients


rho <-cor(residuals.lm(lm3) , residuals.lm(lm4))
sig_3<- sd(residuals.lm(lm3))
sig_4<- sd(residuals.lm(lm4))

parameter_0<- c(beta_3, beta_4, sig_3, sig_4, rho )
 
Case_1 <- Case_2_continuous(data_3, data_4, beta_3, beta_4, sig_3, sig_4, rho_34)


# Newton's test  


eta_0 <- c(beta_3, beta_4, sig_3, sig_4, rho_34)
l_return <- Case_2_continuous(data_3, data_4, beta_3, beta_4, sig_3, sig_4, rho_34)
for (i in (1:100))
{
  eta_new <-   eta_0 - solve (attr(l_return,"hessian")) %*% attr(l_return,"gradient")
  diff <- sum (abs( eta_new - eta_0 ))
  print (diff ) 
  eta_0 <- eta_new
  beta_3 <- eta_0[1:N_3]
  beta_4 <- eta_0[(N_3+1): (N_3 + N_4)]
  sig_3 <- eta_0[N_3 + N_4+1]   
  sig_4 <- eta_0[N_3 + N_4+2]   
  rho_34 <- eta_0[N_3 + N_4+3]   
  l_return <- Case_2_continuous(data_3, data_4, beta_3, beta_4, sig_3, sig_4, rho_34)
}

para_true <- c(t(Beta[3,]) , t(Beta[4,]) , 1, 1, Rho[3,4] )
parameters <- cbind(para_true , parameter_0, eta_new )
colnames(parameters) <- c("True" , "Marginal", "Method")
parameters




