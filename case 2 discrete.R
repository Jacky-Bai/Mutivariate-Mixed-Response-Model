
Case_2_discrete <- function(data_1, data_2, beta_1, beta_2, rho_12 ) {
  
  N_1 <- length(beta_1)
  N_2 <- length(beta_2)
  
  rho <- rho_12
  
  z_1<- data_1[,1]
  z_2<- data_2[,1]
  X_1<- as.matrix (data_1[,-1])
  X_2<- as.matrix (data_2[,-1])
  
  mu_1 <- X_1 %*% beta_1
  mu_2 <- X_2 %*% beta_2
  q_1 <- 2 * z_1 - 1 
  q_2 <- 2 * z_2 - 1 
  rho_tulda <- q_1 * q_2 * rho
  w_1 <- q_1 * mu_1
  w_2 <- q_2 * mu_2
  
  Q_1 <- w_1/sig_1
  Q_2 <- w_2/sig_2
  delta <- 1/sqrt(1- rho_tulda^2)
  V_1 <- delta * (Q_2 - rho_tulda * Q_1)
  V_2 <- delta * (Q_1 - rho_tulda * Q_2)
  G_1 <- dnorm(Q_1) * pnorm(V_1)
  G_2 <- dnorm(Q_2) * pnorm(V_2)
  sigma_tulda <- array(0, dim=c(N, 2, 2)) 
  
  phi_2 <- rep(0,N)
  Phi_2 <- rep(0,N)
  for (ids in c(1:N))
  {
    sigma_tulda[ids,,] <- matrix (c(1,rho_tulda[ids] ,rho_tulda[[ids]], 1 ),nrow=2,  ncol=2,  byrow = TRUE)  
    phi_2[ids] <- dmvnorm( c(Q_1[ids] , Q_2[ids] ), mean = rep(0, 2), sigma = sigma_tulda[ids,,] )
    Phi_2[ids] <- pmvnorm(lower=-Inf, upper= c(Q_1[ids] , Q_2[ids] ) , mean=rep(0, 2),  sigma=sigma_tulda[ids,,] )
  }
  
  A_beta_1 <- - ( Q_1 * G_1 / Phi_2 + rho_tulda * phi_2/Phi_2 + G_1^2 / Phi_2^2   ) / sig_1 ^2  
  A_beta_2 <- - ( Q_2 * G_2 / Phi_2 + rho_tulda * phi_2/Phi_2 + G_2^2 / Phi_2^2   ) / sig_2 ^2
  A_beta_12 <- ( phi_2/Phi_2 - G_1 * G_2 / Phi_2^2 )/(sig_1 * sig_2 )
  
  Arr_beta_1 <- array(0, dim = c(N_1,N_1,N))
  Arr_beta_2 <- array(0, dim = c(N_2,N_2,N))
  Arr_beta_12 <- array(0, dim = c(N_1,N_2,N))
  
  
  d_l_d_beta_1 <- t(X1) %*% ( q_1/sig_1 * G_1 /Phi_2 )
  d_l_d_beta_2 <- t(X2) %*% ( q_2/sig_2 * G_2 /Phi_2 ) 
  d_l_d_rho <- sum (q_1 * q_2 * phi_2 / Phi_2 )
  
  d_l_d_beta_1_point <- X1 * rep ( q_1/sig_1 * G_1 /Phi_2 , N_1 )
  d_l_d_beta_2_point <- X2 * rep ( q_2/sig_2 * G_2 /Phi_2 , N_2 ) 
  d_l_d_rho_point <- q_1 * q_2 * phi_2 / Phi_2 
  
  
  for (ids in c(1:N))
  {
    Arr_beta_1[,,ids] <- X_1[ids,] %*% t(X_1[ids,]) * A_beta_1 [ids]
    Arr_beta_2[,,ids] <- X_2[ids,] %*% t(X_2[ids,]) * A_beta_2 [ids]
    Arr_beta_12[,,ids] <- X_2[ids,] %*% t(X_1[ids,]) * A_beta_12 [ids]
  }
  
  d2_l_d_beta1_beta1 <- matrix( c(sum(Arr_beta_1[1,1, ]) , sum(Arr_beta_1[1,2, ]), sum(Arr_beta_1[1,3, ]), 
                                  sum(Arr_beta_1[2,1, ]), sum(Arr_beta_1[2,2, ]), sum(Arr_beta_1[2,3, ]),
                                  sum(Arr_beta_1[3,1, ]), sum(Arr_beta_1[3,2, ]), sum(Arr_beta_1[3,3, ]) ), 3, 3)
  d2_l_d_beta2_beta2 <- matrix( c(sum(Arr_beta_2[1,1, ]) , sum(Arr_beta_2[1,2, ]),sum(Arr_beta_2[1,3, ]),
                                  sum(Arr_beta_2[2,1, ]), sum(Arr_beta_2[2,2, ]) , sum(Arr_beta_2[2,3, ]),
                                  sum(Arr_beta_2[3,1, ]), sum(Arr_beta_2[3,2, ]) , sum(Arr_beta_2[3,3, ]) ), 3, 3)
  d2_l_d_beta1_beta2 <- matrix( c(sum(Arr_beta_12[1,1, ]) , sum(Arr_beta_12[1,2, ]), sum(Arr_beta_12[1,3, ]) ,
                                  sum(Arr_beta_12[2,1, ]), sum(Arr_beta_12[2,2, ]), sum(Arr_beta_12[2,3, ]) ,
                                  sum(Arr_beta_12[3,1, ]), sum(Arr_beta_12[3,2, ]), sum(Arr_beta_12[3,3, ]) ), 3, 3)
    
  d2_l_d_beta1_rho <-   t(X1) %*%  ( q_2 /sig_1 * phi_2/Phi_2* (rho_tulda * delta * V_1 - Q_1 -  G_1 / Phi_2))
  d2_l_d_beta2_rho <-   t(X2) %*%  ( q_1 /sig_2 * phi_2/Phi_2* (rho_tulda * delta * V_2 - Q_2 -  G_2 / Phi_2))
  d2_l_d_rho_rho <- sum( phi_2/Phi_2 * ( delta^2 * rho_tulda * ( 1 - delta^2 * (Q_1^2 + Q_2^2 - 2*rho_tulda* Q_1 * Q_2 )) + delta^2*Q_1*Q_2 -phi_2/Phi_2 ))
  
  
  
  grad<- c(t(d_l_d_beta_1), t(d_l_d_beta_2), d_l_d_rho )
  grad_point<- cbind(d_l_d_beta_1_point, d_l_d_beta_2_point, d_l_d_rho_point )
  
  hessian<-matrix(0,(N_1 + N_2 + 1),(N_1 + N_2 + 1))
  
  
  hessian [1:N_1,1:N_1] <- d2_l_d_beta1_beta1
  hessian [(N_1+1):(N_1+N_2),(N_1+1):(N_1+N_2)] <- d2_l_d_beta2_beta2
  hessian [1:N_1,(N_1+1):(N_1+N_2)] <- d2_l_d_beta1_beta2
  hessian [(N_1+1):(N_1+N_2),1:N_1] <- d2_l_d_beta1_beta2
  hessian [(N_1+N_2+1),(N_1+N_2+1)] <- d2_l_d_rho_rho
  hessian [(1:N_1),(N_1+N_2+1)] <- d2_l_d_beta1_rho
  hessian [(N_1+N_2+1),1:N_1] <- d2_l_d_beta1_rho
  hessian [(N_1+1):(N_1+N_2),(N_1+N_2+1)] <- d2_l_d_beta2_rho
  hessian [(N_1+N_2+1),(N_1+1):(N_1+N_2)] <- d2_l_d_beta2_rho
  
  
  
  l_12<-0
  attr(l_12,"gradient")<- grad
  attr(l_12,"hessian")<- hessian
  attr(l_12,"grad_point")<- grad_point
  l_12
  
  return(l_12)
}



glm1 <- glm(formula = Y1 ~  X11 + X12 , family = binomial(link = logit),     data = data_1)
glm2 <- glm(formula = Y2 ~  X21 + X22 , family = binomial(link = logit),     data = data_2) 

beta_1<-glm1$coefficients
beta_2<-glm2$coefficients


rho <-cor(residuals.glm(glm1) , residuals.glm(glm2))
sig_1<- sd(residuals.glm(glm1))
sig_2<- sd(residuals.glm(glm2))


parameter_0<- c(beta_1, beta_2,  rho )

Y<-cbind(data_1$Y1,data_1$Y2   )

Case_3<- Case_2_discrete(data_1, data_2, beta_1, beta_2, rho_12 )


# Newton's test  


eta_0 <- parameter_0  
l_return <-  Case_2_discrete(data_1, data_2, beta_1, beta_2, rho_12 )
for (i in (1:20))
{
  eta_new <-   eta_0 - solve (attr(l_return,"hessian")) %*% attr(l_return,"gradient")
  diff <- sum (abs( eta_new - eta_0 ))
  print (diff ) 
  eta_0 <- eta_new
  
  beta_1 <- eta_0[1:N_1]
  beta_2 <- eta_0[(N_1+1): (N_1 + N_2)]
  rho_12 <- eta_0[N_1 + N_2 +1]   
  
  l_return <- Case_2_discrete(data_1, data_2, beta_1, beta_2, rho_12 )
}

para_true <- c(t(Beta[1,]) , t(Beta[2,]) ,  Rho[1,2] )
parameters <- cbind(para_true , parameter_0, eta_new )
colnames(parameters) <- c("True" , "Marginal", "Method")
parameters


