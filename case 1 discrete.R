
data_k <- data_3
data_j <- data_1
beta_k <- beta_3
beta_j <- beta_1
sig_k <- sig_3
rho_kj <- rho_13


Case_1_discrete <- function(data_k, data_j, beta_k, beta_j, sig_k, rho_kj) {
  
  rho <- rho_kj
  N_k <- length(beta_k)
  N_j <- length(beta_j)
  
  z_j<- data_j[,1]
  z_k<- data_k[,1]
  X_j<- as.matrix (data_j[,-1])
  X_k<- as.matrix (data_k[,-1])
  
  s_k<-(z_k- X_k %*% beta_k )/sig_k
  o_j <- -( X_j %*% beta_j + rho * s_k )/(  sqrt(1-rho^2 ))
  phi_oj <- dnorm (o_j)
  
  Phi_oj <- rep(0, nrow(o_j))
  Phi_oj_2 <- rep(0, nrow(o_j))
  for ( val in c(1:nrow(o_j) ) )   {
    Phi_oj[val] <- if (z_j[val] == 0 )  1/pnorm(o_j[val]) else - 1/ ( 1- pnorm(o_j[val])) 
    Phi_oj_2[val] <- if (z_j[val] == 0 )  - 1/(pnorm(o_j[val]) **2) else  - 1/(( 1- pnorm(o_j[val]))**2 )  
  }
  
  
  d_s_beta_k <- - 1/sig_k * X_k
  d_s_sig_k <- - 1/sig_k * s_k 
  
  d_o_beta_k <-  X_k * rho/(sig_k * sqrt (1- rho^2))   
  d_o_beta_j <-  - X_j / (  sqrt (1- rho^2))
  d_o_sig_k <-   rho * s_k / (sig_k * sqrt (1- rho^2))    
  d_o_d_rho <-  -s_k/sqrt (1- rho^2) - X_j %*% beta_j * rho / (  (1- rho^2)^(3/2) ) - rho^2 * s_k /(1- rho^2)^(3/2) 
  
  d2_s_sig_k_2 <-  2*s_k /sig_k^2
  d2_s_beta_k_sig_k <-  1/sig_k^2 * X_k
  
  d2_o_rho_sig_k <- s_k/sig_k * (  1/ sqrt(1-rho^2 ) + rho^2 * (1- rho^2)^(-3/2)   )
  d2_o_beta_j_rho <- -  X_j * rho /  (   (1- rho^2)^(3/2) )
  d2_o_beta_k_rho <-  X_k /sig_k * (  1/ sqrt(1-rho^2 ) + rho^2 * (1- rho^2)^(-3/2) )
  d2_o_sig_k_2 <- - 2 * rho * s_k  / ( sig_k^2 * sqrt(1-rho^2 )  )
  d2_o_beta_k_sig_k <- - X_k * rho / ( sig_k^2 * sqrt(1-rho^2 )  ) 
  d2_o_rho_2 <- - rho * s_k /((1-rho^2)^(3/2)) -  X_j %*% beta_j /   (  1/ ((1-rho^2)^(3/2)) + 3* rho ^2 / ((1-rho^2)^(5/2)) ) - 
                 s_k * (2*rho/ ((1-rho^2)^(3/2)) + 3* rho ^3 / ((1-rho^2)^(5/2)) )
  
  
  d_sig_k <- matrix( 0,nrow=nrow(o_j),ncol= (N_k+N_j+2))
  d_s <- matrix( 0,nrow=nrow(o_j),ncol= (N_k+N_j+2))
  d_o <- matrix( 0,nrow=nrow(o_j),ncol= (N_k+N_j+2))
  
  d_sig_k[,(N_k+N_j+1)]<- 1
  d_s[,1:N_k] <- d_s_beta_k
  d_s[,(N_k+N_j+1)] <- d_s_sig_k
  d_o[,1:N_k] <- d_o_beta_k
  d_o[,(N_k + 1):( N_k + N_j )] <- d_o_beta_j
  d_o[,(N_k+N_j+1)] <- d_o_sig_k
  d_o[,(N_k+N_j+2)] <- d_o_d_rho  
  
  

  s_k_M <- rep( s_k, (N_k+N_j+2))
  Phi_oj_M <- rep (Phi_oj, (N_k+N_j+2))
  phi_oj_M <- rep (phi_oj, (N_k+N_j+2))
  
  d_l_d_eta_N <- -1/sig_k * d_sig_k - s_k_M * d_s + Phi_oj_M * phi_oj_M * d_o
 # colnames(d_l_d_eta_N) <- c("beta_k0","beta_k1","beta_j0","beta_j1","sig_k" ,"rho" )
  d_l_d_eta<-colSums( d_l_d_eta_N  )
  d_l_d_eta
  
  
  # Hessian
  
  d2_l_eta_etat <- matrix(0,(N_k+N_j+2),(N_k+N_j+2))
  
  d2_s_M <- array(0, dim=c( (N_k+N_j+2), (N_k+N_j+2), N))
  for ( ids in c(1:N))
  {
    d2_s_M[(N_k+N_j+1),(N_k+N_j+1), ids] <-  d2_s_sig_k_2[ids]
    d2_s_M[1:N_k,(N_k+N_j+1), ids ] <- d2_s_beta_k_sig_k[ids,]
    d2_s_M[(N_k+N_j+1),1:N_k, ids] <- d2_s_beta_k_sig_k[ids,]
  }
  d2_o_M <- array(0, dim=c((N_k+N_j+2), (N_k+N_j+2), N ))
  for ( ids in c(1:N))
  {
    d2_o_M[(N_k+N_j+2),(N_k+N_j+2), ids ] <- d2_o_rho_2[ids]
    d2_o_M[(N_k+N_j+1),(N_k+N_j+1), ids ] <-  d2_o_sig_k_2[ids]
    d2_o_M[1:N_k,(N_k+N_j+1), ids ] <- d2_o_beta_k_sig_k[ids,]
    d2_o_M[(N_k+N_j+1),1:N_k, ids ] <- d2_o_beta_k_sig_k[ids,]
    d2_o_M[1:N_k,(N_k+N_j+2),ids ] <- d2_o_beta_k_rho[ids,]
    d2_o_M[(N_k+N_j+2),1:N_k,ids]<- d2_o_beta_k_rho[ids,]
    d2_o_M[(N_k+1):(N_k+N_j),(N_k+N_j+2),ids]<- d2_o_beta_j_rho[ids,]
    d2_o_M[(N_k+N_j+2),(N_k+1):(N_k+N_j),ids]<- d2_o_beta_j_rho[ids,]
    
    d2_o_M[(N_k+N_j+2),(N_k+N_j+1),ids]<- d2_o_rho_sig_k[ids]
    d2_o_M[(N_k+N_j+1),(N_k+N_j+2),ids]<- d2_o_rho_sig_k[ids]
  }
  
  

  
  
  
  for ( rows in c(1: (N_k+N_j+2)  ))   {
    for ( cols in c(1:(N_k+N_j+2))) {
      d2_l_eta_etat [rows, cols] <- sum(  1/sig_k^2 * d_sig_k[, rows] * d_sig_k[,cols] - d_s[,rows] * d_s[,cols] - s_k * d2_s_M[rows,cols, ]
                                          +  + Phi_oj_2 * phi_oj^2  * d_o[,rows] * d_o[,cols]
                                          +  Phi_oj * ( - o_j * phi_oj * d_o[,rows] * d_o[,cols] + phi_oj * d2_o_M [rows,cols, ] )
                                          )
                                             
      
    }
  }
  
  l_jk <-0
  
  attr(l_jk,"gradient")<- d_l_d_eta
  attr(l_jk,"hessian")<- d2_l_eta_etat
  attr(l_jk,"grad_point")<-  d_l_d_eta_N
  l_jk
  
  
  return(l_jk)
}


glm1 <- glm(formula = Y1 ~  X11 + X12 , family = binomial(link = logit),     data = data_1)
lm3 <- lm(formula = Y3 ~  X31 + X32  ,   data = data_3)

beta_3<-lm3$coefficients
beta_1<-glm1$coefficients

rho_13 <-cor(residuals.lm(lm3) , residuals.glm(glm1))
sig_3<- sd(residuals.lm(lm3))
sig_1<- sd(residuals.glm(glm1))


parameter_0<- c(beta_k, beta_j, sig_k, rho )

Case_2<- Case_1_discrete (data_3, data_1, beta_3, beta_1, sig_3, rho_13) 


# Newton's test  


eta_0 <- parameter_0  
l_return <-  Case_1_discrete (data_3, data_1, beta_3, beta_1, sig_3, rho_13) 
for (i in (1:20))
{
  eta_new <-   eta_0 - solve (attr(l_return,"hessian")) %*% attr(l_return,"gradient")
  diff <- sum (abs( eta_new - eta_0 ))
  print (diff ) 
  eta_0 <- eta_new  
  
  beta_3 <- eta_0[1:N_1]
  beta_1 <- eta_0[(N_1+1): (N_1 + N_2)]
  sig_3  <- eta_0[(N_1 + N_2 + 1)]
  rho_13 <- eta_0[(N_1 + N_2 +2)]   
  
  l_return <-  Case_1_discrete (data_3, data_1, beta_3, beta_1, sig_3, rho_13) 
}

para_true <- c(t(Beta[3,]) , t(Beta[1,]) , 1,  Rho[1,3] )
parameters <- cbind(para_true , parameter_0, eta_new )
colnames(parameters) <- c("True" , "Marginal", "Method")
parameters
