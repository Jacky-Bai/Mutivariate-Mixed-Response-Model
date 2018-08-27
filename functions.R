

L_return1 <-function(data_1, data_2, data_3, data_4, beta_1, beta_2, beta_3, beta_4, sig_3, sig_4, rho_12, rho_13, rho_14, rho_23, rho_24, rho_34)
{
  
  N_1 <- length(beta_1)
  N_2 <- length(beta_2)
  N_3 <- length(beta_3)
  N_4 <- length(beta_4)
  
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
  
  l_12_grad_point <- attr(l_12,"grad_point")
  l_13_grad_point <- attr(l_13,"grad_point")
  l_14_grad_point <- attr(l_14,"grad_point")
  l_23_grad_point <- attr(l_23,"grad_point")
  l_24_grad_point <- attr(l_24,"grad_point")
  l_34_grad_point <- attr(l_34,"grad_point")
  
  l_12_hess <- attr(l_12,"hessian")
  l_13_hess <- attr(l_13,"hessian")
  l_14_hess <- attr(l_14,"hessian")
  l_23_hess <- attr(l_23,"hessian")
  l_24_hess <- attr(l_24,"hessian")
  l_34_hess <- attr(l_34,"hessian")
  
  
  
  
  l_beta_1 <-  l_12_grad[1:N_1]                 + l_13_grad[(N_3 + 1):(N_3 + N_1)] + l_14_grad[(N_4 + 1):(N_4 + N_1)]
  l_beta_2 <-  l_12_grad[(N_1 + 1):(N_1 + N_2)] + l_23_grad[(N_3 + 1):(N_3 + N_2)] + l_24_grad[(N_4 + 1):(N_4 + N_2)]
  l_beta_3 <-  l_13_grad[1:N_3]                 + l_23_grad[ 1: N_3 ]              + l_34_grad[ 1: N_3 ]
  l_beta_4 <-  l_14_grad[1:N_4]                 + l_24_grad[1:N_4]                 + l_34_grad[(N_3 + 1):(N_3 + N_4)]
  l_sig_3 <-   l_13_grad[(N_3 + N_1 +1)]        + l_23_grad[(N_3 + N_2 +1)]        + l_34_grad[(N_3 + N_4 +1)]
  l_sig_4  <-  l_14_grad[(N_4 + N_1 +1)]        + l_24_grad[(N_4 + N_2 +1)]        + l_34_grad[(N_3 + N_4 +2)]
  l_rho  <- c(   l_12_grad [ (N_1 + N_2+1) ], l_13_grad[(N_3 + N_1 +2)], l_14_grad[(N_4 + N_1 +2)], l_23_grad[(N_3 + N_2 +2)], l_24_grad[(N_4 + N_2 +2)],l_34_grad[(N_3 + N_4 +3)] )
  
  l_gradient <- c(l_beta_1, l_beta_2, l_beta_3, l_beta_4, l_sig_3, l_sig_4, l_rho )
  names(l_gradient) <- cbind( "Beta_10", "Beta_11", "Beta_12", "Beta_20", "Beta_21", "Beta_22", "Beta_30", "Beta_31", "Beta_32", "Beta_40", "Beta_41", "Beta_42",
                              "sig_3",   "sig_4",   "r_12", "r_13", "r_14","r_23", "r_24", "r_34")   
  l_gradient
  
  
  
  l_beta_1_point <-  l_12_grad_point[,1:N_1]                 + l_13_grad_point[,(N_3 + 1):(N_3 + N_1)] + l_14_grad_point[,(N_4 + 1):(N_4 + N_1)]
  l_beta_2_point <-  l_12_grad_point[,(N_1 + 1):(N_1 + N_2)] + l_23_grad_point[,(N_3 + 1):(N_3 + N_2)] + l_24_grad_point[,(N_4 + 1):(N_4 + N_2)]
  l_beta_3_point <-  l_13_grad_point[,1:N_3]                 + l_23_grad_point[ ,1: N_3 ]              + l_34_grad_point[ ,1: N_3 ]
  l_beta_4_point <-  l_14_grad_point[,1:N_4]                 + l_24_grad_point[,1:N_4]                 + l_34_grad_point[,(N_3 + 1):(N_3 + N_4)]
  l_sig_3_point <-   l_13_grad_point[,(N_3 + N_1 +1)]        + l_23_grad_point[,(N_3 + N_2 +1)]        + l_34_grad_point[,(N_3 + N_4 +1)]
  l_sig_4_point  <-  l_14_grad_point[,(N_4 + N_1 +1)]        + l_24_grad_point[,(N_4 + N_2 +1)]        + l_34_grad_point[,(N_3 + N_4 +2)]
  l_rho_point  <- cbind(   l_12_grad_point [ ,(N_1 + N_2+1) ], l_13_grad_point[,(N_3 + N_1 +2)], l_14_grad_point[,(N_4 + N_1 +2)], l_23_grad_point[,(N_3 + N_2 +2)], l_24_grad_point[,(N_4 + N_2 +2)],l_34_grad_point[,(N_3 + N_4 +3)] )
  
  
  l_gradient_point <- cbind(l_beta_1_point, l_beta_2_point, l_beta_3_point, l_beta_4_point, l_sig_3_point, l_sig_4_point, l_rho_point )
  
  ll_hessian <- matrix(0,length(l_gradient),length(l_gradient))
  
  #####  Hessian
  
  
  
  ll_hessian[1:N_1,1:N_1] <- l_12_hess[1:N_1,1:N_1] + l_13_hess[(N_3+1):(N_3+N_1),(N_3+1):(N_3+N_1)] + l_14_hess[(N_4+1):(N_4+N_1),(N_4+1):(N_4+N_1)]
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+1):(N_1+N_2)] <- l_12_hess[(N_1+1):(N_1+N_2),(N_1+1):(N_1+N_2)] + l_23_hess[(N_3+1):(N_3+N_2),(N_3+1):(N_3+N_2)] + l_24_hess[(N_4+1):(N_4+N_2),(N_4+1):(N_4+N_2)]
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3), (N_1+N_2+1): (N_1+N_2+N_3)] <- l_13_hess[1:N_3,1:N_3] + l_23_hess[1:N_3,1:N_3] + l_34_hess[1:N_3,1:N_3]
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4), (N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_14_hess[1:N_4,1:N_4] + l_24_hess[1:N_4,1:N_4] + l_34_hess[(N_3+1):(N_3+N_4),(N_3+1):(N_3+N_4)]
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+N_4+1)] <- l_13_hess[(N_3+N_1+1),(N_3+N_1+1)] + l_23_hess[(N_3+N_2+1),(N_3+N_2+1)] + l_34_hess[(N_3+N_4+1),(N_3+N_4+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+N_4+2)] <- l_14_hess[(N_4+N_1+1),(N_4+N_1+1)] + l_24_hess[(N_4+N_2+1),(N_4+N_2+1)] + l_34_hess[(N_3+N_4+2),(N_3+N_4+2)]
  
  
  ll_hessian[1:N_1,(N_1+1):(N_1+N_2)] <- l_12_hess[1:N_1,(N_1+1):(N_1+N_2)]
  ll_hessian[(N_1+1):(N_1+N_2),1:N_1] <- t(l_12_hess[1:N_1,(N_1+1):(N_1+N_2)])
  
  ll_hessian[1:N_1,(N_1+N_2+1):(N_1+N_2+N_3)] <- l_13_hess[(N_3+1):(N_3+N_1), 1:N_3]
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),1:N_1] <- t(ll_hessian[1:N_1,(N_1+N_2+1):(N_1+N_2+N_3)])
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_14_hess[(N_4+1):(N_4+N_1), 1:N_4]
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),1:N_1] <- t ( ll_hessian[1:N_1,(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)]  )
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+1)] <-  l_13_hess[(N_3+1):(N_3+N_1),(N_3+N_1+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+1),1:N_1] <-  t(ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+1)])
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+2)] <- l_14_hess[(N_4+1):(N_4+N_1),(N_4+N_1+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),1:N_1] <- t ( ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+2)] )
  
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+3)] <- l_12_hess[1:N_1,(N_1+N_2+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+3),1:N_1] <- t( ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+3)]  )
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+4)] <- l_13_hess[(N_3+1):(N_3+N_1) ,(N_3+N_1+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+4),1:N_1] <- t( ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+4)] )
  
  ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+5)] <- l_14_hess[(N_4+1):(N_4+N_1) ,(N_4+N_1+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+5),1:N_1] <- t( ll_hessian[1:N_1,(N_1+N_2+N_3+N_4+5)] )
  
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+1):(N_1+N_2+N_3)] <- l_23_hess[(N_3+1):(N_3+N_2),1:N_3]
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+1):(N_1+N_2)] <-t ( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+1):(N_1+N_2+N_3)] )  
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_24_hess[(N_4+1):(N_4+N_2),1:N_4]
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+1):(N_1+N_2)] <- t( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] )
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+1)] <- l_23_hess[(N_3+1):(N_3+N_2), (N_3+N_2+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+1):(N_1+N_2)] <- t( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+1)] )
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+2)] <- l_24_hess[(N_4+1):(N_4+N_2), (N_4+N_2+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+1):(N_1+N_2)] <- t( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+2)]  )
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+3)] <- l_12_hess[(N_1+1):(N_1+N_2), (N_1+N_2+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+3),(N_1+1):(N_1+N_2)] <- t ( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+3)] )
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+6)] <- l_23_hess[(N_3+1):(N_3+N_2), (N_3+N_2+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+6),(N_1+1):(N_1+N_2)] <- t ( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+6)]  )
  
  ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+7)] <- l_24_hess[(N_4+1):(N_4+N_2), (N_4+N_2+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+7),(N_1+1):(N_1+N_2)] <- t ( ll_hessian[(N_1+1):(N_1+N_2),(N_1+N_2+N_3+N_4+7)] )
  
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_34_hess[1:N_3,(N_3+1):(N_3+N_4)]
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+1):(N_1+N_2+N_3)] <- t ( ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)]  )
  
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+N_4+1)] <- l_13_hess[1:N_3,(N_3+N_1+1)] + l_23_hess[1:N_3,(N_3+N_2+1)] + l_34_hess[1:N_3,(N_3+N_4+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+1):(N_1+N_2+N_3)] <- l_13_hess[(N_3+N_1+1),1:N_3] + l_23_hess[(N_3+N_2+1),1:N_3] + l_34_hess[(N_3+N_4+1),1:N_3]
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+N_4+2)] <- l_34_hess[1:N_3,(N_3+N_4+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+1):(N_1+N_2+N_3)] <- l_34_hess[(N_3+N_4+2),1:N_3]
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+N_4+4)] <- l_13_hess[1:N_3,(N_3+N_1+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+4),(N_1+N_2+1):(N_1+N_2+N_3)] <- l_13_hess[(N_3+N_1+2),1:N_3]
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),(N_1+N_2+N_3+N_4+6)] <- l_23_hess[1:N_3,(N_3+N_2+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+6),(N_1+N_2+1):(N_1+N_2+N_3)] <- l_23_hess[(N_3+N_2+2),1:N_3]
  
  ll_hessian[(N_1+N_2+1):(N_1+N_2+N_3),16] <- l_34_hess[1:N_3,(N_3+N_4+3)]
  ll_hessian[16,(N_1+N_2+1):(N_1+N_2+N_3)] <- l_34_hess[(N_3+N_4+3),1:N_3]
  
  
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+N_3+N_4+1)] <- l_34_hess[(N_3+1:N_3+N_4),(N_3+N_4+1)]
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_34_hess[(N_3+N_4+1),(N_3+1):(N_3+N_4)]
  
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+N_3+N_4+2)] <- l_14_hess[1:N_4,(N_4+N_1+1)] + l_24_hess[1:N_4,(N_4+N_2+1)] + l_34_hess[(N_3+1):(N_3+N_4),(N_3+N_4+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_14_hess[(N_4+N_1+1),1:N_4] + l_24_hess[(N_4+N_2+1),1:N_4] + l_34_hess[(N_3+N_4+2),(N_3+1):(N_3+N_4)]
  
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+N_3+N_4+5)] <- l_14_hess[1:N_4,(N_1+N_4+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+5),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_14_hess[(N_1+N_4+2),1:N_4]
  
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+N_3+N_4+7)] <- l_24_hess[1:N_4,(N_2+N_4+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+7),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_24_hess[(N_2+N_4+2),1:N_4]
  
  ll_hessian[(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4),(N_1+N_2+N_3+N_4+8)] <- l_34_hess[1:N_3,(N_3+N_4+3)]
  ll_hessian[(N_1+N_2+N_3+N_4+8),(N_1+N_2+N_3+1):(N_1+N_2+N_3+N_4)] <- l_34_hess[(N_3+N_4+3),1:N_3]
  
  
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+N_4+2)] <- l_34_hess[(N_3+N_4+1),(N_3+N_4+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+N_4+1)] <- l_34_hess[(N_3+N_4+2),(N_3+N_4+1)]
  
  
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+N_4+4)] <- l_13_hess[(N_3+N_1+1),(N_3+N_1+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+4),(N_1+N_2+N_3+N_4+1)] <- l_13_hess[(N_3+N_1+2),(N_3+N_1+1)]
  
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+N_4+6)] <- l_23_hess[(N_3+N_2+1),(N_3+N_2+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+6),(N_1+N_2+N_3+N_4+1)] <- l_23_hess[(N_3+N_2+2),(N_3+N_2+1)]
  
  ll_hessian[(N_1+N_2+N_3+N_4+1),(N_1+N_2+N_3+N_4+8)] <- l_34_hess[(N_3+N_4+1),(N_3+N_4+3)]
  ll_hessian[(N_1+N_2+N_3+N_4+8),(N_1+N_2+N_3+N_4+1)] <- l_34_hess[(N_3+N_4+3),(N_3+N_4+1)]
  
  
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+N_4+5)] <- l_14_hess[(N_4+N_1+1),(N_4+N_1+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+5),(N_1+N_2+N_3+N_4+2)] <- l_14_hess[(N_4+N_1+2),(N_4+N_1+1)]
  
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+N_4+7)] <- l_24_hess[(N_4+N_2+1),(N_4+N_2+2)]
  ll_hessian[(N_1+N_2+N_3+N_4+7),(N_1+N_2+N_3+N_4+2)] <- l_24_hess[(N_4+N_2+2),(N_4+N_2+1)]
  
  ll_hessian[(N_1+N_2+N_3+N_4+2),(N_1+N_2+N_3+N_4+8)] <- l_34_hess[(N_3+N_4+2),(N_3+N_4+3)]
  ll_hessian[(N_1+N_2+N_3+N_4+8),(N_1+N_2+N_3+N_4+2)] <- l_34_hess[(N_3+N_4+3),(N_3+N_4+2)]
  
  
  ll_rho  <- diag( c(   l_12_hess[(N_1+N_2+1),(N_1+N_2+1)], l_13_hess[(N_3+N_1+2),(N_3+N_1+2)], l_14_hess[(N_4+N_1+2),(N_3+N_4+2)], l_23_hess[(N_3+N_2+2),(N_3+N_2+2)],
                        l_24_hess[(N_4+N_2+2),(N_4+N_2+2)], l_34_hess[(N_3+N_4+3),(N_3+N_4+3)] ),6,6)
  
  
  
  colnames(ll_hessian) <- c( names(l_gradient))
  rownames(ll_hessian) <- c( names(l_gradient))
  ll_hessian[(N_1+N_2+N_3+N_4+3):(N_1+N_2+N_3+N_4+8),(N_1+N_2+N_3+N_4+3):(N_1+N_2+N_3+N_4+8)] <-  ll_rho
  
  
  l_gradient
  ll_hessian
  l_gradient_point
  L_return <- 0 
  attr(L_return,"gradient")<- l_gradient
  attr(L_return,"hessian")<- ll_hessian
  attr(L_return,"gradient_point")<- l_gradient_point
  return(L_return)
}


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


performance <- function(df){
  df[3] <- mean(df[,1])
  names(df)[1] <- "y"
  names(df)[2] <- "y_hat"
  names(df)[3] <- "y_bar"
  attach(df)
  R_squared <- 1 - sum((y - y_hat)^2) /sum ( (y - y_bar)^2 )
  RMSE <- sqrt( sum ( (y-y_hat)^2 )/ nrow(data_1) )
  MAE <- sum(abs( y-y_hat )) / nrow(data_1)
  Performance <- c (R_squared,RMSE, MAE )
  names(Performance) <- c("R_squared","RMSE", "MAE")
  detach(df)
  return (Performance )
}




detach(data_1)
data_1 <- data.frame(cbind( xx,zz ) )
data_1[3] <- mean(data1[,1])
names(data_1)[1] <- "y"
names(data_1)[2] <- "y_hat"
names(data_1)[3] <- "y_bar"
attach(data_1)
R_squared <- 1 - sum((y - y_hat)^2) /sum ( (y - y_bar)^2 )
RMSE <- sqrt( sum ( (y-y_hat)^2 )/ nrow(data_1) )
MAE <- sum(abs( y-y_hat )) / nrow(data_1)
Performance <- c (R_squared,RMSE, MAE )
print (Performance )




theta=c(30,40,22,33,40)
funk=function(param,theta){
  R=param[1]
  DP=param[2]
  fit=sum( log (sqrt(1-R^2)* exp( -1/(2*R^2) * (( sqrt(1-R^2) * 1 /(( pnorm(theta))^2) -DP ))^2+ 1/2 * 1
                                  /(( pnorm(theta))^2) )) )
  return(fit)
}
optim(par=c(0.1 , 0.1), fn=funk, theta=theta)
