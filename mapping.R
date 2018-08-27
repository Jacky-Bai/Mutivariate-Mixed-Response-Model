

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

L_return <- L_return1(data_1, data_2, data_3, data_4, beta_1, beta_2, beta_3, beta_4, sig_3, sig_4, rho_12, rho_13, rho_14, rho_23, rho_24, rho_34)

 
