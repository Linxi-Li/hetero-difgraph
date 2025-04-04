
############################# Functions for main algorithms ############################
# Initialize function (Singular Value Decomposition method)
Initialize_fuc_svd=function(data_0,lambda_1,r){
  Sigma_hat0=cov(data_0)
  p=nrow(Sigma_hat0)
  
  Theta_hat0 = solve(Sigma_hat0)
  S_0=Mcp_prox(Theta_hat0,lambda_1)
  SVD=eigen(S_0-Theta_hat0)
  D_r=diag(sqrt(abs(SVD$values)))[,1:r]
  U_0=(SVD$vectors)%*%D_r
  L_0=U_0%*%t(U_0)
  Theta_0=S_0-L_0
  return(list(S_0=S_0,U_0=U_0,L_0=L_0,Theta_0=Theta_0))
}

# Gradient calculation function for baseline group
Grad_0=function(Sigma_hat0,S_0,L_0,grad_1){
  p=ncol(Sigma_hat0)
  grad_0=Sigma_hat0-solve(S_0-L_0)-grad_1
  return(grad_0)
}
# Gradient calculation function for target heterogeneous group
Grad_1=function(data_K,estimate_K){
  n1=nrow(data_K)
  p=ncol(data_K)
  
  prob=estimate_K$prob
  Mu=estimate_K$Mu
  Theta=estimate_K$Theta
  K=estimate_K$K
  
  f.mat = matrix(0,n1,K)
  G.mat = matrix(0,n1,K)
  grad_1= matrix(0,p,p)
  
  
  for(k.ind in 1:K){
    f.mat[,k.ind]=f_den_vec(data_K, as.numeric(Mu[k.ind,]), Theta[,,k.ind])
  }    
  for(k.ind in 1:K){
    for(i in 1:n1) {
      G.mat[i,k.ind] =prob[k.ind]*f.mat[i,k.ind]/(prob%*%f.mat[i,])
    }
  }
  G.mat[is.na(G.mat)] = 0
  for(k.ind in 1:K){
    for(i in 1:n1) {
      grad_1=grad_1+G.mat[i,k.ind]*(solve(Theta[,,k.ind]+diag(1e-5,p))-(data_K[i,] - as.numeric(Mu[k.ind,]))%*% t(data_K[i,] - as.numeric(Mu[k.ind,])))
    }
  }
  return(2*grad_1/n1)
}

# Alternating gradient descent update function
Update_ATGD=function(data_0,S_0,U_0,eta1=0.1,eta2=0.1,lambda_1,r,eps=1e-2,grad_1){
  
  Sigma_hat0=cov(data_0)
  S_old=S_0
  U_old=U_0
  L_old=U_old%*%t(U_old)
  iter=0
  err_S=err_L=1
  diff_0=1
  
  while ((iter<=500)&(diff_0 >eps)){
    iter=iter+1

    S_new= S_old-eta1*Grad_0(Sigma_hat0,S_old,L_old,grad_1)
    S_new=Mcp_prox(S_new,lambda_1)
    
    U_new= U_old-eta2*(-2*Grad_0(Sigma_hat0,S_new,L_old,grad_1)%*%U_old)
    L_new= U_new%*%t(U_new)
    

    diff_0=norm(L_new-L_old,type="2")/(norm(L_old,type="2")+0.001)+norm(S_new-S_old,type="2")/(norm(S_old,type="2")+0.001)
    L_old=L_new
    
    U_old=U_new
    S_old=S_new
  }
  Theta_0=S_new-L_new
  Theta_0[abs(Theta_0) < 1e-2] <- 0
  return(list(S_0=S_new,L_0=L_new,U_0=U_new,Theta_0=Theta_0))
}
# MCP proximal pperator function
Mcp_prox= function(M,lambda_1,a=3){
  p= nrow(M)
  M_diag=diag(diag(M))
  M_tri=triu(M,1)
  M_tri=apply(M_tri, 1:2,MCP_soft,lambda=lambda_1,a)
  M_prox=M_tri+t(M_tri)+M_diag
  return(M_prox)
}

# Soft threshold function
S_soft = function(z,lambda){
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*sign(z))
}

# MCP threshold operation
MCP_soft = function(z,lambda,a=3){
  return( S_soft(z,lambda)/(1-1/a) * (abs(z) - a*lambda <= 0) + z * (abs(z) - a*lambda > 0) )
}
# Calculate MCP derivative 
mcp_d = function(x,lambda,a=3){
  
  if(lambda!=0){
    rho <- lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

# Density Calculation functions
f_den_vec = function(data, Mu, Theta){ 
  p = length(Mu)
  fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(Mu)) %*% Theta %*% (t(data) - as.numeric(Mu)))) )
  return(fden)
}

f_den = function(data, Mu, Theta){ 
  p = length(Mu)
  fden = as.numeric((det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(Mu)) %*% Theta %*% (t(data) - as.numeric(Mu)))) )
  return(fden)
}

# Matrix symmetrization function
Symmetrize = function(X){
  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}
# Cluster difference truncation function
cut_diff_ama = function(V_kk,K_c,K,cutoff=0.01){
  V_kk_num = which(V_kk < cutoff)
  K_group_final = list()
  if(length(V_kk_num) > 0){
    K_group = list()
    for (j in 1:length(V_kk_num)) {
      K_group[[j]] = K_c[,V_kk_num[j]]
    }
    outnum = setdiff(1:K,Reduce(union,K_group))
    if(length(outnum) > 0){
      for (j in 1:length(outnum)) {
        K_group[[length(V_kk_num)+j]] = outnum[j]
      }
      K_group[[length(V_kk_num)+j+1]] = K
    } else{
      K_group[[length(V_kk_num)+1]] = K
    }
    
    kk = 1
    repeat{
      repeat{
        K_group_old=K_group
        k_del = NULL
        for (kkk in setdiff(1:length(K_group),1) ) {
          if(length(Reduce(intersect,list(K_group[[1]],K_group_old[[kkk]]))) > 0){
            K_group[[1]] = sort(unique(c(K_group[[1]],K_group_old[[kkk]])))
            k_del = c(k_del,kkk)
          }
        }
        if(length(k_del) > 0){
          for (j in sort(k_del,decreasing = T)) {
            K_group[[j]] = NULL
          }
        }
        if(length(K_group_old) == length(K_group)){break}
      }
      K_group_final[[kk]] = K_group[[1]]
      if(kk==1 && length(K_group) == 1){print("Warning: Only one cluster!");break}
      if(length(K_group) == 2){K_group_final[[kk+1]] = K_group[[2]];break}
      if(kk>1 && length(K_group) == 1){break}
      K_group[[1]] = NULL
      kk = kk+1
    }
  }else {
    for (k in 1:K) {
      K_group_final[[k]] = k
    }
  }
  return(K_group_final)
}

# K-means initialization function
Initialize_fuc = function(data, K, n.start = 100){  
  n <- dim(data)[1]
  p <- dim(data)[2]
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, K, nstart = n.start)
  member <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  Sam_cov <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[member == k, , drop = FALSE]) )
    Sam_cov[,,k]  <- cov(data[member == k, , drop = FALSE])
    Theta[,,k] <- solve(Sam_cov[,,k]+diag(1e-5,p))
  }
  
  reslut <- list()
  reslut$prob <-  prob
  reslut$Mu <-  Mu
  reslut$Theta <- Theta
  reslut$Sam_cov<-Sam_cov
  reslut$member <- member
  reslut$K <- K
  return(reslut)
}


#  FGGM update function with fixed estimate_0
Update_FGGM = function(data_K,Theta_0,estimate_K_old,lambda_2, lambda_3, a, rho, 
                       eps, maxiter_EM, maxiter_ADMM, maxiter_AMA, average){
  
  n1 <- nrow(data_K)
  p <- ncol(data_K)
  K <- estimate_K_old$K
  
  if(K==1){
    group_final = 1
    K_final =1
    Mu_final =  estimate_K_old$Mu
    Xi_final = array(0, dim = c(p, p, K_final))
    prob_final = estimate_K_old$prob
    G_final = estimate_K_old$G
    t <- 1
    member_final =rep(1,n1)
    
    
    # update precision matrices
    Theta_final = array(0, dim = c(p, p, 1))
    Delta_final = array(0, dim = c(p, p, 1))
    
    Delta = estimate_K_old$Theta[,,1] -Theta_0
    
    Delta_final[,,1]= MCP_soft(Delta,lambda_2,a)
    Theta_final[,,1]= Theta_0+Delta_final[,,1]
    
    Theta_final[abs(Theta_final) < 1e-3] <- 0
    Delta_final[abs(Delta_final) < 1e-3] <- 0
    
  }else{
    K_c <- combn(K,2)
    
    Theta = estimate_K_old$Theta
    Mu = estimate_K_old$Mu
    prob = estimate_K_old$prob
    K=estimate_K_old$K
    
    # EM algorithm 
    f_mat = matrix(0,n1,K)
    G_mat = matrix(0,n1,K)
    t = 0
    diff_K = 10
    
    while((t<maxiter_EM)&(diff_K>=eps))
    {  
      prob_old = prob
      Mu_old = Mu
      Theta_old = Theta
      G_mat_old = G_mat
      
      # calculate pdfs
      for(k_ind in 1:K){
        f_mat[,k_ind]=as.numeric(det(Theta_old[,,k_ind]))^(1/2) * exp(-1/2*diag(t(t(data_K) - as.numeric(Mu_old[k_ind,])) %*% Theta_old[,,k_ind] %*% (t(data_K) - as.numeric(Mu_old[k_ind,]))))
      }
      # update G and pi
      for(k_ind in 1:K) {
        for(i in 1:n1) {
          G_mat[i,k_ind] = prob_old[k_ind] * f_mat[i,k_ind] / prob_old %*% f_mat[i,]
        }
      }
      G_mat[is.na(G_mat)] = 0
      prob =apply(G_mat,2,mean) 
      nK = apply(G_mat,2,sum)
      
      Theta_kk_diff <- rep(0,dim(K_c)[2])
      for (l in 1:dim(K_c)[2]) {
        Theta_kk_diff[l] <- sum((Theta_old[,,K_c[1,l]] - Theta_old[,,K_c[2,l]])^2)
      }
      
      # update mean vectors                       
      for(j in 1:p){   
        for(k_ind in 1:K) {  
          tmp = t(t(data_K) - as.numeric(Mu[k_ind,])) %*% Theta_old[,j,k_ind] + Mu[k_ind,j] * Theta_old[j,j,k_ind]
          
          hj = t(G_mat[,k_ind])%*%tmp
          tau_k = sqrt(apply((t(Mu[k_ind,] - t(Mu[-k_ind,])))^2,1,sum) + Theta_kk_diff[ceiling(which(K_c == k_ind)/2)])
          v_k = sum(mcp_d(tau_k, lambda_3, a) / (tau_k+0.00001) * Mu[-k_ind,j])
          v_k_hat = sum(mcp_d(tau_k,lambda_3, a) / (tau_k+0.00001))
          Mu[k_ind,j] = (hj + n1*v_k) / (nK[k_ind] * Theta_old[j,j,k_ind] + n1*v_k_hat+0.00001)
        }
      }
      # update precision matrices
      Mu_kk_diff <- rep(0,dim(K_c)[2])
      for (l in 1:dim(K_c)[2]) {Mu_kk_diff[l] <- sum((Mu[K_c[1,l],] - Mu[K_c[2,l],])^2)}
      Sam_cov = array(0, dim = c(p, p, K))
      for (k_ind in 1:K) {
        L_ikx = sqrt(G_mat[,k_ind])*t(t(data_K) - Mu[k_ind,])
        Sam_cov[,,k_ind] = t(L_ikx) %*% L_ikx / nK[k_ind]
      }
      Sam_cov[is.na(Sam_cov)] = 0
      
      Theta_out = Update_ThetaK(Sam_cov,Theta_0,prob,lambda_2,lambda_3, Mu_kk_diff, K_c, a, rho , maxiter_ADMM, maxiter_AMA, eps)
      Theta = Theta_out$Theta
      Xi = Theta_out$Xi
      V_kk = Theta_out$V_kk
      
      t = t + 1
      diff_Mu = norm(Mu_old-Mu,type="2")/(norm(Mu,type="2")+0.001)
      diff_thetaK = norm(Theta_old-Theta,type="2")/(norm(Theta,type="2")+0.001)
      diff_K=diff_Mu+diff_thetaK 
    }
    
    group_final = cut_diff_ama(V_kk,K_c,K,cutoff=0.01)
    K_final = length(group_final)
    Mu_final = matrix(0, K_final, p)
    Theta_final = array(0, dim = c(p, p, K_final))
    Delta_final = array(0, dim = c(p, p, K_final))
    Xi_final = array(0, dim = c(p, p, K_final))
    prob_final = rep(0,K_final)
    G_final = matrix(0,n1,K_final)
    for (l in 1:K_final){
      gg = group_final[[l]]
      prob_final[l] = sum(prob[gg])
      if(length(gg) > 1){
        G_final[,l] = apply(G_mat[,gg],1,sum)
        Mu_final[l,] = apply(Mu[gg,],2,mean)
        
        if(!average){
          Theta_final[,,l] = Theta[,,gg[which.max(nK[gg])]] 
          Xi_final[,,l] = Xi[,,gg[which.max(nK[gg])]]       
        } else {
          Theta_final[,,l] = Theta[,,gg[1]]/length(gg)
          Xi_final[,,l] = Xi[,,gg[1]]/length(gg)
          for (gi in gg[-1]) {
            Theta_final[,,l] = Theta_final[,,l] + Theta[,,gi]/length(gg)
            Xi_final[,,l] = Xi_final[,,l] + Xi[,,gi]/length(gg)
          }
        }
        
      }else{ Mu_final[l,] = Mu[gg,]; Theta_final[,,l] = Theta[,,gg]; Xi_final[,,l] = Xi[,,gg]; G_final[,l]=G_mat[,gg]}
    }
    
    for(k in 1:K_final) {
      Theta_final[,,k] = Symmetrize(Theta_final[,,k])
      Xi_final[,,k] = Symmetrize(Xi_final[,,k])
    }
    Theta_final[abs(Theta_final) < 1e-3] <- 0
    Xi_final[abs(Xi_final) < 1e-3] <- 0
    member_final = apply(G_final,1,which.max)
  }
  
  estimate_K <- list()
  estimate_K$K <- K_final
  estimate_K$Mu <-Mu_final
  estimate_K$Theta <- Theta_final
  estimate_K$Delta <- Xi_final
  estimate_K$prob<- prob_final
  estimate_K$G<- G_final
  estimate_K$member <- member_final
  estimate_K$maxiter_EM <- t
  return(estimate_K)
}
# Update Theta^k via ADMM algorithm
Update_ThetaK= function(Sam_cov,Theta_0, prob, lambda_2, lambda_3, Mu_kk_diff, K_c, a = 3, rho=1, maxiter_ADMM, maxiter_AMA, eps){
  
  p = dim(Sam_cov)[1]
  K = dim(Sam_cov)[3]
  
  prob=prob+1e-5
  
  # Initialize Theta:
  Theta = array(0, dim = c(p, p, K))
  for (k in 1:K) { Theta[,,k] =  diag(1,p)}
  # Initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  # Initialize Phi:
  Phi = array(0, dim = c(p, p, K))
  
  iter=0
  diff_value = 10
  diff_value_Xi = 10
  while( iter<maxiter_ADMM && !(diff_value < eps || diff_value_Xi < eps) )
  {
    Theta_prev = Theta
    Xi_prev = Xi
    Theta=Theta_prev
    # Update Theta:
    for(k_ind in 1:K){
      Sa=Sam_cov[,,k_ind] - rho*Theta_0/(prob[k_ind]) - rho*Xi[,,k_ind]/(prob[k_ind]) + rho*Phi[,,k_ind]/(prob[k_ind])
      Sa[is.na(Sa)] = 0
      edecomp = eigen(Sa)
      D = edecomp$values
      if(is.complex(D)){Sa=Symmetrize(Sa);edecomp = eigen(Sa);D = edecomp$values}
      V = edecomp$vectors
      D2 = prob[k_ind]/(2*rho) * ( -D + sqrt(D^2 + 4*(rho/prob[k_ind])))
      Theta[,,k_ind] = V %*% diag(D2) %*% t(V)
    }
    # Update Xi:
    # Define B matrices:
    B = array(0, dim = c(p, p, K))
    for(k in 1:K){ B[,,k] = Theta[,,k] + Phi[,,k]-Theta_0}
    Xi_out_list = AMA_XI(B,K_c,lambda_2,lambda_3,Mu_kk_diff,a,kappa=rho,maxiter_AMA,eps)
    Xi = Xi_out_list$Xi
    V = Xi_out_list$V
    V_kk = round(apply(V^2,3,sum),4)
    
    # Update the dual variable Phi:
    for(k in 1:K){ Phi[,,k] = Phi[,,k]+Theta[,,k] -Theta_0-Xi[,,k]}
    
    iter = iter+1
    diff_value = sum(abs(Theta - Theta_prev)) / sum(abs(Theta_prev)+0.001)
    diff_value_Xi = sum(abs(Xi - Xi_prev)) / (sum(abs(Xi_prev))+0.001)
  }
  diff = sum(abs(Theta-Xi))
  Theta_out = list(Theta=Theta,Xi=Xi,V_kk=V_kk,diff=diff,iters=iter)
  return(Theta_out)
}
# AMA algorithm for updating Xi
AMA_XI = function(B, K_c, lambda_2, lambda_3, Mu_kk_diff, a, kappa=1,maxiter_AMA,eps){
  
  p = dim(B)[1]
  K = dim(B)[3]
  num_Kc = dim(K_c)[2]
  # Initialize Xi:
  Xi = array(0, dim = c(p, p, K))
  # Initialize V:
  V = array(0, dim = c(p, p, num_Kc))
  # Initialize Delta:
  Delta = array(0, dim = c(p, p, num_Kc))
  Z = array(0, dim = c(p, p, K))
  Omega = array(0, dim = c(p, p, num_Kc))
  e_k12 = matrix(0,K,num_Kc)
  for (i in 1:K) {e_k12[i,which(K_c[1,] == i)] = 1;e_k12[i,which(K_c[2,] == i)] = -1}
  # Iterations
  iter=0
  diff_value = 10
  while( iter<maxiter_AMA && diff_value > eps )
  {
    V_prev = V
    Xi_prev = Xi
    Delta_prev = Delta
    for (i in 1:p) {
      for (j in 1:p) {
        Z[i,j,] = B[i,j,] + apply(t(Delta[i,j,] * t(e_k12)),1,sum)
      }
    }
    
    Xi = S_soft(Z,lambda_2)/(1-1/a) * (abs(Z) <= a*lambda_2) + Z * (abs(Z) > a*lambda_2)
    
    for (l in 1:num_Kc) {
      Omega[,,l] = Xi[,,K_c[1,l]] - Xi[,,K_c[2,l]] - Delta[,,l]/kappa
      Omegal2 = sum(Omega[,,l]^2)
      V[,,l] = Omega[,,l] * MCP_soft(sqrt(Omegal2 + Mu_kk_diff[l]),lambda_3/kappa) / (sqrt(Omegal2+Mu_kk_diff[l])+1e-8) 
      Delta[,,l] = Delta[,,l] + kappa * ( V[,,l] - Xi[,,K_c[1,l]] + Xi[,,K_c[2,l]] )* (sum(V[,,l]^2) > 0 )
    }
    
    iter = iter+1
    diff_value = sum(abs(Xi - Xi_prev)^2) / (sum(abs(Xi_prev)^2)+0.001)
    
  }
  Xi_out = list(Xi=Xi,V=V,Delta=Delta,diff=diff_value,iters=iter)
  return(Xi_out)
}

# Refitting function for FGGM model
FGGM_refit = function(data_K,Theta_0,K,lambda_2, lambda_3, a , rho , 
                      eps, maxiter_EM, maxiter_ADMM, maxiter_AMA,average){
  estimate_ini= Initialize_fuc(data_K,K)
  estimate_K= Update_FGGM (data_K,Theta_0,estimate_ini,lambda_2,lambda_3=0,a,rho , 
                           eps, maxiter_EM, maxiter_ADMM, maxiter_AMA,average)
  return(estimate_K)
}

# Main function for joint estimation 
Update_joint= function(data_0,data_K,K,r,lambda_1,lambda_2 ,lambda_3,eta1,eta2,a = 3, rho = 1, 
                       eps = 1e-2, maxiter_BGD =5 ,maxiter_EM = 20, maxiter_ADMM=10, maxiter_AMA=5, average=F){
  n0 <- nrow(data_0)
  n1 <- nrow(data_K)
  p <- ncol(data_0)
  
  estimate_0=Initialize_fuc_svd(data_0,lambda_1,r)
  S_ini=estimate_0$S_0
  U_ini=estimate_0$U_0
  estimate_K=Initialize_fuc(data_K,K)
  
  t = 0
  diff_loss = 10
  
  while((t<=maxiter_BGD)&(diff_loss>eps))
  {
    grad_1=Grad_1(data_K,estimate_K)
    estimate_0=Update_ATGD(data_0,S_ini,U_ini,eta1,eta2,lambda_1,r,eps,grad_1)
    loss_old=AIC_all(data_0,data_K,estimate_0,estimate_K)$fit_error
    
    estimate_fused=Update_FGGM(data_K,estimate_0$Theta_0,estimate_K,lambda_2, lambda_3, a, rho, 
                               eps, maxiter_EM, maxiter_ADMM, maxiter_AMA, average)
    
    estimate_K=FGGM_refit(data_K,estimate_0$Theta_0,estimate_fused$K,lambda_2, lambda_3 ,a, rho, 
                          eps, maxiter_EM, maxiter_ADMM, maxiter_AMA, average)
    a<<-estimate_K
    loss_new=AIC_all(data_0,data_K,estimate_0,estimate_K)$fit_error
    
    diff_loss=abs(loss_new-loss_old)/(abs(loss_old)+0.001)
    diff_loss[is.na(diff_loss)]=0
    
    t = t + 1
  }
  Estimate=list()
  Estimate$estimate_0=estimate_0
  Estimate$estimate_K=estimate_K
  return(Estimate)
}

# AIC calculation functions 
AIC_0=function(data_0,estimate_0){
  
  p=ncol(data_0)
  r=ncol(estimate_0$U_0)
  n0=nrow(data_0)
  
  Theta_hat0=estimate_0$Theta_0
  S_hat0=estimate_0$S_0
  
  
  fit0 = f_den_vec(data_0,rep(0,p), Theta_hat0) 
  fit_error0 = sum(log(fit0 + min(fit0[fit0>0])))
  fit_error = -2*fit_error0/n0
  
  df = 2/n0*(length(which(S_hat0!= 0))+p*r)
  
  return(fit_error + df)
}

AIC_K=function(data_K,estimate_K){
  
  K =estimate_K$K
  n1=nrow(data_K)
  
  Mu_hat=estimate_K$Mu
  Theta_hat=estimate_K$Theta
  prob_hat=estimate_K$prob
  
  fit_error_mat = matrix(0, n1, K)
  for(k in 1:K) {
    fit_error_mat[,k] = prob_hat[k] * f_den_vec( data_K, as.numeric(Mu_hat[k,]), Theta_hat[,,k] )                  
  }
  fit1 = apply(fit_error_mat, 1, sum)
  fit_error = sum(log( fit1 + min(fit1[fit1>0])))
  fit_error = -2*fit_error/n1
  df = 2/n1*length(which(Mu_hat!= 0)) +2/n1*(length(which(Theta_hat!= 0)))
  return(fit_error + df)
}

AIC_all=function(data_0,data_K,estimate_0,estimate_K){
  
  K =estimate_K$K
  p=ncol(data_0)
  r=ncol(estimate_0$U_0)
  n0=nrow(data_0)
  n1=nrow(data_K)
  
  Mu_hat=estimate_K$Mu
  Theta_hat=estimate_K$Theta
  prob_hat=estimate_K$prob
  Delta_hat=estimate_K$Delta
  
  Sigma_hat0=cov(data_0)
  Theta_hat0=estimate_0$Theta_0
  S_hat0=estimate_0$S_0
  
  fit_error_mat = matrix(0, n1, K)
  for(k in 1:K) {
    fit_error_mat[,k] = prob_hat[k] * f_den_vec( data_K, as.numeric(Mu_hat[k,]), Theta_hat[,,k] )                  
  }
  fit1= apply(fit_error_mat, 1, sum)
  fit0 = f_den_vec(data_0,rep(0,p), Theta_hat0) 
  fit_error1 = sum(log(fit1 + min(fit1[fit1>0])))
  fit_error0 = sum(log(fit0 + min(fit0[fit0>0])))
  fit_error = -2*fit_error1/n1+ -2*fit_error0/n0
  
  df = 2/n1*(length(which(Mu_hat != 0)))+ 2/n1*length(which(Delta_hat != 0))+2/n0*(length(which(S_hat0!= 0))+p*r)
  
  AIC = fit_error + df
  P = list()
  P$fit_error = fit_error
  P$df = df
  P$AIC = AIC
  return(P)
}
# Lambda selection via line search
Select_lambda<-function(data_0,data_K,K,r,lambda1_seq,lambda2_seq,lambda3_seq,eta1,eta2,a = 3, rho = 1, 
                        eps = 1e-2, maxiter_BGD =5 ,maxiter_EM = 20, maxiter_ADMM=10, maxiter_AMA=5, average=F)
{
  
  n_lambda1<- length(lambda1_seq)
  n_lambda2<- length(lambda2_seq)
  n_lambda3<- length(lambda3_seq)
  
  lambda_1 = median(lambda1_seq)
  lambda_2 = median(lambda2_seq)
  
  aAIC=c()
  for (l in 1:n_lambda3){
    lambda_3=lambda3_seq[l]
    fit=Update_joint(data_0,data_K,K,r,lambda_1,lambda_2 ,lambda_3,eta1,eta2,a , rho, 
                     eps, maxiter_BGD ,maxiter_EM , maxiter_ADMM, maxiter_AMA, average)
    aAIC[l] = AIC_all(data_0,data_K,fit$estimate_0,fit$estimate_K)$AIC
  }
  lambda_3 = lambda3_seq[which.min(aAIC)]
  
  aAIC=c()  
  
  for (l in 1:n_lambda2){
    lambda_2=lambda2_seq[l]
    fit=Update_joint(data_0,data_K,K,r,lambda_1,lambda_2 ,lambda_3,eta1,eta2,a , rho, 
                     eps, maxiter_BGD ,maxiter_EM , maxiter_ADMM, maxiter_AMA, average)
    aAIC[l] = AIC_all(data_0,data_K,fit$estimate_0,fit$estimate_K)$AIC
    
  }
  lambda_2 = lambda2_seq[which.min(aAIC)]
  
  
  aAIC=c()
  result=list()
  for (l in 1:n_lambda1){
    lambda_1=lambda1_seq[l]
    fit=Update_joint(data_0,data_K,K,r,lambda_1,lambda_2 ,lambda_3,eta1,eta2,a , rho, 
                     eps, maxiter_BGD ,maxiter_EM , maxiter_ADMM, maxiter_AMA, average)
    aAIC[l] = AIC_all(data_0,data_K,fit$estimate_0,fit$estimate_K)$AIC
    result[[l]]=fit
  }
  lambda_1 = lambda1_seq[which.min(aAIC)]
  result_final=result[[which.min(aAIC)]] 
  print(c(cat(paste("lam1 lam2 lam3 ="),lambda_1,lambda_2,lambda_3,":"),paste("K =",result_final$estimate_K$K)))
  return(result_final)
}
