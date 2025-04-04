############################# Functions for numerical experiment ############################
#Function to generate random symmetric matrix
Random_symmat <- function(p, sparsity_Theta) {
  mat <- matrix(runif(p**2, min=-1,max=1) * rbinom(n=p**2, size=1,prob = sparsity_Theta/2),p)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  return(mat)
}
#Data generation function
Generate_data = function(n0,n1,prob,p,r,K_true,sparsity_Theta,mue,seed){
  set.seed(seed)
  Para_true=Generate_theta(p, r,K_true,sparsity_Theta,mue)
  
  data_0=mvrnorm(n0,matrix(0,p),solve(Para_true$Theta_0))
  member_K=sample(1:K_true,n1,replace=T,prob=prob)
  data_K=matrix(0,n1,p)
  for (k in 1:K_true){
    data_K[member_K==k,] = mvrnorm(length(which(member_K==k)),Para_true$Mu_K[k,],solve(Para_true$Theta_K[,,k]))
  }
  Para_true$member_K=member_K
  return(list(data_0=data_0,data_K=data_K,Para_true=Para_true))
}
# Parameter generation function
Generate_theta <- function(p,r,K_true,sparsity_Theta=0.02,mue=1.5) {
  
  p_all=p+r
  
  Theta0_tilde=Random_symmat(p_all,sparsity_Theta)
  Theta0_tilde<- Theta0_tilde+abs(min(eigen(Theta0_tilde)$val)) * diag(p_all) + 0.1*diag(p_all)
  
  
  S_0= Theta0_tilde[1:p,1:p]
  L_0 = Theta0_tilde[1:p,(p+1):p_all]%*%solve((Theta0_tilde[(p+1):p_all,(p+1):p_all]))%*%Theta0_tilde[(p+1):p_all,1:p]
  Theta_0 = S_0-L_0
  
  Theta_K=array(0, dim = c(p,p,K_true))
  Delta_K=array(0, dim = c(p,p,K_true))
  
  Mu_K=matrix(0,K_true,p)
  Mu_K[1,]= c(rep(mue,4),rep(-mue,4),rep(0,p-2*4))
  Mu_K[2,]= c(rep(mue,2*4),rep(0,p-2*4))
  Mu_K[3,]= c(rep(-mue,2*4),rep(0,p-2*4))
  
  for (k in 1:K_true){
    
    Dif_k=Random_symmat(p,sparsity_Theta)
    Theta_k=Dif_k+Theta_0
    
    Theta_K[,,k]=Theta_k+abs(min(eigen(Theta_k)$val)) * diag(p) + 0.1*diag(p)
    Delta_K[,,k]=Theta_K[,,k]-Theta_0
    
  }
  return(list(S_0=S_0,L_0=L_0,Theta_0=Theta_0,Theta_K=Theta_K,Delta_K=Delta_K,Mu_K=Mu_K,K_true=K_true))
}
# Metric calculation function
Metric = function(data_0,data_K,estimate_0,estimate_K,Para_true){
  
  Theta0_true=Para_true$Theta_0
  S0_true=Para_true$S_0
  L0_true=Para_true$L_0
  
  S0_hat=estimate_0$S_0
  L0_hat=estimate_0$L_0
  Theta0_hat=estimate_0$Theta_0
  
  K_true=Para_true$K_true
  Mu_true=Para_true$Mu_K
  DeltaK_true=Para_true$Delta_K
  ThetaK_true=Para_true$Theta_K
  memb_true=Para_true$member_K
  memb_hat=estimate_K$member
  
  
  Mu_hat=estimate_K$Mu
  K_hat= estimate_K$K
  DeltaK_hat=estimate_K$Delta
  ThetaK_hat=estimate_K$Theta
  prob_hat=estimate_K$prob
  
  n1 = nrow(data_K)
  p=ncol(data_K)
  
  #  MSE_S0 & MSE_L0
  MSE_S0 = norm(S0_hat- S0_true,type="F")
  MSE_L0 = norm(L0_hat- L0_true,type="F")
  MSE_Theta0 = norm(Theta0_hat- Theta0_true,type="F")
  
  
  
  if(K_hat == K_true){
    per=1
    num = rep(0,K_true)
    numk = NULL
    for (k in 1:K_true) {
      errork = apply((t(Mu_hat) - Mu_true[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu_hat = Mu_hat[num,]
    DeltaK_hat = DeltaK_hat[,,num]
    #  MSE_Mu & MSE_DeltaK & MSE_ThetaK
    MSE_Mu = mean(sqrt(diag((Mu_hat - Mu_true)%*%t(Mu_hat - Mu_true))))
    MSE_DeltaK =mean(apply(DeltaK_hat-DeltaK_true,3,norm,type="F"))
    
    #  TPR & FPR of Delta
    TPR_DeltaK = sum(((apply(DeltaK_hat, 3, function(x) x[upper.tri(x)])!=0) + (apply(DeltaK_true, 3, function(x) x[upper.tri(x)])!=0))== 2)/ (sum(apply(DeltaK_true, 3, function(x) x[upper.tri(x)])!=0))
    FPR_DeltaK =sum(((apply(DeltaK_hat, 3, function(x) x[upper.tri(x)])!=0) + (apply(DeltaK_true, 3, function(x) x[upper.tri(x)])==0))== 2)/ (sum(apply(DeltaK_true, 3, function(x) x[upper.tri(x)])==0))
    
    
    #  CE
    aa =memb_true
    cap_matrix0 = matrix(0,n1,n1)
    for(i in 1:(n1-1)){
      for (j in (i+1):(n1)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = memb_hat
    cap_matrix1 = matrix(0,n1,n1)
    for(i in 1:(n1-1)){
      for (j in (i+1):(n1)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n1*(n1-1)/2)
  } else{
    per=0
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      Mu_hatk = Mu_hat[k,]
      errork.Mu = apply((t(Mu_true) - Mu_hatk)^2,2,sum)
      DeltaK_hatk = DeltaK_hat[,,k]
      errork.Delta = apply((DeltaK_true - rep(DeltaK_hatk,K_true))^2,3,sum)
      errork = errork.Mu + errork.Delta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    
    Mu_true.re = Mu_true[num,]
    DeltaK_true.re = DeltaK_true[,,num]
    ThetaK_true.re = ThetaK_true[,,num]
    if(K_hat == 1){
      Mu_true.re = as.matrix(t(Mu_true.re))
      DeltaK_true.re = array(DeltaK_true.re,dim=c(p,p,1))
      ThetaK_true.re = array(ThetaK_true.re,dim=c(p,p,1))
    }
    
    #  MSE_Mu & MSE_DeltaK
    MSE_Mu = mean(sqrt(diag((Mu_hat - Mu_true.re)%*%t(Mu_hat - Mu_true.re))))
    MSE_DeltaK =mean(apply(DeltaK_hat-DeltaK_true.re,3,norm,type="F"))
    MSE_ThetaK =mean(apply(ThetaK_hat-ThetaK_true.re,3,norm,type="F"))
    
    #  TPR & FPR of Delta
    
    TPR_DeltaK =sum(((apply(DeltaK_hat, 3, function(x) x[upper.tri(x)])!=0) + (apply(DeltaK_true.re, 3, function(x) x[upper.tri(x)])!=0))== 2)/ (sum(apply(DeltaK_true.re, 3, function(x) x[upper.tri(x)])!=0))
    FPR_DeltaK = sum(((apply(DeltaK_hat, 3, function(x) x[upper.tri(x)])!=0) + (apply(DeltaK_true.re, 3, function(x) x[upper.tri(x)])==0))== 2)/ (sum(apply(DeltaK_true.re, 3, function(x) x[upper.tri(x)])==0))
    
    
    #  CE
    aa =memb_true
    cap_matrix0 = matrix(0,n1,n1)
    for(i in 1:(n1-1)){
      for (j in (i+1):(n1)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = memb_hat
    cap_matrix1 = matrix(0,n1,n1)
    for(i in 1:(n1-1)){
      for (j in (i+1):(n1)) {
        
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n1*(n1-1)/2)
  }
  index = as.data.frame(t(c(MSE_S0,MSE_L0,MSE_Theta0,per,K_hat,CE,MSE_Mu,MSE_DeltaK,TPR_DeltaK,FPR_DeltaK)))
  names(index) = c("MSE_S0","MSE_L0","MSE_Theta0","per","K_hat","CE","MSE_Mu","MSE_DeltaK","TPR_DeltaK","FPR_DeltaK")
  return(index)
}

# Function to plot the network of subgroup k based on the parameters delta 
Plot_network <- function(Delta,subgroup) {
  Delta.graph <- (abs(Delta[,,subgroup]) > 1e-5) * 1
  rownames(Delta.graph) <- colnames(Delta.graph) <- colnames(data_K)
  diag(Delta.graph) <- 0
  net.dtrace <- graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  l.dtrace <- layout_with_fr(net.dtrace)
  degrees <- degree(net.dtrace)
  node_colors <-  "#9ECAE1"
  plot(net.dtrace, layout = l.dtrace, vertex.size =5, vertex.color = node_colors, 
       vertex.label.cex =1, edge.width = 2, edge.color = "#E64B357F", vertex.frame.color = "transparent")
}