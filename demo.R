

data=Generate_data(n0=1500,n1=1500,prob=c(1/3,1/3,1/3),p=100,r=2,sparsity_Theta=0.02,K_true=3,mue=1.5,seed=123)
data_0=data$data_0
data_K=data$data_K
Para_true=data$Para_true
K=6
r=2

fit=Update_joint(data_0,data_K,K,r,lambda_1=0.01,lambda_2=0.01,lambda_3=2,eta1=0.01,eta2=0.01,a = 3, rho = 1, 
                        eps = 1e-2, maxiter_BGD =5 ,maxiter_EM = 10, maxiter_ADMM=10, maxiter_AMA=5, average=F)
Metric(data_0,data_K,fit$estimate_0,fit$estimate_K,Para_true)


lambda1_seq =seq(0.01,0.1,len=5)
lambda2_seq =seq(0.01,0.1,len=5)
lambda3_seq =seq(1,4,len=5)
fit_opt=Select_lambda(data_0,data_K,K,r,lambda1_seq,lambda2_seq,lambda3_seq,eta1=0.01,eta2=0.001,a = 3, rho = 1, 
                        eps = 1e-2, maxiter_BGD =5 ,maxiter_EM = 10, maxiter_ADMM=10, maxiter_AMA=5, average=F)
Metric(data_0,data_K,fit_opt$estimate_0,fit_opt$estimate_K,Para_true)

Plot_network(Delta=fit_opt$estimate_K$Delta,subgroup = 2)

