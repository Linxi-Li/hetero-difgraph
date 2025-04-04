# Hetero-difgraph
This repository contains our R implementation  of "Subgroup Analysis of Differential Networks with Latent Variables"

The following .R files are submitted.

main.R: This file includes all main functions of proposed methods to support numerical simulation studies and real data analysis. Detailed descriptions and manuals of these functions can be found in the notes to this R file.

data.gene.R: This file includes functions for generating simulated data and evaluating performances of competitors, which are used to support numerical simulation studies. Detailed descriptions and manuals of these functions can be found in the notes to this R file.

demo.R: This file contains example code and usage demonstrations for implementing the proposed methods.

# Useful parameters


data_0: n*p data set of baseline group

data_K: n*p data set of target heterogeneous group

K: the initial number of subgroups  
 
r: Rank constraint  

lambda_1, lambda_2, lambda_3: Regularization parameters

eta1, eta2: Learning rates for S_0 and U_0
 
a, rho: Tuning parameters for MCP and ADMM, dedault as a=3, rho=1

 eps: Convergence threshold
 
maxiter_*: Maximum iterations for each step

average: Logical value, indicating the calculation form after fusion 


# An example in demo.R
```
# Generate simulated dataset

data = Generate_data(n0=1500, n1=1500, prob=c(1/3,1/3,1/3), p=100, r=2, sparsity_Theta=0.02, K_true=3, mue=1.5, seed=123)
data_0 = data$data_0  # baseline group data
data_K = data$data_K  # heterogeneous group data
Para_true = data$Para_true  # True parameters for evaluation

# Initialize model parameters
K = 6  # Number of subgroups (larger than true K)
r = 2  # Latent factor dimension

# Fit joint model with fixed regularization parameters

fit = Update_joint(data_0, data_K, K, r, lambda_1=0.01, lambda_2=0.01, lambda_3=2,  eta1=0.01, eta2=0.01, a=3, rho=1, eps=1e-2,  maxiter_BGD=5, maxiter_EM=10, maxiter_ADMM=10, maxiter_AMA=5, average=F)

# Evaluate model performance against true parameters
mt=Metric(data_0, data_K, fit$estimate_0, fit$estimate_K, Para_true)
print(mt)

# Parameter tuning sequence setup
lambda1_seq = seq(0.01, 0.1, len=5)  # Lambda1 candidates
lambda2_seq = seq(0.01, 0.1, len=5)  # Lambda2 candidates
lambda3_seq = seq(1, 4, len=5)       # Lambda3 candidates

# Automatically select optimal regularization parameters using AIC
fit_opt = Select_lambda(data_0, data_K, K, r, lambda1_seq, lambda2_seq, lambda3_seq, eta1=0.01, eta2=0.001, a=3, rho=1, eps=1e-2, maxiter_BGD=5, maxiter_EM=10, maxiter_ADMM=10, maxiter_AMA=5, average=F)

# Evaluate optimized model
mt=Metric(data_0, data_K, fit_opt$estimate_0, fit_opt$estimate_K, Para_true)
print(mt)

# Visualize estimated network structure for subgroup k
Plot_network(Delta=fit_opt$estimate_K$Delta, subgroup=1)
```
