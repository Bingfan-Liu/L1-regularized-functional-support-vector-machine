############################################################################
#                     One standard error rule illustration

one_se_data = readRDS("one_se_data_for_l1fsvm_n=200_p=20.rds")
one_se_data

dev.new(width=4, height=3, unit="in")
plot(1:length(one_se_data$accuracy), one_se_data$accuracy,
     xaxt = "n", xlab='cost parameter lambda',
     ylab="average accuracy", ylim=c(0.6, 0.8), type = "o",
     pch=19)
axis(1, at=1:length(one_se_data$accuracy), labels=one_se_data$test_cases[,1])
title(main="Average accuracy vs cost parameter lambda")
for (i in 1:length(one_se_data$accuracy)) {
  segments(i,
           one_se_data$accuracy[i]-one_se_data$accuracy_standard_error[i],
           i,
           one_se_data$accuracy[i]+one_se_data$accuracy_standard_error[i])
  segments(i-0.1,
           one_se_data$accuracy[i]+one_se_data$accuracy_standard_error[i],
           i+0.1,
           one_se_data$accuracy[i]+one_se_data$accuracy_standard_error[i])
  segments(i-0.1,
           one_se_data$accuracy[i]-one_se_data$accuracy_standard_error[i],
           i+0.1,
           one_se_data$accuracy[i]-one_se_data$accuracy_standard_error[i])
}





##########################################################################
#                        Simulation Prediction Accuracy
#l1fsvm_files = list.files(pattern="*l1fsvm.rds")
#flg_files = list.files(pattern="*grplFlogit.rds")
files = list.files(pattern="*.rds")
# prediciton accuracy comparison
x11(width=8, height=14)
par(mar=c(3, 2, 3, 1), mfrow=c(2,3))
for(i in 1:length(files)){
  df = readRDS(files[i])
  #l1fsvm_data = readRDS(l1fsvm_files[i])
  #flg_data = readRDS(flg_files[i])
  
  data = data.frame(L1-fSVM = df[5,],
                    grplFlogit = df[10,])
  boxplot(data, main=paste(strsplit(files[i]," ")[[1]][2],
                           strsplit(strsplit(files[i]," ")[[1]][3], ".rds")[1]),
          ylim=c(0.7,1))
  
}




# FP FN rates comparison
x11(width=8, height=14)
#par(mar=c(2,2,2,1), mfrow=c(3,4))
files = list.files(pattern="*.rds")

cal_one_FPFN <- function(f){
  # return a vector of 2 FP and 2 FN rates under one simulation scenario
  #f=files[[1]]
  df = readRDS(f)#[5,]
  
  data = data.frame("l1fsvm FP rate" = df[1,],
                    "lg FP rate" = df[6,],
                    "l1fsvm FN rate" = df[2,],
                    "lg FN rate" = df[7,])#l1fSVM = l1fsvm_data
  result <- c()
  for (column in 1:dim(data)[2]){
    result <- c(result, paste0(round(median(data[,column]),4), " (", round(sd(data[,column]),4),")"))  
  }
  
  return (result)
}



# vertical version
result_mat = matrix(nrow = 9,ncol=4)
result_df = data.frame(result_mat)
for(i in 1:length(files)){
  f = files[[i]]
  
  
  row_name = paste(strsplit(f," ")[[1]][2],
                      strsplit(strsplit(f," ")[[1]][3], ".rds")[1])
  
  rownames(result_df)[i] <- row_name
  result_df[i,]<- cal_one_FPFN(f)
}
result_df
metrics = c("L1fSVM FP rate", "grplFlogit FP rate",
            "L1fSVM FN rate", "grplFlogit FN rate")
result_df = rbind(metrics, result_df)

write.csv(result_df,"FPFN_vertical.csv", row.names = TRUE)



















################################################################################
#                             Real Data Study
real_data_result = c()
files = list.files(pattern="*.rds")
for(f in files){
  real_data_result = cbind(real_data_result , readRDS(f))
}

real_data_result
par(mar=c(2, 2, 2, 0), mfrow=c(2,1))
hist(real_data_result[1,], main ="Prediction accuracy L1-fSVM", xlim=c(0.2,1), ylim=c(0,50))
abline(v=median(real_data_result[1,]), col="red", lwd=3)
#legend("topleft", legend=c(paste0("median = ",round(median(real_data_result[1,]), 2))),
#       lty=1, lwd=3, col="red")

hist(real_data_result[4,], main ="Prediction accuracy grplFlogit", xlim=c(0.2,1), ylim=c(0,50))
abline(v=median(real_data_result[4,]), col="red", lwd=3)
#legend("topleft", legend=c(paste0("median = ",round(median(real_data_result[4,]), 2))),
#       lty=1, lwd=3, col="red")


l1fsvm <- function(X, y, tps, cost=10, nbasis=15, lr=0.01,
                   svm_type=5, tolerance=0.002, maxiter=10000,
                   track=FALSE, ...){
  
  ######## initialize C, Alpha, Z and scalar svm data input ############
  
  n <- length(y)
  p <- length(tps)
  Basis <- vector(mode = "list", length = p)
  
  C <- matrix(1, nrow=nbasis, ncol=p) 
  C_new <- matrix(1, nrow=nbasis, ncol=p)
  
  
  Z <- vector(mode = "list", length = p)
  
  Alpha0 <- ifelse(sum(y==1)>=sum(y==-1), 1, -1)
  Alpha <- rep(0, p)
  
  converged <- FALSE
  loss_record <- rep(NA, maxiter)
  
  
  
  # create bspline basis on [0,1]
  compute_Zj <- function(j, tps,nbasis, X, Basis){
    
    max_tps <- max(tps[[j]])
    min_tps <- min(tps[[j]])
    length_tps <- length(tps[[j]])
    bBasis_obj <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis)
    
    if(min_tps != 0){
      # transfer tps to [0,1]
      tps_trans <- c(0, (tps[[j]]/(tps[[j]][length_tps-1]))[-length_tps]) 
      delta_T <- tps_trans[2] - tps_trans[1]
      Basis[[j]] <- eval.basis(tps_trans, bBasis_obj) # T by nbasis(q) matrix 
    }else{
      tps_trans <- tps[[j]]/max_tps # transfer tps to [0,1]
      delta_T <- tps_trans[2] - tps_trans[1]
      Basis[[j]] <- eval.basis(tps_trans, bBasis_obj)
    }
    
    X[[j]]%*%Basis[[j]]*delta_T 
  }
  
  
  Z <- lapply(1:p, FUN=function(j){compute_Zj(j, tps,nbasis, X, Basis)})
  
  
  ####################     Start the main loop   ###########################
  
  iter <- 1
  
  while(TRUE){
    
    ################# step 1 ######################
    
    # Z[[j]] here is the (Z_ij)^T in the paper
    # Z[[j]] n by nbasis(q), C[,j] nbasis(q) by 1
    #   ==> Z[[j]]%*%C[,j] n by 1 matrix
    #   ==> ZC an n by p matrix 
    
    ZC <- sapply(1:p, FUN = function(j) Z[[j]]%*%C[,j]) 
    
    W <- LiblineaR(data=ZC, target=y, type=svm_type, cost=cost)$W 
    
    Alpha0_new <- W[p+1] # the last element is the intercept
    Alpha_new <- W[-(p+1)]
    Ind <- Alpha_new != 0 # indicator of wether the coefficient is included in the trained_model (non-zero)
    
    num_coef <- sum(Ind) # number of coefficients that are included by the trained_model except for intercept
    
    
    
    ################### step 2 #######################
    
    ####################### Gradient Descent Start #################################
    
    iter_gd <- 1
    
    if(num_coef != 0){ 
      
      dL_Cj_mat <- matrix(NA,nrow=(nbasis-1),ncol=sum(Ind)) # (q-1) rows by sum(num_j!=0) columns
      
      while(TRUE){
        # C[,Ind] # nbasis (q) by sum(num_j!=0)
        # Z[Ind] # list of length sum(num_j!=0). Each ele is an n by nbasis(q) matrix 
        # Alpha_new[Ind] 
        
        # ZC_ind a matrix of n by sum(num_j!=0) (current num_coef <= p),
        #     each column is Z_jC_j
        # Z[Ind][[j]] n by nbasis (q) 
        # C[,Ind][,j] nbasis(q) by 1
        
        
        ZC_ind <- sapply(1:num_coef,
                         FUN=function(j) Z[Ind][[j]]%*%as.matrix(as.matrix(C[,Ind])[,j])) 
        
        dL_Cj <- function(j, y, Alpha_new, Ind, Z, Alpha0_new, ZC_ind, C, C_new, lr){
          b <- 1- y * (Alpha0_new + ZC_ind %*% matrix(Alpha_new[Ind], ncol=1))
          #print(dim(b))
          #print(dim(- 2 * cost * y * Alpha_new[Ind][j] * Z[Ind][[j]][,-1]))
          dL_Cij <- (- 2 * cost * y * Alpha_new[Ind][j] * Z[Ind][[j]][,-1] *
                       c( b > 0 ) * c(b) )
          return(colSums(dL_Cij))# 1 by (nbasis-1)
        }
        
        dL_Cj_mat <- sapply(1:sum(Ind),FUN=function(j){
          dL_Cj(j, y, Alpha_new, Ind, Z, Alpha0_new, ZC_ind, C, C_new, lr)})
        
        
        if(sum(Ind) != 1){
          C_new[,Ind][-1,] <- C[,Ind][-1,] - lr * dL_Cj_mat
        }else{
          C_new[,Ind][-1] <- C[,Ind][-1] - lr * dL_Cj_mat
        }
        
        
        
        
        if(max(abs((C_new-C))) <= tolerance){
          if(track){
            cat("Gradient Descent converged.\n")
          }
          break
        }
        if(iter_gd >= maxiter){
          if(track){
            cat("Reached maxiter. Gradient Descent didn't converge.\n")
          }
          break
        }
        
        # go for the next iteration
        C <- C_new
        iter_gd <- iter_gd + 1
        
      } 
      
    }else{
      # the alphas shrink to zero
      # or the alphas were zero in the first place
      
      result <- list("Alpha0" = ifelse(sum(y==1) >= sum(y==-1), 1, -1),
                     "Alpha" = Alpha_new,
                     "C" = C_new,
                     "Iteration" = iter,
                     "Convergence" = TRUE,
                     "LossRecord"=sum(sapply(1-y*Alpha0_new, FUN = max, 0)**2)*cost)
      
      return(result) 
      
    }
    
    #################  Gradient Descent Finished   ######################
    
    
    # compute the Total Loss
    ZC <- sapply(1:p, FUN = function(j) Z[[j]]%*%C_new[,j,drop=FALSE]) #n by p matrix
    f <- (rep(Alpha0_new, dim(ZC)[1]) + ZC%*%Alpha_new)
    l <- 1 - (y*f)
    
    total_loss <- sum(( sapply(l,FUN = max,0) )**2)*cost + sum(abs(Alpha_new)) 
    loss_record[iter] <- total_loss 
    
    
    if(track){
      cat("Total Loss: ", total_loss, "\n")
    }
    
    
    
    
    if(max(abs( c(Alpha0_new-Alpha0, Alpha_new-Alpha, C_new-C) )) <= tolerance){  
      
      if(track){
        cat("Algorithm converged.\n")
      }
      converged <- TRUE
      
      break
    }
    
    if(iter >= maxiter){
      cat("Reached maxiter. Algorithm didn't converge.\n")
      break
    }
    
    # go for the next iteration in the main loop
    
    C <- C_new
    Alpha0 <- Alpha0_new
    Alpha <- Alpha_new
    iter <- iter + 1
    
  }
  
  ######################### End of the Loop ###################################
  
  result <- list("Alpha0" = Alpha0_new,
                 "Alpha" = Alpha_new,
                 "C" = C_new,
                 "Iteration" = iter,
                 "Convergence" = converged,
                 "LossRecord"= loss_record[1:iter])
  
  return(result) 
  
}


grplFlogit <- function(Y, X, Tps, lambda, phi, dfs = 20,
                       adapt1 = NULL, adapt2 = NULL, ...){
  
  nsub = length(Y) 
  nfunc = length(Tps)
  
  if (length(dfs) == 1)
    dfs = rep(dfs, nfunc) 
  if (length(dfs) != nfunc)
    stop("length of dfs does not match number of predictors")
  
  B <- Psi <- Omega <- K <- iR <- eK <- list() 
  delt <- rep(NA, nfunc)
  
  for (jj in 1:nfunc){
    
    spj = diff(range(Tps[[jj]]))
    bknj = c(min(Tps[[jj]]) - spj, max(Tps[[jj]]) + spj) 
    B[[jj]] = bs(Tps[[jj]], df=dfs[jj], Boundary.knots=bknj) 
    delt[jj] = Tps[[jj]][2] - Tps[[jj]][1] 
    Psi[[jj]] = delt[jj] * t(B[[jj]]) %*% B[[jj]] 
    
    if (length(adapt1) == nfunc)
      Psi[[jj]] = adapt1[jj]*Psi[[jj]]
    
    dBj <- matrix(NA,nrow(B[[jj]]),ncol(B[[jj]]))
    for (k in 1:ncol(B[[jj]])) 
    {
      iS <- interpSpline(Tps[[jj]],B[[jj]][,k])
      dBj[,k] <- predict(iS, Tps[[jj]], deriv = 2)$y
    }
    Omega[[jj]] = delt[jj] * t(dBj) %*% dBj 
    
    if (length(adapt2) == nfunc)
      Omega[[jj]] = adapt2[jj]*Omega[[jj]]
    
    K[[jj]] = Psi[[jj]] + phi * Omega[[jj]] 
    eK[[jj]] <- eigen(K[[jj]])
    iR[[jj]] = backsolve(chol(K[[jj]]), x = diag(ncol(K[[jj]]))) 
  }
  
  Z = 1
  for (jj in 1:nfunc)
  {
    tmp = delt[jj]*(X[[jj]]%*%B[[jj]])
    Z = cbind(Z, tmp%*%iR[[jj]]) 
  }
  
  
  index = c(NA,rep(1:nfunc,dfs))   
  grpl = grplasso(x = Z, y = Y, index = index, model = LogReg(), lambda = lambda, standardize = F)
  
  intercept = grpl$coef[1,]
  Coef <- list()
  index[1] = 0
  for (jj in 1:nfunc)
  {
    Coef[[jj]] <- B[[jj]]%*%iR[[jj]]%*%grpl$coef[index == jj,]
  }
  
  out = list("intercept" = intercept, "Coef" = Coef)
  return(out)
}



predict_lg <- function(X_valid, Tps, trained_model){
  n <- dim(X_valid[[1]])[1]
  
  eta <- trained_model$intercep + rowSums(sapply(1:length(Tps), FUN=function(j){X_valid[[j]]%*%trained_model$Coef[[j]]*(Tps[[j]][2] - Tps[[j]][1])}))
  
  predicted_prob <- exp(eta)/(1+exp(eta))
  
  result <- rep(NA,n)
  
  result[predicted_prob>0.5] <- 1
  result[predicted_prob<=0.5] <- 0
  
  return(result)
}


############################################################
#################### For L1-fSVM ########################

set.seed(123456)
require(fda)
require(LiblineaR)
require(parallel)

load("UCI_data.RData")
y[y==0] <- -1
pos_ind <- sample(which(y==1), min(sum(y==1), sum(y==-1)))
neg_ind <- sample(which(y==-1), min(sum(y==1), sum(y==-1)))
y <- y[c(pos_ind, neg_ind)]

X <- lapply(1:length(X), function(i){X[[i]][c(pos_ind, neg_ind),]})



X_smoothed <- X
t = 0:255
for(j in 1:length(X)){
  for(i in 1:dim(X[[1]])[1]){
    fit = ksmooth(t, X[[j]][i,], kernel = c("box"), bandwidth = 10)
    X_smoothed[[j]][i,] <- fit$y
  }
}

X_smoothed_central <- X_smoothed
for(j in 1:length(X)){
  mean_curve = colMeans(X_smoothed[[j]])
  for(i in 1:dim(X[[1]])[1]){
    X_smoothed_central[[j]][i,] <- X_smoothed[[j]][i,] - mean_curve
  }
}
result = l1fsvm(X_smoothed_central, y, tps, cost=0.05, nbasis=5, lr=0.01,
                svm_type=5, tolerance=0.0001, maxiter=10000,
                track=FALSE)
result
load("UCI_data.RData")
names(X)[which(result$Alpha != 0)]







############################################################
#################### For grplFlogit ########################

require(grplasso)
set.seed(123456)
load("UCI_data.RData")
pos_ind <- sample(which(y==1), min(sum(y==1), sum(y==0)))
neg_ind <- sample(which(y==0), min(sum(y==1), sum(y==0)))
y <- y[c(pos_ind, neg_ind)]

X <- lapply(1:length(X), function(i){X[[i]][c(pos_ind, neg_ind),]})



X_smoothed <- X
t = 0:255
for(j in 1:length(X)){
  for(i in 1:dim(X[[1]])[1]){
    fit = ksmooth(t, X[[j]][i,], kernel = c("box"), bandwidth = 10)
    X_smoothed[[j]][i,] <- fit$y
  }
}

X_smoothed_central <- X_smoothed
for(j in 1:length(X)){
  mean_curve = colMeans(X_smoothed[[j]])
  for(i in 1:dim(X[[1]])[1]){
    X_smoothed_central[[j]][i,] <- X_smoothed[[j]][i,] - mean_curve
  }
}



result_lg = grplFlogit(Y=y, X=X_smoothed_central, Tps=tps, lambda=20, phi=.1)
selected_features_lg <- c()
for(i in 1:length(result_lg$Coef)){
  
  if(!all(result_lg$Coef[[i]] == 0)){
    selected_features_lg <- c(selected_features_lg, i)
  }
  
}
load("UCI_data.RData")
names(X)[selected_features_lg]
