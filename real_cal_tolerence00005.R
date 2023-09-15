set.seed(123456)
load("UCI_data.RData")

require(fda)
require(LiblineaR)
require(parallel)
require(grplasso)

train_test_split <- function(X,y,train_size=n,test_size=test_size, mask=mask){
  
  p <- length(X)
  
  #train_ind <- sample(which(y==1), size=train_size/2)
  #train_ind <- c(train_ind, sample(which(y==-1), size=train_size/2))
  train_ind <- sample(seq(1,length(y), by=1), size=train_size)
  #print(mask)
  cond1 <- y==1
  cond2 <- y==-1
  cond3 <- !(seq(1,length(y), by=1)%in% train_ind)
  #cond4 <- !(seq(1,length(y), by=1)%in% mask)
  
  #test_ind <- sample(which(cond1&cond3), size=test_size/2)
  #test_ind <- c(test_ind, sample(which(cond2&cond3), size=test_size/2))
  test_ind <- sample(which(cond3), size=test_size)
  
  y_train <- y[train_ind]
  y_test <- y[test_ind]
  
  X_train <- vector(mode = "list", length=p)
  X_test <- vector(mode = "list", length=p)
  
  for(j in 1:p){
    X_train[[j]] <- X[[j]][train_ind,]
    X_test[[j]] <- X[[j]][test_ind,]
  } 
  result <- list("X_train"=X_train, "X_test"=X_test,
                 "y_train"=y_train,"y_test"=y_test)
  return(result)
}



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
        
        
        if(max(colSums(abs(C_new-C))) <= tolerance){
          
          break
        }
        if(iter_gd >= maxiter){
          
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
    
    
    convergence_c1 <- abs(Alpha0_new-Alpha0) <= tolerance
    convergence_c2 <- max(abs(Alpha_new - Alpha)) <= tolerance
    convergence_c3 <- max(colSums(abs(C_new-C))) <= tolerance
    # max(abs( c(Alpha0_new-Alpha0, Alpha_new-Alpha, C_new-C) )) <= tolerance
    if(convergence_c1 & convergence_c2 & convergence_c3){  
      
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






predict_l1fsvm <- function(X, tps, trained_model){
  
  Alpha0_new <- trained_model$Alpha0
  Alpha_new <- trained_model$Alpha
  C_new <- trained_model$C
  
  q <- dim(C_new)[1]
  p <- dim(C_new)[2]
  Basis <- vector(mode = "list", length = p)
  Z <- vector(mode = "list", length = p)
  
  compute_Zj <- function(j, tps,nbasis, X, Basis){
    
    max_tps <- max(tps[[j]])
    min_tps <- min(tps[[j]])
    length_tps <- length(tps[[j]])
    bBasis_obj <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis)
    
    if(min_tps != 0){
      # transfer tps to [0,1]
      tps_trans <- c(0, (tps[[j]]/tps[[j]][length_tps-1])[-length_tps]) 
      delta_T <- tps_trans[2] - tps_trans[1]
      # T by nbasis(q) matrix 
      Basis[[j]] <- eval.basis(tps_trans, bBasis_obj) 
    }else{
      tps_trans <- tps[[j]]/max_tps # transfer tps to [0,1]
      delta_T <- tps_trans[2] - tps_trans[1]
      Basis[[j]] <- eval.basis(tps_trans, bBasis_obj)
    }
    
    X[[j]]%*%Basis[[j]]*delta_T 
  }
  
  
  Z <- lapply(1:p, FUN=function(j){compute_Zj(j, tps, nbasis=q, X, Basis)})
  
  ZC <- sapply(1:p, FUN = function(j) Z[[j]]%*%C_new[,j,drop=FALSE]) #n by p matrix
  f <- (Alpha0_new + ZC%*%Alpha_new)
  sign(f)
  
}





cv_l1fsvm <- function(para_ind,test_cases,n,p,tps,nfold,folds,...){
  
  cost <- test_cases[para_ind,1]
  nbasis <- test_cases[para_ind,2]
  scores <- rep(NA,nfold) 
  
  ############################ evaluate the model ##################################
  
  cv_eva <- function(i, l1fsvm, n, p, folds, nfold, tps, cost, nbasis,...){
    # for cv validation set:
    X_valid <- folds[[i]]$X # a list of p, each ele is n by T 
    y_valid <- folds[[i]]$y # a vector of n
    
    
    
    # for cv training set:
    X_training <- vector(mode="list", length=p)
    n_cvtrain <- n - length(y_valid)
    y_training <- rep(NA,n_cvtrain)
    
    ender = 0
    for(f in 1:(nfold-1)){ 
      X_training <- lapply(1:p, FUN=function(j){
        
        
        rbind(X_training[[j]], folds[-i][[f]]$X[[j]])})
      
      starter = 1 + ender
      ender = ender + length(folds[-i][[f]]$y)
      
      y_training[starter:ender] <- folds[-i][[f]]$y
    }
  
    
    # evaluation:
    m <- l1fsvm(X_training, y_training, tps, cost=cost, nbasis=nbasis, ...)
    sum(predict_l1fsvm(X_valid, tps, trained_model=m) == y_valid)/length(y_valid)
  }
  
  
  scores <- sapply(1:nfold, FUN=function(i){
    cv_eva(i, l1fsvm, n, p, folds, nfold, tps, cost, nbasis,...)})
  
  cat(paste0("cost: ",cost,", nbasis: ",nbasis),"\n")
  cat("Prediction accuracy in k-folds: ", scores,"\n")
  return(c(mean(scores), sd(scores)/sqrt(nfold)))
}



gridsearch_cv_l1fsvm <- function(X,y,tps,nfold=nfold,cost_paras=c(1,5),
                          nbasis_paras=c(15,20),OneStandardError=TRUE,...){
  
  ############################ Split the data ################################
  n <- length(y)
  
  if(nfold<2 | nfold>n){stop("k is not in an appropriate range.\n")}
  
  p <- length(X)
  
  fold_size <- floor(n/nfold)
  folds <- vector(mode="list", length=nfold) 
  
  ind_pool <- 1:n
  select_ind <- vector(mode="list", length=nfold)
  
  for(i in 1:nfold){
    
    if(i == nfold){
      X_select <- lapply(1:p, FUN=function(j){X[[j]][ind_pool,,drop=FALSE]})
      y_select <- y[ind_pool]
      
      folds[[i]] <- list("X"=X_select, "y"=y_select)
      break
    }
    
    select_ind[[i]] <- sample(ind_pool, fold_size)
    ind_pool <- ind_pool[!ind_pool %in% select_ind[[i]]]
    
    X_select <- vector(mode="list", length=p)
    X_select <- lapply(1:p, FUN=function(j){X[[j]][select_ind[[i]],,drop=FALSE]})
    
    y_select <- y[select_ind[[i]]]
    
    folds[[i]] <- list("X"=X_select, "y"=y_select)
    
  }
  
  ####################### Evaluation each test case #########################
  
  test_cases <- as.matrix(expand.grid(cost_paras=sort(cost_paras),
                                      nbasis_paras=sort(nbasis_paras)))
  n_test_cases <- nrow(test_cases)
  accuracy_se <- vector(mode="list", length=n_test_cases)
  
  if(Sys.info()['sysname'] == "Windows"){
    accuracy_se <- lapply(1:nrow(test_cases),
                          FUN =function(para_ind){
                            cv_l1fsvm(para_ind,test_cases,n=n,p=p,tps=tps,
                               nfold=nfold,folds=folds,...)})
  }else{
    ncores <- detectCores()
    if(ncores >= n_test_cases){
      accuracy_se <- mclapply(1:nrow(test_cases),
                              FUN =function(para_ind){
                                cv_l1fsvm(para_ind,test_cases,n=n,p=p,tps=tps,
                                   nfold=nfold,folds=folds,...)},
                              mc.cores=n_test_cases)
    }else{
      accuracy_se <- mclapply(1:nrow(test_cases),
                              FUN =function(para_ind){
                                cv_l1fsvm(para_ind,test_cases,n=n,p=p,tps=tps,
                                   nfold=nfold,folds=folds,...)},
                              mc.cores=ncores)
    }
  }
  
  accuracy_se <- matrix(unlist(accuracy_se), nrow=2)
  
  if(OneStandardError==TRUE){
    accuracy <- accuracy_se[1,]
    
    standard_error <- accuracy_se[2,]
    
    within_se <- accuracy >= (max(accuracy) - standard_error[which.max(accuracy)])
    
    # the nbasis will tend to be the smaller one since sorted
    best_case <- which.min(test_cases[within_se, 1]) 
    
    return(list("accuracy"=accuracy,
                "accuracy_standard_error"=standard_error,
                "test_cases"=test_cases,
                "best_para"=test_cases[within_se,,drop=FALSE][best_case,]))
  }else{
    accuracy <- accuracy_se[1,]
    return(list("accuracy"=accuracy,
                "test_cases"=test_cases,
                "best_para"=test_cases[which.max(accuracy),]))
  }
}



#####################################################################################
#####################    lg


# grplFlogit author: J. Gertheiss  A. Maity  A.a Staicu (2013)
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



cv_lg <- function(para_ind,test_cases,n,p,tps,nfold,folds,...){
  
  lambda <- test_cases[para_ind,1]
  phi <- test_cases[para_ind,2]
  scores <- rep(NA,nfold) 
  
  ############################ evaluate the model ##################################
  
  cv_eva <- function(i, grplFlogit, n, p, folds, nfold, tps, lambda, phi,...){
    # for cv validation set:
    X_valid <- folds[[i]]$X # a list of p, each ele is n by T 
    y_valid <- folds[[i]]$y # a vector of n
    
    
    # for cv training set:
    X_training <- vector(mode="list", length=p)
    n_cvtrain <- n - length(y_valid)
    y_training <- rep(NA,n_cvtrain)
    
    ender = 0
    for(f in 1:(nfold-1)){ 
      X_training <- lapply(1:p, FUN=function(j){
        rbind(X_training[[j]], folds[-i][[f]]$X[[j]])})
      
      starter = 1 + ender
      ender = ender + length(folds[-i][[f]]$y)
      
      y_training[starter:ender] <- folds[-i][[f]]$y
    }
    
    
    # evaluation:
    m <- grplFlogit(Y=y_training, X=X_training, Tps=tps, lambda=lambda, phi=phi)
    
    
    sum(predict_lg(X_valid=X_valid, Tps=tps, trained_model=m) == y_valid)/length(y_valid)
  }
  
  scores <- sapply(1:nfold, FUN=function(i){
    cv_eva(i,grplFlogit, n, p, folds, nfold, tps, lambda, phi,...)})
  
  cat(paste0("lambda: ",lambda,", phi: ",phi),"\n")
  cat("Prediction accuracy in k-folds: ", scores,"\n")
  return(mean(scores))
}



gridsearch_cv_lg <- function(X,y,tps,nfold,lambda_paras=10^seq(3,0,by=-1),
                             phi_paras=10^c(10,8,6,4,2),...){
  
  ############################ Split the data ################################
  n <- length(y)
  
  if(nfold<2 | nfold>n){stop("k is not in an appropriate range.\n")}
  
  p <- length(X)
  
  fold_size <- floor(n/nfold)
  folds <- vector(mode="list", length=nfold) 
  
  ind_pool <- 1:n
  select_ind <- vector(mode="list", length=nfold)
  
  for(i in 1:nfold){
    
    if(i == nfold){
      X_select <- lapply(1:p, FUN=function(j){X[[j]][ind_pool,,drop=FALSE]})
      y_select <- y[ind_pool]
      
      folds[[i]] <- list("X"=X_select, "y"=y_select)
      break
    }
    
    select_ind[[i]] <- sample(ind_pool, fold_size)
    ind_pool <- ind_pool[!ind_pool %in% select_ind[[i]]]
    
    X_select <- vector(mode="list", length=p)
    X_select <- lapply(1:p, FUN=function(j){X[[j]][select_ind[[i]],,drop=FALSE]})
    
    y_select <- y[select_ind[[i]]]
    
    folds[[i]] <- list("X"=X_select, "y"=y_select)
    
  }
  
  ####################### Evaluation each test case #########################
  
  test_cases <- as.matrix(expand.grid(lambda_paras=sort(lambda_paras),
                                      phi_paras=sort(phi_paras)))
  n_test_cases <- nrow(test_cases)
  
  if(Sys.info()['sysname'] == "Windows"){
    accuracy <- lapply(1:nrow(test_cases),
                       FUN =function(para_ind){
                         cv_lg(para_ind,test_cases,n=n,p=p,tps=tps,
                               nfold=nfold,folds=folds,...)})
  }else{
    ncores <- detectCores()
    if(ncores >= n_test_cases){
      accuracy <- mclapply(1:nrow(test_cases),
                           FUN =function(para_ind){
                             cv_lg(para_ind,test_cases,n=n,p=p,tps=tps,
                                   nfold=nfold,folds=folds,...)},
                           mc.cores=n_test_cases)
    }else{
      accuracy <- mclapply(1:nrow(test_cases),
                           FUN =function(para_ind){
                             cv_lg(para_ind,test_cases,n=n,p=p,tps=tps,
                                   nfold=nfold,folds=folds,...)},
                           mc.cores=ncores)
    }
  }
  
  return(list("accuracy"=accuracy,
              "test_cases"=test_cases,
              "best_para"=test_cases[which.max(accuracy),]))
  
}




##################################################################################
##############################  run  #############################################


one_run <- function(rep, X,y,tps,nfold, tolerance=0.001, maxiter=100000,
                    cost_paras, 
                    nbasis_paras,lambda_paras, phi_paras,train_size=98, test_size=24,...){
  cat("rep == ", rep, "\n")
  split_result <- train_test_split(X=X,y=y,train_size=train_size, test_size=test_size)
  X_train <- split_result$X_train 
  X_test <- split_result$X_test
  y_train <- split_result$y_train
  y_test <- split_result$y_test
  
  
  #   l1fsvm
  r <- gridsearch_cv_l1fsvm(X=X_train, y=y_train, tps=tps, nfold=nfold,
                     cost_paras=cost_paras,
                     nbasis_paras=nbasis_paras,
                     tolerance=tolerance, maxiter=maxiter,...)
  
  best_cost <- r$best_para[1]
  best_nbasis <- r$best_para[2]
  
  m <- l1fsvm(X=X_train, y=y_train, tps=tps, cost=best_cost, nbasis=best_nbasis,
              tolerance=tolerance)
  
  test_accuracy <- sum(predict_l1fsvm(X=X_test, tps=tps, trained_model=m) == y_test)/length(y_test)
  
  result <- c("test_accuracy_fsvm"=test_accuracy,
              "best_cost_fsvm"=best_cost,
              "best_nbasis_fsvm"=best_nbasis)
  
  #grpllg
  y_train[y_train==-1] <- 0
  y_test[y_test==-1] <- 0
  r <- gridsearch_cv_lg(X=X_train, y=y_train, tps=tps, nfold=nfold,
                        lambda_paras=lambda_paras,
                        phi_paras=phi_paras,...)
  
  best_lambda <- r$best_para[1]
  best_phi <- r$best_para[2]
  
  m <- grplFlogit(X=X_train, Y=y_train, Tps=tps, lambda=best_lambda, phi=best_phi)
  
  test_accuracy <- sum(predict_lg(X=X_test, Tps=tps, trained_model=m) == y_test)/length(y_test)
  
  result <- c(result, c("test_accuracy_lg"=test_accuracy,
                        "best_lambda_lg"=best_lambda,
                        "best_phi_lg"=best_phi))
  return(result)
  
}

rep_run <- function(X, y, tps,
                    njobs=1, nrep=100,nfold=2,tolerance=0.0001,
                    maxiter=1000000,...){
  
  rep_run_result <- mclapply(1:nrep,
                             FUN=function(rep){
                               one_run(rep, X, y, tps, tolerance=tolerance, maxiter=maxiter,
                                      cost_paras=c(0.05,0.07,0.09, 0.12, 0.15, 0.22),  
                                      nbasis_paras=c(5, 10, 15),nfold=nfold,
                                      lambda_paras=c(0.1,1,10,15,20),   
                                      phi_paras=c(0.1,1,10),...)},
                             mc.cores=njobs)
  
  rep_run_result <- matrix(unlist(rep_run_result),nrow=6)
  rownames(rep_run_result) <- c("test_accuracy_fsvm","best_cost_fsvm","best_nbasis_fsvm",
                                "test_accuracy_lg","best_lambda_lg","best_phi_lg")
  #saveRDS(rep_run_result, file=paste0("real_rep_run_result"))
  return(rep_run_result)
}


y[y==0] <- -1
#pos_ind <- sample(which(y==1), min(sum(y==1), sum(y==-1)))
#neg_ind <- sample(which(y==-1), min(sum(y==1), sum(y==-1)))
#y <- y[c(pos_ind, neg_ind)]

#X <- lapply(1:length(X), function(i){X[[i]][c(pos_ind, neg_ind),]})



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

n_cores <- detectCores(logical=FALSE)
set.seed(1)
realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
                                   njobs=10, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=100000)
saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_1.rds")
set.seed(2)
realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
                                  njobs=10, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=100000)
saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_2.rds")
set.seed(3)
realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
                                  njobs=10, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=100000)
saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_3.rds")
set.seed(4)
realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
                                  njobs=10, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=100000)
saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_4.rds")
# set.seed(5)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_5.rds")
# set.seed(6)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_6.rds")
# set.seed(7)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_7.rds")
# set.seed(8)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_8.rds")
# set.seed(9)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_9.rds")
# set.seed(10)
# realdata_calculation_fsvm= rep_run(X=X_smoothed_central,y=y,tps=tps,
#                                   njobs=n_cores, nrep=10,OneStandardError=TRUE,nfold = 5, lr=0.01, tolerance=0.0005, maxiter=1000000)
# saveRDS(realdata_calculation_fsvm, file="realdata_calculation_fsvm_and_lg_tl00005_10.rds")
# 
# 




