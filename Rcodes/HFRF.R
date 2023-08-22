

gen.minkowski.dist <- function(x1, x2, range, p){  sum((abs(x1 - x2)/range)^p)^(1/p)  }

HFRF <- function(formula, m, c = NULL, p, train, test, train_x, test_x, train_y, test_y, 
                 train_control, contVars = NULL, ranges = NULL, genDist = NULL, svmHypers = NULL){
  
  # Clustering 
  nomVars <- setdiff(1:ncol(train_x),contVars)
  range <- ranges[-3]

  contVarsD <- c(contVars,contVars[length(contVars)]+1)
  nomVarsD <- setdiff(1:ncol(train),contVarsD)
  trainScaled <- cbind(scale(train[,contVarsD]),train[,nomVarsD])
    
  if (is.null(c)){
      kumeList <- list()
      hierAglo <- list()
      maxClustNum <- 6 
      SIL <- array(NA,(maxClustNum-1))
      for (nClust in 2:maxClustNum){
        kumeList[[nClust-1]] <- agnes(genDist, diss = TRUE, method = "ward")
        hierAglo[[nClust-1]] <- cutree(kumeList[[nClust-1]], k = nClust)
        ss <- silhouette(hierAglo[[nClust-1]], genDist)
        SIL[nClust-1] <- mean(ss[,3])
        # print(SIL)
      }
      kume <- hierAglo[[which.max(SIL)]]
      center <- aggregate(train/ranges, list(kume=kume), mean)
      new_center <- center[,-c(1,4)]
      c <- nrow(center)
  } else {
      kume <- agnes(genDist, diss = TRUE,  method = "ward")
      kume <- cutree(kume, k = c)
      center <- aggregate(train/ranges, list(kume=kume), mean)
      new_center <- center[,-c(1,4)]
      # c <- nrow(center)
  }
  
  
  # Prepare training and test data 
  dist <- data.frame(matrix(nrow=nrow(train_x), ncol=c))
  mu_train <- data.frame(matrix(nrow=nrow(train_x), ncol=c))
  genMinkowskiDist <- matrix(NA,nrow = nrow(train_x), ncol = c)
  for (j in 1:c){
    for (i in 1:nrow(train_x)){
      genMinkowskiDist[i,j] <-  gen.minkowski.dist(as.numeric(train_x[i,]), as.numeric(new_center[j,]), 
                                             range = range, p = p) # euc.dist_train(as.numeric(train_x[i,contVars]), as.numeric(new_center[j,contVars]))
      
    }
  }
  
  for (i in 1:c){
    dist[,i] <- genMinkowskiDist[,i] 
  }
  
  for (i in 1:nrow(dist)){
    for (j in 1:c){
      tvalue=0
      for (k in 1:c){
        value <- (dist[i,j]/dist[i,k])^(2/(m-1))
        tvalue=tvalue+value
      }
      mu_train[i,j]=1/tvalue
    }
  }
  
  Utrain_kume <- list()
  trainNew <- list()
  for ( i in 1:c){
    Utrain_kume[[i]] <- data.frame(m1 = (mu_train[,i]))
    trainNew[[i]] <- cbind(train_xy, Utrain_kume[[i]])
  }
  
  dist <- data.frame(matrix(nrow=nrow(test_x), ncol=c))
  mu_test <- data.frame(matrix(nrow=nrow(test_x), ncol=c))
  
  genMinkowskiDist <- matrix(NA,nrow = nrow(test_x), ncol = c)
  for (j in 1:c){
    for (i in 1:nrow(test_x)){
      genMinkowskiDist[i,j] <-  gen.minkowski.dist(as.numeric(test_x[i,]), as.numeric(new_center[j,]), 
                                             range = range, p = p) 
    }
  }
  
  for (i in 1:c){
    dist[,i] <- genMinkowskiDist[,i] 
  }
  
  for (i in 1:nrow(dist)){
    for (j in 1:c){
      tvalue=0
      for (k in 1:c){
        value <- (dist[i,j]/dist[i,k])^(2/(m-1))
        tvalue=tvalue+value
      }
      mu_test[i,j]=1/tvalue
    }
  }
  
  Utest_kume <- list()
  testNew <- list()
  for ( i in 1:c){
    Utest_kume[[i]] <- data.frame(m1 = (mu_test[,i])) 
    testNew[[i]] <- cbind(test_x, Utest_kume[[i]])
  }
  
  # Fit the model
  HFRF_model <- list()
  fitted.results <- list()
  
  for ( i in 1:c){
    if (is.null(svmHypers)){
      HFRF_model[[i]] <- train(formula, data = trainNew[[i]],  trControl = trainControl(method = "cv", search = "random"), method = "svmRadial",
                              preProcess = c("center", "scale"), seed=12342) #C and sigma are optimised automatically wth RBF kernel
     } else {
      HFRF_model[[i]] <- train(formula, data = trainNew[[i]],  trControl = trainControl(method = "cv"), method = "svmRadial",
                              preProcess = c("center", "scale"),tuneGrid = data.frame(C = svmHypers[i,2],sigma = svmHypers[i,1]), seed=12342)
     }
  
    fitted.results[[i]] <- predict(HFRF_model[[i]], newdata=trainNew[[i]])
    if ( i == 1){
      fit <- fitted.results[[i]]
    } else {
      fit <- cbind(fit, fitted.results[[i]])
    }
  }
  fit <- fit*mu_train
  fitted.results_train <- rowSums(fit)
  RMSE_train <- (sum((fitted.results_train - train_y)^2)/length(train_y))^0.5
  R_squared_train <- 1 - (sum((fitted.results_train - train_y)^2) / sum((train_y - mean(train_y))^2))
  MAE_train <- mean(abs(fitted.results_train - train_y))
  rMAE_train <- MAE_train / mean(train_y)
  rRMSE_train <- RMSE_train / mean(train_y)
  GOF_train <- c(RMSE_train, rRMSE_train, MAE_train, rMAE_train, R_squared_train)
  
  fitted.results <- list()
  for ( i in 1:c){
    fitted.results[[i]] <- predict(HFRF_model[[i]], newdata=testNew[[i]])
    if ( i == 1){
      fit <- fitted.results[[i]]
    } else {
      fit <- cbind(fit, fitted.results[[i]])
    }
  }
  
  fit_test <- fit*mu_test
  fitted.results <- rowSums(fit_test)
  
  RMSE_test <- (sum((fitted.results - test_y)^2)/length(test_y))^0.5
  MAE_test <- mean(abs(fitted.results - test_y))
  rMAE_test <- MAE_test / mean(test_y)
  rRMSE_test <- RMSE_test / mean(test_y)
  GOF_test <- c(RMSE_test, rRMSE_test, MAE_test, rMAE_test)
  
  return(list(GOF_train = GOF_train, GOF_test = GOF_test, fits.train = fitted.results_train, nClust = c))
}


